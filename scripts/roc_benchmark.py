"""ROC benchmark for cablastp vs hs-cablastp.

Steps:
1. Sample N query sequences from an input FASTA.
2. Build ground-truth set: vanilla blastp of queries against the raw FASTA
   at evalue <= GROUND_TRUTH_EVALUE. Every surviving (query, subject) pair is
   a "true positive".
3. Run cablastp-search and hs-cablastp-search on the same batched query
   FASTA against their respective compressed DBs.
4. Normalize hits to (qseqid, subject-original-id) pairs with a bitscore.
5. Compute ROC over (query x dataset) pair universe for each pipeline; report
   AUC and the count of missed ground-truth pairs.
6. Plot ROC and save bench/roc.png.

Usage:
    python scripts/roc_benchmark.py \
        --fasta data/medium.fasta \
        --cab-db bench/cab_db \
        --hs-db  bench/hs_db \
        --n-queries 50 \
        --out bench

`bench/raw_blast_db.*` is built automatically if missing.
"""
from __future__ import annotations

import argparse
import os
import random
import shutil
import subprocess
import sys
import time
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Set, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from hs_cablastp.search import search as hs_search


Pair = Tuple[str, str]


# -------------------- FASTA --------------------


def parse_fasta(path: Path) -> List[Tuple[str, str]]:
    name = None
    chunks: List[str] = []
    out: List[Tuple[str, str]] = []
    with open(path, "r", encoding="ascii", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                if name is not None:
                    out.append((name, "".join(chunks)))
                name = line[1:].split()[0] if len(line) > 1 else ""
                chunks = []
            elif line:
                chunks.append(line)
    if name is not None:
        out.append((name, "".join(chunks)))
    return out


def write_fasta(path: Path, records: Iterable[Tuple[str, str]]) -> None:
    with open(path, "w", encoding="ascii") as fh:
        for name, seq in records:
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")


# -------------------- BLAST helpers --------------------


_TAB_COLS = (
    "qseqid sseqid pident length mismatch gapopen "
    "qstart qend sstart send evalue bitscore"
)


def ensure_raw_blastdb(fasta: Path, db_prefix: Path, makeblastdb: str) -> None:
    if (db_prefix.parent / (db_prefix.name + ".psq")).exists():
        return
    db_prefix.parent.mkdir(parents=True, exist_ok=True)
    proc = subprocess.run(
        [makeblastdb, "-dbtype", "prot", "-in", str(fasta), "-out", str(db_prefix)],
        capture_output=True, text=True,
    )
    if proc.returncode != 0:
        raise RuntimeError(f"makeblastdb failed:\n{proc.stderr}")


def run_blastp_tabular(
    blastp: str, query: Path, db_prefix: Path, evalue: float,
) -> List[List[str]]:
    proc = subprocess.run(
        [
            blastp, "-query", str(query), "-db", str(db_prefix),
            "-evalue", str(evalue), "-outfmt", f"6 {_TAB_COLS}",
        ],
        capture_output=True, text=True,
    )
    if proc.returncode != 0:
        raise RuntimeError(f"blastp failed:\n{proc.stderr}")
    return [ln.split("\t") for ln in proc.stdout.splitlines() if ln.strip()]


def parse_tabular_rows(rows: Iterable[List[str]]) -> Dict[Pair, float]:
    """Return {(q, s): best_bitscore} from tabular rows."""
    best: Dict[Pair, float] = {}
    for r in rows:
        if len(r) < 12:
            continue
        try:
            score = float(r[11])
        except ValueError:
            continue
        key = (r[0], r[1])
        if key not in best or score > best[key]:
            best[key] = score
    return best


# -------------------- Pipeline runners --------------------


def run_cablastp_pipeline(
    query_fasta: Path, db_dir: Path, evalue: float,
) -> Dict[Pair, float]:
    """cablastp-search is a CLI script. We call it and parse tabular output from stdout."""
    cmd = [
        "cablastp-search", str(db_dir), str(query_fasta), "--quiet",
        "--blast-args", "-outfmt", "6", "-evalue", str(evalue),
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(
            f"cablastp-search failed (exit {proc.returncode}):\n{proc.stderr}"
        )
    rows = [ln.split("\t") for ln in proc.stdout.splitlines() if ln.strip()]
    return parse_tabular_rows(rows)


def run_hs_cablastp_pipeline(
    query_fasta: Path, db_dir: Path, fine_evalue: float,
) -> Dict[Pair, float]:
    """Call hs_cablastp.search.search() directly; map node hits to original ids."""
    result = hs_search(query_fasta, db_dir, fine_evalue=fine_evalue)
    best: Dict[Pair, float] = {}
    for h in result["fine_hits"]:
        # fasta_ref looks like "<orig_id>:<start>-<end>" -- keep only orig_id.
        # Strip everything from the first colon (yeast/UniProt ids don't contain ':').
        ref = h.fasta_ref.split(":", 1)[0].split()[0] if h.fasta_ref else ""
        if not ref:
            continue
        key = (h.qseqid, ref)
        if key not in best or h.bitscore > best[key]:
            best[key] = h.bitscore
    return best


# -------------------- ROC --------------------


def compute_roc(
    positives: Set[Pair], universe: Set[Pair], scores: Dict[Pair, float],
) -> Tuple[List[float], List[float], float]:
    """Return (fpr_curve, tpr_curve, auc). Pairs missing from scores have score=0."""
    P = len(positives)
    N = len(universe) - P
    if P == 0 or N == 0:
        return [0.0, 1.0], [0.0, 1.0], float("nan")

    pairs = sorted(universe, key=lambda p: scores.get(p, 0.0), reverse=True)
    fpr_curve = [0.0]
    tpr_curve = [0.0]
    tp = fp = 0
    prev_score: float | None = None
    for p in pairs:
        s = scores.get(p, 0.0)
        if prev_score is not None and s != prev_score:
            fpr_curve.append(fp / N)
            tpr_curve.append(tp / P)
        if p in positives:
            tp += 1
        else:
            fp += 1
        prev_score = s
    fpr_curve.append(fp / N)
    tpr_curve.append(tp / P)

    # Trapezoidal AUC.
    auc = 0.0
    for i in range(1, len(fpr_curve)):
        auc += (fpr_curve[i] - fpr_curve[i - 1]) * (
            tpr_curve[i] + tpr_curve[i - 1]
        ) / 2.0
    return fpr_curve, tpr_curve, auc


# -------------------- Driver --------------------


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fasta", type=Path, required=True,
                    help="Source FASTA (also the search target).")
    ap.add_argument("--cab-db", type=Path, required=True,
                    help="cablastp database directory.")
    ap.add_argument("--hs-db", type=Path, required=True,
                    help="hs-cablastp database directory.")
    ap.add_argument("--n-queries", type=int, default=50)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--ground-truth-evalue", type=float, default=1e-3)
    ap.add_argument("--pipeline-evalue", type=float, default=1e-3,
                    help="E-value passed to cablastp / hs-cablastp fine search.")
    ap.add_argument("--blastp", default="blastp")
    ap.add_argument("--makeblastdb", default="makeblastdb")
    ap.add_argument("--out", type=Path, default=Path("bench"),
                    help="Output directory for the queries FASTA, raw blast db, "
                         "miss lists, ROC PNG, and summary.")
    args = ap.parse_args()

    args.out.mkdir(parents=True, exist_ok=True)

    # 1) Sample queries.
    records = parse_fasta(args.fasta)
    if len(records) < args.n_queries:
        raise SystemExit(f"only {len(records)} records in {args.fasta}")
    rng = random.Random(args.seed)
    queries = rng.sample(records, args.n_queries)
    query_ids: Set[str] = {q[0] for q in queries}
    subjects: List[str] = [r[0] for r in records]
    subject_set: Set[str] = set(subjects)

    queries_path = args.out / "queries.fasta"
    write_fasta(queries_path, queries)
    print(f"[1/5] sampled {len(queries)} queries -> {queries_path}")

    # 2) Ground truth: vanilla blastp.
    raw_db = args.out / "raw_blast_db"
    ensure_raw_blastdb(args.fasta, raw_db, args.makeblastdb)
    t = time.perf_counter()
    gt_rows = run_blastp_tabular(args.blastp, queries_path, raw_db,
                                 args.ground_truth_evalue)
    gt_scores = parse_tabular_rows(gt_rows)
    gt_positives: Set[Pair] = set(gt_scores.keys())
    print(f"[2/5] vanilla blastp: {len(gt_rows)} rows, "
          f"{len(gt_positives)} unique positive pairs "
          f"({time.perf_counter() - t:.1f}s)")

    # 3) Run cablastp.
    t = time.perf_counter()
    cab_scores = run_cablastp_pipeline(
        queries_path, args.cab_db, args.pipeline_evalue,
    )
    print(f"[3/5] cablastp: {len(cab_scores)} unique hit pairs "
          f"({time.perf_counter() - t:.1f}s)")

    # 4) Run hs-cablastp.
    t = time.perf_counter()
    hs_scores = run_hs_cablastp_pipeline(
        queries_path, args.hs_db, args.pipeline_evalue,
    )
    print(f"[4/5] hs-cablastp: {len(hs_scores)} unique hit pairs "
          f"({time.perf_counter() - t:.1f}s)")

    # Constrain universe to (sampled queries) x (all subjects in fasta).
    universe: Set[Pair] = {(q, s) for q in query_ids for s in subject_set}

    # 5) Missed-hit report + ROC.
    cab_missed = sorted(gt_positives - set(cab_scores.keys()))
    hs_missed = sorted(gt_positives - set(hs_scores.keys()))
    cab_extra = sorted(set(cab_scores.keys()) - gt_positives)
    hs_extra = sorted(set(hs_scores.keys()) - gt_positives)

    def dump(path: Path, label: str, pairs: List[Pair],
             scores: Dict[Pair, float] | None = None) -> None:
        with open(path, "w", encoding="ascii") as fh:
            fh.write(f"# {label}: {len(pairs)} pairs\n")
            fh.write("# qseqid\tsseqid\tscore\n")
            for q, s in pairs:
                bs = scores.get((q, s)) if scores else None
                fh.write(f"{q}\t{s}\t{bs if bs is not None else ''}\n")

    dump(args.out / "missed_cablastp.tsv",
         "ground-truth pairs missed by cablastp", cab_missed, gt_scores)
    dump(args.out / "missed_hs_cablastp.tsv",
         "ground-truth pairs missed by hs-cablastp", hs_missed, gt_scores)
    dump(args.out / "extra_cablastp.tsv",
         "cablastp hits not in ground truth", cab_extra, cab_scores)
    dump(args.out / "extra_hs_cablastp.tsv",
         "hs-cablastp hits not in ground truth", hs_extra, hs_scores)

    cab_fpr, cab_tpr, cab_auc = compute_roc(gt_positives, universe, cab_scores)
    hs_fpr, hs_tpr, hs_auc = compute_roc(gt_positives, universe, hs_scores)

    # Plot.
    fig, (ax_full, ax_zoom) = plt.subplots(1, 2, figsize=(12, 5.5))
    for ax, title, xlim, ylim in (
        (ax_full, "ROC (full)", (0, 1), (0, 1)),
        (ax_zoom, "ROC (low FPR zoom)", (0, 0.05), (0, 1)),
    ):
        ax.plot(cab_fpr, cab_tpr, label=f"cablastp  (AUC={cab_auc:.4f})",
                linewidth=1.6)
        ax.plot(hs_fpr, hs_tpr, label=f"hs-cablastp (AUC={hs_auc:.4f})",
                linewidth=1.6)
        ax.plot([0, 1], [0, 1], color="gray", linestyle=":", linewidth=1,
                label="random")
        ax.set_xlabel("False Positive Rate")
        ax.set_ylabel("True Positive Rate")
        ax.set_xlim(xlim); ax.set_ylim(ylim)
        ax.set_title(title)
        ax.legend(loc="lower right")
        ax.grid(True, alpha=0.3)
    fig.suptitle(
        f"ROC: {args.fasta.name}  ({args.n_queries} queries, "
        f"GT evalue<={args.ground_truth_evalue}, pipeline evalue<="
        f"{args.pipeline_evalue})", fontsize=11,
    )
    fig.tight_layout()
    roc_png = args.out / "roc.png"
    fig.savefig(roc_png, dpi=130)
    print(f"[5/5] saved ROC plot to {roc_png}")

    # Summary
    print()
    print("=" * 64)
    print(f"Dataset: {args.fasta.name} ({len(subjects)} subjects)")
    print(f"Queries: {args.n_queries}  Seed: {args.seed}")
    print(f"Universe: {len(universe)} (query, subject) pairs")
    print(f"Ground-truth positives: {len(gt_positives)}")
    print("-" * 64)
    print(f"{'pipeline':<14}  {'hits':>6}  {'TP':>6}  {'missed':>7}  "
          f"{'extra':>6}  {'recall':>7}  {'AUC':>7}")
    cab_tp = len(set(cab_scores) & gt_positives)
    hs_tp = len(set(hs_scores) & gt_positives)
    print(f"{'cablastp':<14}  {len(cab_scores):>6d}  {cab_tp:>6d}  "
          f"{len(cab_missed):>7d}  {len(cab_extra):>6d}  "
          f"{cab_tp / max(len(gt_positives), 1):>7.4f}  {cab_auc:>7.4f}")
    print(f"{'hs-cablastp':<14}  {len(hs_scores):>6d}  {hs_tp:>6d}  "
          f"{len(hs_missed):>7d}  {len(hs_extra):>6d}  "
          f"{hs_tp / max(len(gt_positives), 1):>7.4f}  {hs_auc:>7.4f}")
    print("=" * 64)
    print(f"Missed-hit tables:")
    print(f"  {args.out / 'missed_cablastp.tsv'}")
    print(f"  {args.out / 'missed_hs_cablastp.tsv'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
