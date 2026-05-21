"""ROC benchmark for cablastp only (fallback when hs-cablastp is unavailable)."""
from __future__ import annotations

import argparse
import csv
import os
import random
import subprocess
import sys
import time
from pathlib import Path
from typing import Dict, Iterable, List, Set, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

Pair = Tuple[str, str]


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


def run_blastp_tabular(blastp: str, query: Path, db_prefix: Path, evalue: float) -> List[List[str]]:
    proc = subprocess.run(
        [blastp, "-query", str(query), "-db", str(db_prefix),
         "-evalue", str(evalue), "-outfmt", f"6 {_TAB_COLS}"],
        capture_output=True, text=True,
    )
    if proc.returncode != 0:
        raise RuntimeError(f"blastp failed:\n{proc.stderr}")
    return [ln.split("\t") for ln in proc.stdout.splitlines() if ln.strip()]


def parse_tabular_rows(rows: Iterable[List[str]]) -> Dict[Pair, float]:
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


def run_cablastp_pipeline(query_fasta: Path, db_dir: Path, evalue: float) -> Dict[Pair, float]:
    cmd = ["cablastp-search", str(db_dir), str(query_fasta), "--quiet",
           "--blast-args", "-outfmt", "6", "-evalue", str(evalue)]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"cablastp-search failed (exit {proc.returncode}):\n{proc.stderr}")
    rows = [ln.split("\t") for ln in proc.stdout.splitlines() if ln.strip()]
    return parse_tabular_rows(rows)


def compute_roc(positives: Set[Pair], universe: Set[Pair], scores: Dict[Pair, float]):
    P = len(positives)
    N = len(universe) - P
    if P == 0 or N == 0:
        return [0.0, 1.0], [0.0, 1.0], float("nan")
    pairs = sorted(universe, key=lambda p: scores.get(p, 0.0), reverse=True)
    fpr_curve = [0.0]
    tpr_curve = [0.0]
    tp = fp = 0
    prev_score = None
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
    auc = sum(
        (fpr_curve[i] - fpr_curve[i - 1]) * (tpr_curve[i] + tpr_curve[i - 1]) / 2.0
        for i in range(1, len(fpr_curve))
    )
    return fpr_curve, tpr_curve, auc


def _dir_size_kb(path: Path) -> float:
    total = 0
    for root, _dirs, files in os.walk(path):
        for f in files:
            try:
                total += (Path(root) / f).stat().st_size
            except OSError:
                pass
    return total / 1024.0


def _append_csv(path: Path, rows: List[Dict]) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    write_header = not path.exists() or path.stat().st_size == 0
    with open(path, "a", newline="", encoding="ascii") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        if write_header:
            writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", type=Path, required=True)
    ap.add_argument("--cab-db", type=Path, required=True)
    ap.add_argument("--n-queries", type=int, default=50)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--ground-truth-evalue", type=float, default=1e-3)
    ap.add_argument("--pipeline-evalue", type=float, default=1e-3)
    ap.add_argument("--blastp", default="blastp")
    ap.add_argument("--makeblastdb", default="makeblastdb")
    ap.add_argument("--out", type=Path, default=Path("bench"))
    ap.add_argument("--summary-csv", type=Path, default=None)
    ap.add_argument("--roc-csv", type=Path, default=None)
    ap.add_argument("--cab-compress-seconds", type=float, default=None)
    ap.add_argument("--label", default=None)
    args = ap.parse_args()

    args.out.mkdir(parents=True, exist_ok=True)

    records = parse_fasta(args.fasta)
    rng = random.Random(args.seed)
    queries = rng.sample(records, args.n_queries)
    query_ids: Set[str] = {q[0] for q in queries}
    subjects = [r[0] for r in records]
    subject_set: Set[str] = set(subjects)

    queries_path = args.out / "queries.fasta"
    write_fasta(queries_path, queries)
    print(f"[1/4] sampled {len(queries)} queries -> {queries_path}")

    raw_db = args.out / "raw_blast_db"
    ensure_raw_blastdb(args.fasta, raw_db, args.makeblastdb)
    t = time.perf_counter()
    gt_rows = run_blastp_tabular(args.blastp, queries_path, raw_db, args.ground_truth_evalue)
    gt_scores = parse_tabular_rows(gt_rows)
    gt_positives: Set[Pair] = set(gt_scores.keys())
    print(f"[2/4] vanilla blastp: {len(gt_rows)} rows, {len(gt_positives)} unique positives "
          f"({time.perf_counter() - t:.1f}s)")

    t = time.perf_counter()
    cab_scores = run_cablastp_pipeline(queries_path, args.cab_db, args.pipeline_evalue)
    cab_search_seconds = time.perf_counter() - t
    print(f"[3/4] cablastp: {len(cab_scores)} unique hit pairs ({cab_search_seconds:.1f}s)")

    universe: Set[Pair] = {(q, s) for q in query_ids for s in subject_set}
    cab_missed = sorted(gt_positives - set(cab_scores.keys()))
    cab_extra = sorted(set(cab_scores.keys()) - gt_positives)

    def dump(path: Path, label: str, pairs, scores=None):
        with open(path, "w", encoding="ascii") as fh:
            fh.write(f"# {label}: {len(pairs)} pairs\n")
            fh.write("# qseqid\tsseqid\tscore\n")
            for q, s in pairs:
                bs = scores.get((q, s)) if scores else None
                fh.write(f"{q}\t{s}\t{bs if bs is not None else ''}\n")

    dump(args.out / "missed_cablastp.tsv", "ground-truth pairs missed by cablastp", cab_missed, gt_scores)
    dump(args.out / "extra_cablastp.tsv", "cablastp hits not in ground truth", cab_extra, cab_scores)

    cab_fpr, cab_tpr, cab_auc = compute_roc(gt_positives, universe, cab_scores)
    cab_tp = len(set(cab_scores) & gt_positives)

    # Plot (cablastp only)
    fig, (ax_full, ax_zoom) = plt.subplots(1, 2, figsize=(12, 5.5))
    for ax, title, xlim, ylim in (
        (ax_full, "ROC (full)", (0, 1), (0, 1)),
        (ax_zoom, "ROC (low FPR zoom)", (0, 0.05), (0, 1)),
    ):
        ax.plot(cab_fpr, cab_tpr, label=f"cablastp  (AUC={cab_auc:.4f})", linewidth=1.6)
        ax.plot([0, 1], [0, 1], color="gray", linestyle=":", linewidth=1, label="random")
        ax.set_xlabel("False Positive Rate")
        ax.set_ylabel("True Positive Rate")
        ax.set_xlim(xlim); ax.set_ylim(ylim)
        ax.set_title(title)
        ax.legend(loc="lower right")
        ax.grid(True, alpha=0.3)
    fig.suptitle(f"ROC (cablastp only, hs-cablastp TIMED OUT): {args.fasta.name} "
                 f"({args.n_queries} queries)", fontsize=11)
    fig.tight_layout()
    roc_png = args.out / "roc.png"
    fig.savefig(roc_png, dpi=130)
    print(f"[4/4] saved ROC plot to {roc_png}")

    print()
    print("=" * 64)
    print(f"cablastp: {len(cab_scores)} hits, {cab_tp} TP, {len(cab_missed)} missed, "
          f"recall={cab_tp / max(len(gt_positives), 1):.4f}, AUC={cab_auc:.4f}")
    print("NOTE: hs-cablastp TIMED OUT — no hs results")
    print("=" * 64)

    label = args.label or args.fasta.name
    cab_db_kb = _dir_size_kb(args.cab_db)

    summary_rows = [
        {
            "dataset": label,
            "pipeline": "cablastp",
            "n_queries": args.n_queries,
            "n_subjects": len(subjects),
            "universe_pairs": len(universe),
            "ground_truth_positives": len(gt_positives),
            "compress_seconds": args.cab_compress_seconds if args.cab_compress_seconds is not None else "",
            "db_size_kb": f"{cab_db_kb:.2f}",
            "search_seconds": f"{cab_search_seconds:.3f}",
            "total_hit_pairs": len(cab_scores),
            "true_positives": cab_tp,
            "missed_pairs": len(cab_missed),
            "extra_pairs": len(cab_extra),
            "recall": f"{cab_tp / max(len(gt_positives), 1):.6f}",
            "precision": f"{cab_tp / max(len(cab_scores), 1):.6f}",
            "auc": f"{cab_auc:.6f}",
        },
        {
            "dataset": label,
            "pipeline": "hs-cablastp",
            "n_queries": "TIMEOUT",
            "n_subjects": "TIMEOUT",
            "universe_pairs": "TIMEOUT",
            "ground_truth_positives": "TIMEOUT",
            "compress_seconds": "TIMEOUT",
            "db_size_kb": "TIMEOUT",
            "search_seconds": "TIMEOUT",
            "total_hit_pairs": "TIMEOUT",
            "true_positives": "TIMEOUT",
            "missed_pairs": "TIMEOUT",
            "extra_pairs": "TIMEOUT",
            "recall": "TIMEOUT",
            "precision": "TIMEOUT",
            "auc": "TIMEOUT",
        },
    ]
    if args.summary_csv:
        _append_csv(args.summary_csv, summary_rows)
        print(f"appended 2 rows to {args.summary_csv}")

    if args.roc_csv:
        roc_rows = []
        for fpr, tpr in zip(cab_fpr, cab_tpr):
            roc_rows.append({"dataset": label, "pipeline": "cablastp",
                             "fpr": f"{fpr:.6f}", "tpr": f"{tpr:.6f}"})
        _append_csv(args.roc_csv, roc_rows)
        print(f"appended {len(roc_rows)} rows to {args.roc_csv}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
