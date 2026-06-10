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
7. Optionally append rows to --summary-csv and --roc-csv for cross-dataset
   aggregation.

Usage:
    python scripts/roc_benchmark.py \
        --fasta data/ecoli_trembl_2k.fasta \
        --cab-db bench/etrembl2k_cab \
        --hs-db  bench/etrembl2k_hs \
        --n-queries 50 \
        --out bench \
        --summary-csv bench/summary.csv \
        --roc-csv bench/roc_points.csv

`bench/raw_blast_db.*` is built automatically if missing.

Compress timings and DB sizes are external inputs (compression isn't done by
this script). Pass them with --cab-compress-seconds / --hs-compress-seconds;
DB size is read from the disk footprint of --cab-db / --hs-db.
"""
from __future__ import annotations

import argparse
import csv
import os
import random
import shutil
import subprocess
import sys
import time
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Set, Tuple

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


def _extract_accession(seqid: str) -> str:
    """Reduce any seqid form we encounter to a bare protein accession.

    Handles:
      - UniProt-style `tr|ACC|NAME` / `sp|ACC|NAME` -> `ACC`
      - hs-cablastp refs like `tr|ACC|NAME:start-end` -> `ACC`
      - Plain accessions left untouched.

    Applied to *both* sides of the comparison (vanilla blastp ground truth
    and hs-cablastp / cablastp pipeline outputs) so a header-format quirk
    on one side can't desync the (query, subject) keys.
    """
    if not seqid:
        return ""
    # Drop hs-cablastp's `:start-end` coord suffix and any trailing whitespace.
    head = seqid.split(":", 1)[0].split()[0]
    # UniProt: db|acc|name -> acc
    if head.startswith(("tr|", "sp|")):
        parts = head.split("|")
        if len(parts) >= 2 and parts[1]:
            return parts[1]
    return head


def parse_tabular_rows(rows: Iterable[List[str]]) -> Dict[Pair, float]:
    """Return {(q_accession, s_accession): best_bitscore} from tabular rows."""
    best: Dict[Pair, float] = {}
    for r in rows:
        if len(r) < 12:
            continue
        try:
            score = float(r[11])
        except ValueError:
            continue
        key = (_extract_accession(r[0]), _extract_accession(r[1]))
        if key not in best or score > best[key]:
            best[key] = score
    return best


# -------------------- Pipeline runners --------------------


# Both pipelines run IN-PROCESS and are timed with the one-time DB load/open
# EXCLUDED, so the reported search time reflects algorithm cost rather than a
# storage-format artifact: hs-cablastp deserializes one pickle (~0.03s) while
# cablastp opens several index files (~0.5s) — a structural difference, not an
# algorithmic one. Interpreter startup / import are likewise excluded (both are
# Python; that cost is identical and not algorithmic). The two pipelines are
# also driven with the SAME coarse e-value (see --coarse-evalue) so neither
# casts a wider coarse net than the other.


def run_cablastp_pipeline(
    query_fasta: Path, db_dir: Path, fine_evalue: float, coarse_evalue: float,
) -> Tuple[float, Dict[Pair, float]]:
    """Run cablastp's search phases in-process; return (search_seconds, hits).

    new_read_db (the multi-file DB open) runs before the timer starts, so the
    timed region is the actual search work: coarse BLAST + hit decompression +
    fine DB build + fine BLAST. -dbsize is pinned to db.BlastDBSize, same as the
    cablastp-search CLI and as hs-cablastp / the ground truth.
    """
    import types
    from cablastp.db import new_read_db, FILE_BLAST_FINE
    from cablastp.commands import search as _cab
    from cablastp.commands._search_common import (
        expand_blast_hits, make_fine_blast_db, read_input_fasta, write_fasta,
    )

    args = types.SimpleNamespace(
        blastp="blastp", makeblastdb="makeblastdb",
        coarse_eval=coarse_evalue, num_workers=os.cpu_count() or 1,
    )
    query_bytes = read_input_fasta(str(query_fasta))
    db = new_read_db(str(db_dir))                       # EXCLUDED from timing
    t0 = time.perf_counter()
    try:
        coarse_xml = _cab._blast_coarse(args, db, query_bytes)
        expanded = expand_blast_hits(db, coarse_xml, args.coarse_eval)
        fasta_bytes = write_fasta(expanded)
        tmp_dir = make_fine_blast_db(args.makeblastdb, fasta_bytes, FILE_BLAST_FINE)
        try:
            proc = subprocess.run(
                ["blastp", "-db", os.path.join(tmp_dir, FILE_BLAST_FINE),
                 "-dbsize", str(db.BlastDBSize),
                 "-num_threads", str(args.num_workers),
                 "-outfmt", f"6 {_TAB_COLS}", "-evalue", str(fine_evalue)],
                input=query_bytes, capture_output=True,
            )
            if proc.returncode != 0:
                raise RuntimeError(
                    "cablastp fine blastp failed:\n"
                    + proc.stderr.decode("utf-8", "replace")
                )
            stdout = proc.stdout.decode("ascii", "replace")
        finally:
            shutil.rmtree(tmp_dir, ignore_errors=True)
    finally:
        db.read_close()
    search_seconds = time.perf_counter() - t0
    rows = [ln.split("\t") for ln in stdout.splitlines() if ln.strip()]
    return search_seconds, parse_tabular_rows(rows)


def run_hs_cablastp_pipeline(
    query_fasta: Path, db_dir: Path, fine_evalue: float, coarse_evalue: float,
) -> Tuple[float, Dict[Pair, float]]:
    """Run hs-cablastp search in-process; return (search_seconds, hits).

    search() reports search_seconds with forest deserialization (load_db)
    already subtracted, so the timer excludes the same DB-load overhead excluded
    for cablastp, and uses the same coarse e-value.
    """
    result = hs_search(
        query_fasta, db_dir,
        coarse_evalue=coarse_evalue, fine_evalue=fine_evalue,
    )
    best: Dict[Pair, float] = {}
    for h in result["fine_hits"]:
        acc = _extract_accession(h.fasta_ref)
        if not acc:
            continue
        key = (_extract_accession(h.qseqid), acc)
        if key not in best or h.bitscore > best[key]:
            best[key] = h.bitscore
    return float(result["search_seconds"]), best


def mean_search_time(
    run_fn, *fn_args, reps: int
) -> Tuple[float, float, Dict[Pair, float]]:
    """Run a pipeline `reps` times; return (mean search_seconds, std, hits).

    Search timing is noisy (BLAST subprocess scheduling), so the single-run
    number is not robust. The hit set is deterministic across reps, so we keep
    the first run's hits and report the mean ± sample std of the per-run search
    times (std = 0.0 for a single rep).
    """
    import statistics

    times: List[float] = []
    hits: Dict[Pair, float] = {}
    for i in range(max(1, reps)):
        secs, h = run_fn(*fn_args)
        times.append(secs)
        if i == 0:
            hits = h
    mean = statistics.mean(times)
    std = statistics.stdev(times) if len(times) > 1 else 0.0
    return mean, std, hits


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
    ap.add_argument("--summary-csv", type=Path, default=None,
                    help="Append a per-pipeline summary row to this CSV "
                         "(created with header if missing).")
    ap.add_argument("--roc-csv", type=Path, default=None,
                    help="Append (fpr, tpr) curve points to this CSV "
                         "(created with header if missing).")
    ap.add_argument("--cab-compress-seconds", type=float, default=None,
                    help="Compression wall time for cablastp (informational, "
                         "written to the summary CSV).")
    ap.add_argument("--hs-compress-seconds", type=float, default=None,
                    help="Compression wall time for hs-cablastp.")
    ap.add_argument("--label", default=None,
                    help="Dataset label written into the CSVs (defaults to "
                         "the fasta filename).")
    ap.add_argument("--timing-reps", type=int, default=3,
                    help="Run each pipeline's search this many times and report "
                         "the median search time (BLAST timing is noisy; default 3). "
                         "Hits are deterministic, so only timing repeats.")
    ap.add_argument("--coarse-evalue", type=float, default=1e-3,
                    help="Coarse-search E-value used by BOTH pipelines (default "
                         "1e-3). cablastp's native default is 5.0 and hs-cablastp's "
                         "is 1e-3; forcing them equal removes a confound where one "
                         "pipeline casts a wider coarse net and so does more/less "
                         "fine-search work (search-time fairness).")
    ap.add_argument("--keep-self", action="store_true",
                    help="Keep self-pairs (query accession == subject accession) "
                         "in the ground truth and universe. By default they are "
                         "excluded: queries are sampled FROM the dataset, so every "
                         "query trivially matches itself at 100%% identity, which "
                         "inflates recall/AUC identically for both pipelines and "
                         "hides their real difference (#2 fairness fix).")
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
    # Original uncompressed DB length. All three pipelines compute e-values
    # against this same effective size (#4): the ground truth searches the full
    # FASTA naturally, cablastp-search pins db.BlastDBSize internally, and
    # hs-cablastp pins params['orig_db_residues']. We compute it here only to
    # report it, so the fairness of the 1e-3 gate is visible in the run log.
    orig_residues = sum(len(seq) for _, seq in records)
    print(f"      original DB residues (e-value -dbsize for all 3): {orig_residues}")

    queries_path = args.out / "queries.fasta"
    write_fasta(queries_path, queries)
    print(f"[1/5] sampled {len(queries)} queries -> {queries_path}")

    # 2) Ground truth: vanilla blastp.
    # Name the raw GT db per source FASTA so switching --fasta never silently
    # reuses a stale ground-truth db built from a different dataset (#1). The old
    # fixed name "raw_blast_db" was the cause of cross-dataset contamination.
    raw_db = args.out / f"raw_blast_db_{args.fasta.stem}"
    ensure_raw_blastdb(args.fasta, raw_db, args.makeblastdb)
    t = time.perf_counter()
    gt_rows = run_blastp_tabular(args.blastp, queries_path, raw_db,
                                 args.ground_truth_evalue)
    gt_scores = parse_tabular_rows(gt_rows)
    gt_positives: Set[Pair] = set(gt_scores.keys())
    print(f"[2/5] vanilla blastp: {len(gt_rows)} rows, "
          f"{len(gt_positives)} unique positive pairs "
          f"({time.perf_counter() - t:.1f}s)")

    # 3) Run cablastp. search_seconds excludes DB open and uses --coarse-evalue;
    # reported as the mean +/- std over --timing-reps runs.
    cab_search_seconds, cab_search_std, cab_scores = mean_search_time(
        run_cablastp_pipeline,
        queries_path, args.cab_db, args.pipeline_evalue, args.coarse_evalue,
        reps=args.timing_reps,
    )
    print(f"[3/5] cablastp: {len(cab_scores)} unique hit pairs "
          f"({cab_search_seconds:.2f}+/-{cab_search_std:.2f}s search, excl. DB load)")

    # 4) Run hs-cablastp. Same timing basis (search excl. DB load) and coarse e-value.
    hs_search_seconds, hs_search_std, hs_scores = mean_search_time(
        run_hs_cablastp_pipeline,
        queries_path, args.hs_db, args.pipeline_evalue, args.coarse_evalue,
        reps=args.timing_reps,
    )
    print(f"[4/5] hs-cablastp: {len(hs_scores)} unique hit pairs "
          f"({hs_search_seconds:.2f}+/-{hs_search_std:.2f}s search, excl. DB load)")

    # Constrain universe to (sampled queries) x (all subjects in fasta).
    # Keys here must match the accession form used by gt_scores / cab_scores /
    # hs_scores, otherwise set differences silently produce empty intersections.
    universe: Set[Pair] = {
        (_extract_accession(q), _extract_accession(s))
        for q in query_ids for s in subject_set
    }

    # Exclude trivial self-pairs unless explicitly kept (#2). A query sampled
    # from the dataset always matches itself; counting that as a recovered
    # positive inflates recall/AUC equally for both pipelines and masks the real
    # gap. Drop q==s from the universe AND every score/positive set so the
    # downstream ROC, recall, missed, and extra counts are all self-free.
    if not args.keep_self:
        def _drop_self(d):
            return {k: v for k, v in d.items() if k[0] != k[1]}
        n_self_gt = sum(1 for q, s in gt_positives if q == s)
        universe = {(q, s) for (q, s) in universe if q != s}
        gt_scores = _drop_self(gt_scores)
        cab_scores = _drop_self(cab_scores)
        hs_scores = _drop_self(hs_scores)
        gt_positives = set(gt_scores.keys())
        print(f"      excluded {n_self_gt} self-pairs from ground truth "
              f"(use --keep-self to retain)")

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

    # 6) CSV output (append rows so multiple datasets accumulate).
    label = args.label or args.fasta.name
    cab_db_kb = _dir_size_kb(args.cab_db)
    hs_db_kb = _dir_size_kb(args.hs_db)

    def _safe_div(a: float, b: float) -> float:
        return a / b if b else 0.0

    # 5b) Headline 3-metric comparison bar chart. This — not the ROC AUC, which
    # saturates near 1.0 on a self-similar corpus — is the figure that actually
    # tells the story: does hs-cablastp give a smaller DB and faster search while
    # holding recall? Each panel is a 2-bar cablastp-vs-hs comparison.
    cab_recall = _safe_div(cab_tp, len(gt_positives))
    hs_recall = _safe_div(hs_tp, len(gt_positives))
    mfig, maxes = plt.subplots(1, 3, figsize=(12, 4.5))
    # Each panel: (title, cab value, hs value, hint, errs). errs is a (cab, hs)
    # std pair for panels with measurement noise (search time), else None.
    panels = [
        ("DB size (KB)", cab_db_kb, hs_db_kb, "lower is better", None),
        ("Search time (s, excl. DB load)", cab_search_seconds, hs_search_seconds,
         f"lower is better; mean±std of {args.timing_reps} runs",
         (cab_search_std, hs_search_std)),
        ("Recall", cab_recall, hs_recall, "higher is better", None),
    ]
    bar_labels = ["cablastp", "hs-cablastp"]
    bar_colors = ["#4C72B0", "#DD8452"]
    for ax, (title, cval, hval, hint, errs) in zip(maxes, panels):
        vals = [cval, hval]
        yerr = list(errs) if errs else None
        bars = ax.bar(bar_labels, vals, color=bar_colors, width=0.6,
                      yerr=yerr, capsize=6, ecolor="#333")
        ax.set_title(f"{title}\n({hint})", fontsize=10)
        ax.grid(True, axis="y", alpha=0.3)
        head = [v + (e or 0.0) for v, e in zip(vals, errs or [0.0, 0.0])]
        top = max(head) if max(head) > 0 else 1.0
        ax.set_ylim(0, top * 1.20)
        for b, v, e in zip(bars, vals, errs or [None, None]):
            if title == "Recall":
                txt = f"{v:.4f}"
            elif e is not None:
                txt = f"{v:.1f}±{e:.1f}"
            else:
                txt = f"{v:.1f}"
            ax.text(b.get_x() + b.get_width() / 2, v + (e or 0.0) + top * 0.02,
                    txt, ha="center", va="bottom", fontsize=9)
    # Annotate the hs/cab ratio for the two "lower is better" panels.
    for ax, (title, cval, hval, _hint, _errs) in list(zip(maxes, panels))[:2]:
        if cval > 0:
            ax.text(0.5, 0.93, f"hs = {hval / cval * 100:.0f}% of cab",
                    transform=ax.transAxes, ha="center", fontsize=8.5,
                    color="#555")
    mfig.suptitle(
        f"hs-cablastp vs cablastp: {args.fasta.name} "
        f"({args.n_queries} queries, self-pairs "
        f"{'kept' if args.keep_self else 'excluded'}, "
        f"coarse e-value {args.coarse_evalue:g}, search excl. DB load)",
        fontsize=11,
    )
    mfig.tight_layout()
    metrics_png = args.out / "metrics.png"
    mfig.savefig(metrics_png, dpi=130)
    print(f"      saved metrics bar chart to {metrics_png}")

    summary_rows = [
        {
            "dataset": label,
            "pipeline": "cablastp",
            "n_queries": args.n_queries,
            "n_subjects": len(subjects),
            "universe_pairs": len(universe),
            "ground_truth_positives": len(gt_positives),
            "compress_seconds": (args.cab_compress_seconds
                                 if args.cab_compress_seconds is not None else ""),
            "db_size_kb": f"{cab_db_kb:.2f}",
            "search_seconds": f"{cab_search_seconds:.3f}",
            "search_seconds_std": f"{cab_search_std:.3f}",
            "total_hit_pairs": len(cab_scores),
            "true_positives": cab_tp,
            "missed_pairs": len(cab_missed),
            "extra_pairs": len(cab_extra),
            "recall": f"{_safe_div(cab_tp, len(gt_positives)):.6f}",
            "precision": f"{_safe_div(cab_tp, len(cab_scores)):.6f}",
            "auc": f"{cab_auc:.6f}",
        },
        {
            "dataset": label,
            "pipeline": "hs-cablastp",
            "n_queries": args.n_queries,
            "n_subjects": len(subjects),
            "universe_pairs": len(universe),
            "ground_truth_positives": len(gt_positives),
            "compress_seconds": (args.hs_compress_seconds
                                 if args.hs_compress_seconds is not None else ""),
            "db_size_kb": f"{hs_db_kb:.2f}",
            "search_seconds": f"{hs_search_seconds:.3f}",
            "search_seconds_std": f"{hs_search_std:.3f}",
            "total_hit_pairs": len(hs_scores),
            "true_positives": hs_tp,
            "missed_pairs": len(hs_missed),
            "extra_pairs": len(hs_extra),
            "recall": f"{_safe_div(hs_tp, len(gt_positives)):.6f}",
            "precision": f"{_safe_div(hs_tp, len(hs_scores)):.6f}",
            "auc": f"{hs_auc:.6f}",
        },
    ]
    if args.summary_csv:
        _append_csv(args.summary_csv, summary_rows)
        print(f"appended 2 rows to {args.summary_csv}")

    if args.roc_csv:
        roc_rows: List[Dict[str, object]] = []
        for fpr, tpr in zip(cab_fpr, cab_tpr):
            roc_rows.append({"dataset": label, "pipeline": "cablastp",
                             "fpr": f"{fpr:.6f}", "tpr": f"{tpr:.6f}"})
        for fpr, tpr in zip(hs_fpr, hs_tpr):
            roc_rows.append({"dataset": label, "pipeline": "hs-cablastp",
                             "fpr": f"{fpr:.6f}", "tpr": f"{tpr:.6f}"})
        _append_csv(args.roc_csv, roc_rows)
        print(f"appended {len(roc_rows)} rows to {args.roc_csv}")

    return 0


def _dir_size_kb(path: Path) -> float:
    total = 0
    for root, _dirs, files in os.walk(path):
        for f in files:
            try:
                total += (Path(root) / f).stat().st_size
            except OSError:
                pass
    return total / 1024.0


def _append_csv(path: Path, rows: List[Dict[str, object]]) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    write_header = not path.exists() or path.stat().st_size == 0
    with open(path, "a", newline="", encoding="ascii") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        if write_header:
            writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    sys.exit(main())
