# HS-CaBLASTP vs CaBLASTP — final benchmark results

Finalized 2026-06-10. Pipeline state: `ABSORB_SHORT_FLANKS=True`, `OVERLAP_RESIDUES=0`
(boundary overlap removed — short-flank absorption now provides the boundary-recall
protection that overlap padding used to, so removing it no longer regresses recall).
Benchmark methodology: `scripts/roc_benchmark.py` with the fairness fixes.

## Protocol

- **Datasets**: `ecoli_trembl_dense_2k.fasta` (2000 near-identical E. coli strain proteins,
  highly redundant) and `ecoli_trembl_2k.fasta` (2000 uniformly-sampled TrEMBL proteins,
  divergent).
- **Queries**: 50 sampled from each set (seed 42).
- **Ground truth**: vanilla `blastp` of queries vs the full original FASTA, e-value ≤ 1e-3.
- **Fairness**:
  - self-pairs (query == subject) excluded;
  - e-value `-dbsize` pinned to the original DB residue count for all three (GT / cablastp / hs);
  - **same coarse e-value** (1e-3) for both pipelines, so neither casts a wider coarse net
    (cablastp's native default is 5.0, hs's is 1e-3 — equalized so the coarse stage does
    comparable work);
  - **search time excludes the one-time DB load/open** and is measured in-process, so the
    number reflects *algorithm* cost rather than a storage-format artifact (hs deserializes one
    pickle in ~0.03 s; cablastp opens several index files in ~0.5 s — structural, not
    algorithmic). Interpreter startup is excluded for the same reason.
  - reported as **mean ± std over 5 runs** (BLAST scheduling is noisy).

## Headline numbers (`summary.csv`)

| dataset | pipeline | DB size (KB) | search (s, excl. DB load) | recall | precision | missed | extra |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| dense_2k   | cablastp    | 1918.1 | 5.9 ± 0.3 | 0.9996 | 1.0000 | 1 | 0 |
| dense_2k   | hs-cablastp | **1546.7** (81%) | 6.4 ± 0.3 (109%) | 0.9978 | 1.0000 | 6 | 0 |
| etrembl_2k | cablastp    | 2467.2 | 4.0 ± 0.2 | 1.0000 | 1.0000 | 0 | 0 |
| etrembl_2k | hs-cablastp | 2392.3 (97%) | **3.3 ± 0.1** (82%) | 0.9884 | 1.0000 | 1 | 0 |

Figures: `dense_metrics.png`, `dense_roc.png`, `etrembl_metrics.png`, `etrembl_roc.png`.
The 4-metric bar chart (`*_metrics.png`: DB size / search time / recall / precision) is the
headline figure; ROC AUC saturates near 1.0 on this self-similar corpus and does not
discriminate the pipelines.

**Precision is 1.0 for both pipelines by construction**: the pipeline fine-search e-value
(1e-3) equals the ground-truth e-value, so neither pipeline reports any hit outside the ground
truth (extra = 0) — there are no false positives at this threshold. The precision panel
documents that, rather than separating the pipelines. (To make precision discriminating you
would relax the pipeline e-value so sub-threshold hits appear and plot a precision-recall
curve; that pass is intentionally not run here.)

## What the result says

- **DB size is the robust win on redundant data**: on dense_2k the hs-cablastp DB is **81% of
  cablastp's** with essentially equal recall (0.9978 vs 0.9996). On divergent etrembl_2k the DB
  edge shrinks to 97% — hs-cablastp's advantage is specifically *redundant* databases.
- **Search time is a data-dependent trade-off, not a uniform hs win.** Under the algorithm-only
  measurement (DB-load excluded, equal coarse e-value), hs-cablastp is **slightly slower on
  dense** (6.5 vs 6.0 s, +8%) because its forest retains more near-identical fine candidates, and
  **faster on etrembl** (3.2 vs 4.1 s, −22%) because its pruning is more selective on divergent
  data. An earlier writeup reported hs as uniformly faster; that was largely cablastp's DB-load
  overhead and its laxer coarse e-value, both of which are now removed — so this is the honest
  picture.
- **Recall is held**, not perfect: hs-cablastp misses 6 dense pairs (all on one query, A6MJP2)
  and 1 etrembl pair. These are hierarchical-fragmentation misses, not an e-value artifact.

## Known limitations / next step

The 6 residual dense misses come from a subject being split across forest nodes so no single
node carries the full-length HSP. The short-flank absorption (commit `37ccdf0`) cut forest
redundancy hard (dense duplicate roots 301→99, roots 928→455) but does not merge these deeper
splits. The next candidate fix is full-length gapped attach (span indels so a near-identical
subject becomes a single child).
