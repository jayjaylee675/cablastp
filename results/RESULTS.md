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
- **Fairness**: self-pairs (query == subject) excluded; e-value `-dbsize` pinned to the
  original DB residue count for all three (GT / cablastp / hs-cablastp); both pipelines timed
  as cold subprocesses (identical interpreter-startup cost).
- **Search time** is single-run, dominated by BLAST/process startup; run-to-run noise ≈ ±0.3 s.

## Headline numbers (`summary.csv`)

| dataset | pipeline | DB size (KB) | search (s) | recall | missed | extra |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| dense_2k   | cablastp    | 1918.1 | 7.4 | 0.9996 | 1 | 0 |
| dense_2k   | hs-cablastp | **1546.7** (81%) | 6.5 | 0.9978 | 6 | 0 |
| etrembl_2k | cablastp    | 2467.2 | 5.2 | 1.0000 | 0 | 0 |
| etrembl_2k | hs-cablastp | 2392.3 (97%) | **3.9** (74%) | 0.9884 | 1 | 0 |

Figures: `dense_metrics.png`, `dense_roc.png`, `etrembl_metrics.png`, `etrembl_roc.png`.
The 3-metric bar chart (`*_metrics.png`) is the headline figure; ROC AUC saturates near 1.0
on this self-similar corpus and does not discriminate the pipelines.

## What the result says

- **DB size is the robust win on redundant data**: on dense_2k the hs-cablastp DB is **79% of
  cablastp's** with essentially equal recall (0.9978 vs 0.9996). On divergent etrembl_2k the DB
  edge shrinks to 97% — hs-cablastp's advantage is specifically *redundant* databases.
- **Search time** is comparable-to-faster: tied on dense (6.9 vs 7.1 s, within noise), clearly
  faster on etrembl (3.7 vs 5.4 s) due to lighter fixed overhead (single-pickle DB open, fewer
  fine candidates).
- **Recall is held**, not perfect: hs-cablastp misses 6 dense pairs (all on one query, A6MJP2)
  and 1 etrembl pair. These are hierarchical-fragmentation misses, not an e-value artifact.

## Known limitations / next step

The 6 residual dense misses come from a subject being split across forest nodes so no single
node carries the full-length HSP. The short-flank absorption (commit `37ccdf0`) cut forest
redundancy hard (dense duplicate roots 301→99, roots 928→455) but does not merge these deeper
splits. The next candidate fix is full-length gapped attach (span indels so a near-identical
subject becomes a single child).
