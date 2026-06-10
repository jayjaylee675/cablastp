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
  - hs-cablastp tree-pruning threshold = **0.3** (the library default).
  - **fine-candidate dedup (align-once / report-all)** is on (see below): identical reconstructed
    candidates are aligned once and the hit is fanned out to every subject sharing that sequence.

## Headline numbers (`summary.csv`)

| dataset | pipeline | DB size (KB) | search (s, excl. DB load) | recall | precision | missed | extra |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| dense_2k   | cablastp    | 1918.1 | 5.9 ± 0.2 | 0.9996 | 1.0000 | 1 | 0 |
| dense_2k   | hs-cablastp | **1546.7** (81%) | **5.6 ± 0.2** (94%) | 0.9978 | 1.0000 | 6 | 0 |
| etrembl_2k | cablastp    | 2467.2 | 4.0 ± 0.2 | 1.0000 | 1.0000 | 0 | 0 |
| etrembl_2k | hs-cablastp | 2392.3 (97%) | **3.1 ± 0.2** (78%) | 0.9884 | 1.0000 | 1 | 0 |

### Fine-candidate dedup (align-once / report-all)

On redundant data many surviving candidates reconstruct to the **same** residue string —
different strain proteins (distinct accessions) that are byte-identical, or a parent and an
empty-diff child. They are *not* a compression failure (those duplicates are stored as 0-byte
diff-scripts, which is why the DB is small) and they cannot be dropped (each is a distinct
subject that a query must be able to find). But re-aligning each identical copy in the fine
BLAST is wasted work: BLAST returns an identical HSP for an identical subject sequence.

`search.py` therefore writes each *unique* reconstructed sequence into the fine BLAST DB once
and, after the search, copies every HSP on that representative to every node sharing its
sequence (each keeps its own `ref_original_seq`). This is **recall-safe by construction** — the
inherited hit is bit-for-bit what re-aligning the duplicate would produce. On dense_2k it cuts
the fine-search set from 2325 to 1606 candidates (−31%) and is what makes hs-cablastp's dense
search *faster* than cablastp rather than slower.

Because dedup removes the candidate redundancy losslessly, the earlier workaround of raising
`--prune-threshold` (which trimmed candidates at a recall cost — and hit a recall cliff on dense
at 0.54: 7→20 missed) is **no longer needed**; the threshold is back at the recall-safe 0.3.

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
- **Search is now faster on both sets** under the algorithm-only measurement (DB-load excluded,
  equal coarse e-value): dense **5.6 vs 5.9 s (94%)** and etrembl **3.1 vs 4.0 s (78%)**. The
  align-once/report-all dedup removed the dense penalty: hs-cablastp previously paid for
  re-aligning many byte-identical strain candidates, and now aligns each unique sequence once.
  (An even earlier writeup over-credited hs's speed because it also included cablastp's DB-load
  overhead; that is excluded here, so this is the conservative, fair number.)
- **Recall is held, just below cablastp**: hs-cablastp misses 6 dense pairs (all on query
  A6MJP2) and 1 etrembl pair. These are hierarchical-fragmentation misses (a subject split across
  forest nodes so no single node carries the full-length HSP), not an e-value artifact, and the
  dedup does not affect them (it is lossless).

## Known limitations / next step

The 6 residual dense misses come from a subject being split across forest nodes so no single
node carries the full-length HSP. The short-flank absorption (commit `37ccdf0`) cut forest
redundancy hard (dense duplicate roots 301→99, roots 928→455) but does not merge these deeper
splits. The next candidate fix is full-length gapped attach (span indels so a near-identical
subject becomes a single child).
