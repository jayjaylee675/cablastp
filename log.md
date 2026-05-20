# Session log — 2026-05-20

Working log of everything done in this session: commits, documentation,
benchmarks, profiling, and findings. Ordered topically (commits, docs,
results) rather than chronologically.

---

## 1. Commits added this session

Nine commits, branch `master`, branch is **9 ahead of `origin/master`** at the
end of the session.

| Commit | Title |
| --- | --- |
| `5fa18a8` | Ignore Python build/cache artifacts and local Claude settings. |
| `085a350` | Port cablastp from Go to Python. |
| `f782933` | Add HS-CaBLASTP hierarchical-subsequence compression module. |
| `fe36bf0` | Add tests and cablastp-run-with-log helper script. |
| `927c22c` | Add usage guide and HS-CaBLASTP algorithm spec. |
| `92f5fa2` | Expand BLAST+ install instructions in RUNNING.md. |
| `de19a82` | Add cProfile + tracemalloc profiling script for compress runs. |
| `9340e76` | Include qseqid on hs-cablastp SearchHit. |
| `ab6f7cb` | Add ROC benchmark for cablastp vs hs-cablastp. |

Not pushed.

---

## 2. Files added or modified

### Source

- `cablastp/**` (26 files) — Python port of BurntSushi/cablastp.
- `hs_cablastp/**` (8 files) — new hierarchical-subsequence module.
- `hs_cablastp/search.py` — `SearchHit` now carries `qseqid` so batched
  searches can attribute hits to specific queries.
- `pyproject.toml` — packages + 7 console scripts (`cablastp-compress`,
  `cablastp-decompress`, `cablastp-search`, `cablastp-psisearch`,
  `cablastp-deltasearch`, `hs-cablastp-compress`, `hs-cablastp-search`).
- `tests/__init__.py`, `tests/test_cablastp.py`, `tests/test_compress.py`.

### Scripts

- `scripts/cablastp-run-with-log.py` — pre-existing helper.
- `scripts/profile_compress.py` — cProfile + tracemalloc wrapper for either
  compress pipeline. Dumps a `.prof` file and prints top time / memory hot
  spots.
- `scripts/roc_benchmark.py` — samples queries, runs vanilla blastp as
  ground truth, runs both compressed pipelines, computes ROC + AUC, lists
  missed hits, and renders a 2-pane PNG (full + low-FPR zoom).

### Documentation

- `RUNNING.md` — install BLAST+ (with Windows / macOS / Linux notes),
  install the Python package, compress + search end-to-end with both
  pipelines.
- `HS_CABLASTP_ALGORITHM.md` — code-derived spec of the HS-CaBLASTP
  algorithm with every data structure, hyperparameter, and control-flow
  branch noted, including divergences from the original `debug.md` sketch.
- `debug.md`, `hs_cablastp.md` — preserved as algorithm specs alongside
  the new code-derived doc.

### Gitignore

- `.gitignore` — added `__pycache__/`, `*.py[cod]`, `*.egg-info/`,
  `.claude/`, and `bench/` (machine-specific profile / benchmark output).

---

## 3. Performance comparison (compress)

Input: `data/medium.fasta` — 955 yeast ORFs, 684.6 KB.

| Metric | cablastp | hs-cablastp | hs / cab |
| --- | ---: | ---: | --- |
| Wall time | **12.3 s** | **434.5 s** | **35×** slower |
| DB size on disk | 2157 KB (315% of input) | 3221 KB (470% of input) | 1.5× larger |
| Records / nodes | 943 sequences | 1015 nodes (965 roots, 50 children) | — |
| Peak RAM (tracemalloc, small.fasta proxy) | **1896 MB** (pre-allocated) | **8.88 MB** (linear in forest) | hs **213×** smaller |

**Why hs-cablastp barely compresses:** out of 1015 nodes, 965 are roots —
yeast ORFs across gene families are too divergent for the default
`MIN_IDENTITY=0.7` / `MIN_LENGTH=40` thresholds to attach children. The
forest essentially mirrors the input plus a small amount of metadata.

**Why cablastp burns 1.9 GB on tiny input:** `seeds.py:62` allocates
~1.1 GB up front (kmer hash table sized to the alphabet) and
`_memory.py:18` allocates ~760 MB (NW DP arena). Neither grows with
input size — the same ~1.9 GB peak holds for any input.

---

## 4. Performance comparison (search)

Two queries were tested.

### Synthetic query (`query.fasta` from the RUNNING.md example)

| | cablastp | hs-cablastp |
| --- | ---: | ---: |
| Wall time | 4.3 s | 1.4 s |
| Fine hits | 0 | 0 |

Both pipelines correctly find nothing because the example query has no
biological relationship to yeast ORFs. Times are dominated by BLAST
process startup, not the algorithms.

### Real yeast ORF query (`data/yeast_query.fasta` = YBL071C)

| | cablastp | hs-cablastp |
| --- | ---: | ---: |
| Wall time | 3.67 s | 3.63 s |
| Top fine hit | YBL071C, 100% id, 209 bits, E 1.63e-73 | node_221 (YBL071C), 100% id, 209 bits, E 3.08e-77 |

Identical recall, near-identical wall time. Single-query searches are
dominated by BLAST startup overhead in both pipelines.

---

## 5. Profiling findings

Profile runs used `scripts/profile_compress.py` on `data/small.fasta`
(107 ORFs, 78 KB) so the slower pipeline finishes in minutes rather than
an hour.

### `cablastp-compress` — `bench/cablastp_profile.txt`

- Wall: 16.4 s · peak RAM: **1896 MB** · 1.5M function calls.
- Top time hotspots (self time, % of total):
  - `seeds.py:124 hash_kmer` — 5.21 s (32%), 91,413 calls.
  - `_nw.py:16 nw_align` — 2.71 s (16%), 157 calls.
  - `seq.py:26 is_low_complexity` — 1.91 s (12%).
  - `seeds.py:85 Add` — 1.37 s self, 7.84 s cumulative.
- Top allocations: `seeds.py:62` 1129 MiB · `_memory.py:18` 763 MiB.

### `hs-cablastp-compress` — `bench/hs_cablastp_profile.txt`

- Wall: 364.7 s · peak RAM: **8.88 MB** · 33.6M function calls.
- Top time hotspots (self time, % of total):
  - `alignment.py:78 needleman_wunsch` — **191.8 s (52.6%)**, 16,054 calls.
  - `builtins.ord` — 66.3 s (18.2%), 16.8M calls (almost all inside NW).
  - `builtins.max` — 62.8 s (17.2%), 12.7M calls (inside NW DP).
  - `list.append` — 24.1 s (6.6%), 1.97M calls (inside NW DP / traceback).
  - `ungapped_extend` — 5.7 s (1.6%), 16,742 calls.
- **~95% of total time is in or under `needleman_wunsch`.** Memory is a
  non-issue.

### Optimisation priorities implied

For **hs-cablastp**, attack `needleman_wunsch`:

1. Vectorise the DP fill with NumPy or numba (expected 20-50× speedup).
2. Skip NW when ungapped extension already covers `MIN_LENGTH`.
3. Hoist `_RES_INDEX[ord(...)]` lookups out of the inner loop.
4. Replace nested Python lists for the DP table with `array.array` or a
   flat `bytearray`.

For **cablastp**, attack memory:

1. Lazily size the seed table and `Memory` arena to the actual input size
   rather than the worst case.
2. The 32% in `hash_kmer` is per-residue `abs(ord(c) - ord('A'))`; same
   vectorisation case.

Profile dumps preserved at `bench/cab_profile.prof` and
`bench/hs_profile.prof`. Open in `snakeviz` for a flame-graph view.

---

## 6. ROC benchmark — `scripts/roc_benchmark.py`

Setup:

- 50 queries sampled (seed = 42) from `data/medium.fasta` (955 subjects).
- Universe = 50 × 955 = **47,750 (query, subject) pairs**.
- Ground truth = vanilla `blastp` against the raw FASTA, evalue ≤ 1e-3.
- Pipeline evalue = 10.0 (lax, to expose the full score distribution).
- Scoring metric = bitscore.

### Headline numbers

| Pipeline | Hits | TP | **Missed** | Extra | Recall | AUC |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| cablastp | 313 | 98 | **1** | 215 | 0.9899 | **0.9949** |
| hs-cablastp | 380 | 99 | **0** | 281 | 1.0000 | **1.0000** |

- 99 ground-truth positive pairs out of 47,750.
- cablastp missed exactly one pair: `YBR240C → YDL170W`, bitscore 38.1
  (near the noise floor). The coarse representation diverged enough to
  drop the coarse-search HSP. hs-cablastp's forest happened to retain
  enough of YDL170W's residues as a root that the pair survived.
- "Extras" (215 / 281) are pipeline hits below the GT evalue cutoff — not
  false positives, just lower-significance hits. With a stricter pipeline
  cutoff most go away. hs-cablastp's extras outnumber cablastp's because
  its fine-search candidate set (every retained tree node) is larger.

### Caveats

The benchmark stacks the deck in both pipelines' favour:

1. Queries are sampled **from** the dataset. Most positives are
   self-matches and very-close paralogs — easy.
2. Yeast proteome is small and self-similar.
3. Cross-organism queries (e.g. human kinases against SwissProt) would
   exercise compression sensitivity much more aggressively.

### Output artifacts (under `bench/`, git-ignored)

- `roc.png` — full + low-FPR-zoom side-by-side plot.
- `queries.fasta` — the 50 sampled queries (deterministic with seed=42).
- `missed_cablastp.tsv`, `missed_hs_cablastp.tsv` — ground-truth pairs
  each pipeline missed.
- `extra_cablastp.tsv`, `extra_hs_cablastp.tsv` — pipeline hits not in
  the ground truth.
- `roc_run.log` — stdout of the run.

---

## 7. Environment and dependencies pinned during this session

- NCBI BLAST+ 2.17.0 found at `C:\Program Files\NCBI\blast-2.17.0+\bin\`.
  Not on PATH inside the harness — every BLAST-calling command prepends
  this directory to `$env:PATH` for the duration of that command.
- Python 3.13, `pip install -e .` performed for the cablastp package
  (cablastp / hs-cablastp console scripts in
  `C:\Users\user\AppData\Local\Programs\Python\Python313\Scripts`).
- `pip install matplotlib` for ROC plotting (matplotlib 3.10.9 +
  transitive deps).

---

## 8. Files left untracked at end of session

Intentionally not committed (data / local artifacts):

- `data/uniprot_sprot.fasta/` (SwissProt, ~274 MB).
- `data/yeast_query.fasta`, `data/query.fasta`, `query.fasta`.
- `bench/` (profile dumps, ROC outputs, etc.).

The `.gitignore` covers `bench/` automatically; the `data/` files are
left for you to decide whether to track.

---

## 9. Next steps you might want

1. Push the 9 commits — but probably not to `BurntSushi/cablastp`. The
   memory note flags `origin` as Brent's repo; a personal fork or
   private remote is the right destination.
2. Re-run the ROC benchmark on a harder dataset (e.g. `orf_trans_all.fasta`
   with cross-organism queries) to expose a meaningful gap between the
   two pipelines.
3. Tackle the hs-cablastp NW hotspot — vectorising that single function
   should move the compress-time ratio from 35× back to something
   tolerable.
4. Add a `hs-cablastp-decompress` CLI. The forest preserves enough
   information to round-trip every input residue (see
   `HS_CABLASTP_ALGORITHM.md` §7) but no decompressor exists yet.
