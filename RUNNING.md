# Running cablastp and hs-cablastp on real data

This guide walks through compressing a real protein database and running a search query against it with both the classical `cablastp` pipeline (Python port of BurntSushi/cablastp) and the new `hs-cablastp` (hierarchical-subsequence) pipeline.

Both pipelines have the same shape: **compress once, search many times**.

```
FASTA database ──► compress ──► database directory ──► search ──► hits
                                       ▲
                                       │
                                  query FASTA
```

## 1. Prerequisites

### NCBI BLAST+

Both pipelines shell out to the BLAST+ binaries `blastp` and `makeblastdb`. Install them and make sure they're on `PATH`.

```bash
blastp -version
makeblastdb -version
```

If they're not on `PATH`, you'll point at them explicitly via `--blastp` / `--makeblastdb` flags below.

Get BLAST+ from <https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/> (Linux/macOS/Windows binaries).

### Python package

From the repo root:

```bash
pip install -e .
```

That installs the `cablastp` and `hs_cablastp` packages and exposes seven console scripts:

| Script | Purpose |
| --- | --- |
| `cablastp-compress` | Classical compression (Go-port). |
| `cablastp-decompress` | Reconstruct sequences from a classical DB. |
| `cablastp-search` | Classical coarse → fine BLAST search. |
| `cablastp-psisearch` | PSI-BLAST variant. |
| `cablastp-deltasearch` | DELTA-BLAST variant. |
| `hs-cablastp-compress` | Hierarchical-subsequence compression. |
| `hs-cablastp-search` | Hierarchical coarse → prune → fine search. |

Verify the install:

```bash
cablastp-compress --help
hs-cablastp-compress --help
```

## 2. Get a real protein database

Any FASTA file of protein sequences works. Examples to start with:

| Source | Size | URL |
| --- | --- | --- |
| UniProt SwissProt (reviewed) | ~250 MB | `https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz` |
| UniProt TrEMBL (unreviewed, huge) | ~100 GB | same directory, `uniprot_trembl.fasta.gz` |
| NCBI nr | ~200 GB | `https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz` |

For a first run, start with SwissProt. It compresses quickly and gives meaningful results.

```bash
mkdir -p data && cd data
curl -O https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
gunzip uniprot_sprot.fasta.gz
cd ..
```

For a query, grab any single sequence in FASTA format. A 200–500 residue protein is a good size to start:

```bash
# Example: a human kinase
cat > query.fasta <<'EOF'
>sp|P12345|EXAMPLE_HUMAN A query protein
MARKLKVLLQFNAGEDRCYLLEEELDRYRKLLDEEAERLELQEEKKKQALADAA
RKYGLSKKDLAQAVMHLNETFDRDPVRFLENTLSCRCKEPLDHWVAFASIQHFA
EOF
```

## 3. Pipeline A: classical `cablastp`

### 3a. Compress

```bash
cablastp-compress \
    cablastp_sprot_db \
    data/uniprot_sprot.fasta
```

Useful flags:

- `-p N` — number of worker processes (default: number of CPUs).
- `--overwrite` — delete `cablastp_sprot_db` first if it exists.
- `--append` — add to an existing database instead of starting fresh.
- `--min-match-len`, `--match-kmer-size`, `--ext-seq-id-threshold`, … — tune the match heuristics (see `--help`; these mirror the Go flags one-for-one).
- `--makeblastdb /path/to/makeblastdb` — if `makeblastdb` isn't on `PATH`.

The output directory contains binary-encoded `coarse`, `compressed`, `seeds`, and `params` files plus a coarse BLAST DB (`blastdb-coarse.*`).

### 3b. Search

```bash
cablastp-search \
    cablastp_sprot_db \
    query.fasta \
    --blast-args -evalue 1e-10 -outfmt 5 -num_alignments 50
```

Anything after `--blast-args` is forwarded verbatim to the **fine** `blastp` run. The default coarse E-value is 5.0 and is set with `--coarse-eval`. Use `--no-cleanup` to keep the temporary fine BLAST database.

The fine BLAST output goes straight to stdout — pipe it to a file or BLAST XML parser of your choice:

```bash
cablastp-search cablastp_sprot_db query.fasta \
    --blast-args -evalue 1e-10 -outfmt 5 \
    > hits.xml
```

### 3c. Decompress (optional)

If you only have the database and want the original sequences back:

```bash
cablastp-decompress cablastp_sprot_db > restored.fasta
```

## 4. Pipeline B: `hs-cablastp` (hierarchical subsequence)

This pipeline factors each input sequence into a **forest of trees** over locally-similar subsequences. Searching the coarse FASTA only touches the root sequences; tree pruning then decides which descendants are worth materialising for a fine search.

### 4a. Compress

```bash
hs-cablastp-compress \
    hs_sprot_db \
    data/uniprot_sprot.fasta
```

Tunable parameters (defaults from `hs_cablastp/compress.py`):

| Flag | Default | Meaning |
| --- | --- | --- |
| `--k` | 4 | Seed k-mer length. |
| `--min-identity` | 0.7 | Minimum identity to attach a child node to its parent. |
| `--min-length` | 40 | Minimum aligned-region length to attach a child. |
| `--max-depth` | 5 | Maximum tree depth (root counts as 0). |
| `--makeblastdb` | `makeblastdb` | Path to `makeblastdb`. |
| `--overwrite` | off | Delete the output directory first if present. |

After compression you'll see:

```
Forest: 138245 total nodes, 41213 roots, 97032 children
Compression took 312.45s
Database written to /abs/path/hs_sprot_db
```

The directory contains:

```
hs_sprot_db/
  forest.pkl          # pickled forest of TreeNodes
  meta.pkl            # k, min_identity, ... and next_id counter
  coarse.fasta        # one FASTA record per ROOT node ("node_<id>")
  blastdb-coarse.*    # makeblastdb output over coarse.fasta
```

### 4b. Search

```bash
hs-cablastp-search hs_sprot_db query.fasta -v
```

Tunable parameters (defaults from `hs_cablastp/search.py`):

| Flag | Default | Meaning |
| --- | --- | --- |
| `--coarse-evalue` | 1e-3 | E-value for the coarse `blastp` against roots. |
| `--fine-evalue` | 1e-10 | E-value for the fine `blastp` against reconstructed candidates. |
| `--prune-threshold` | 0.5 | Heuristic score required to descend into a child. 0 = keep everything, 1 = keep only perfect matches. |
| `--keep-temp` | off | Keep the temporary fine BLAST DB after the run. |
| `-v` / `--verbose` | off | Print intermediate counts. |
| `--blastp`, `--makeblastdb` | from `PATH` | Override executable locations. |

Example output:

```
coarse hits: 412
candidates after pruning: 287
coarse HSPs:           412
matched root nodes:    103
candidates after prune: 287
fine HSPs:             14
total time:            8.21s

# node_id  source_ref            pident  length  qstart  qend  evalue   bits
node_4118   sp|P12345|EXAMPLE_HUMAN  98.5     193      1   193 4.20e-87   355.0
...
```

`source_ref` is `<original_fasta_id>:<start>-<end>` — the slice of the original input protein this node represents.

## 5. Comparing pipelines

The two pipelines answer different questions:

- **`cablastp`** is faithful to the original Go implementation: same database layout (`coarse`/`compressed`/`seeds` binary files), same tuning flags, same coarse-then-fine BLAST loop. Use it when you want a drop-in replacement for the Go `cablastp`.
- **`hs-cablastp`** is a redesign: the coarse DB only contains *roots* of the subsequence forest, and a deterministic edit-script pruning step shrinks the candidate set before the fine BLAST. Use it when you want to experiment with hierarchical compression or when your data has many near-duplicate locally-similar regions.

The two databases are **not** interchangeable. `cablastp-search` cannot read `hs_sprot_db/` and vice versa.

## 6. Troubleshooting

**`blastp: command not found`** — install BLAST+ or pass `--blastp /full/path/to/blastp`.

**`makeblastdb failed (exit 1)`** — usually a malformed input FASTA. Strip `J`, `O`, `U` residues (`cablastp-compress` already maps these to `X`).

**Compression looks stuck** — for >10 GB FASTAs the first few thousand sequences are slow because the seed table is still small. Use `-p N` (classical) or wait it out (`hs-cablastp` is single-threaded today).

**Out of memory during `hs-cablastp-compress`** — the forest lives entirely in memory and is pickled at the end. For very large inputs, split the FASTA and run multiple shards; merging shards is not yet supported.

**Inconsistent results between runs** — order of seed-table insertion affects which root a child attaches to. Reproducibility requires the same input order.

## 7. Programmatic use

Both packages are importable. Skip the CLI when scripting:

```python
from pathlib import Path
from hs_cablastp.compress import compress_fasta
from hs_cablastp.io import save_db
from hs_cablastp.search import search

comp = compress_fasta(["data/uniprot_sprot.fasta"])
save_db(comp.db, Path("hs_sprot_db"), params={"k": 4})

result = search(Path("query.fasta"), Path("hs_sprot_db"), verbose=True)
for h in result["fine_hits"]:
    print(h.fasta_ref, h.pident, h.evalue)
```
