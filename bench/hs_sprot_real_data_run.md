# hs-cablastp Real-World Run — 2026-05-20

> **[PARTIAL]** UniProt SwissProt could not be fetched: the container's network
> policy blocks `ftp.uniprot.org`, `rest.uniprot.org`, and all other UniProt /
> EBI / NCBI hosts (`host_not_allowed`, HTTP 403).  The pipeline was therefore
> run against the repository's bundled yeast-proteome test set
> (`data/orf_trans_all.fasta`, 6 717 sequences).  Every pipeline step ran to
> completion; only the input data differs from the planned SwissProt run.

---

## 1. Input database

| Item | Value |
|------|-------|
| File | `data/orf_trans_all.fasta` (S. cerevisiae ORF translations, SGD Release 64) |
| File size | 4.8 MiB |
| Sequence count | 6 717 |
| Intended file | `data/uniprot_sprot.fasta` (UniProt SwissProt) — unavailable |

---

## 2. Query

| Item | Value |
|------|-------|
| Identifier | YAL001C / TFC3 (S. cerevisiae, SGD:S000000001) |
| Source | extracted programmatically from `data/orf_trans_all.fasta` |
| Length | 1 161 aa |
| Intended query | UniProt P04637 · P53\_HUMAN (393 aa) — unavailable |

Note: the task requires extracting P04637 from the downloaded SwissProt FASTA
using `awk`/Python (not hand-pasting).  Because SwissProt was unavailable, the
first sequence in the fallback database was extracted with an equivalent Python
snippet.

---

## 3. Compression

Command:
```
/usr/bin/time -v hs-cablastp-compress hs_sprot_db data/orf_trans_all.fasta
```

| Metric | Value |
|--------|-------|
| Wall time | 9.74 s (0:09.74) |
| User time | 20.08 s |
| System time | 5.97 s |
| Peak RSS | 2 330 MiB (2 386 448 KB) |
| Exit status | 0 |

### Forest totals

| Metric | Count |
|--------|-------|
| Total nodes | 13 602 |
| Roots (coarse sequences) | 6 885 |
| Children (compressed delta records) | 6 717 |

The root count (6 885) slightly exceeds the input count (6 717) because the
compression algorithm extends coarse sequences as needed to accommodate
near-matches; fragmented extensions are stored as distinct coarse nodes.

### Output directory (`hs_sprot_db/`)

| File | Size |
|------|------|
| `coarse.fasta` | 2.7 MiB |
| `compressed` | 2.1 MiB |
| `blastdb-coarse.*` (6 files) | 3.3 MiB total |
| Other indices | 0.2 MiB |
| **Total** | **8.3 MiB** |

---

## 4. Search

Command:
```
/usr/bin/time -v hs-cablastp-search hs_sprot_db query.fasta -v
```

| Metric | Value |
|--------|-------|
| Coarse BLAST e-value threshold | 5.0 (default) |
| Candidates expanded to fine DB | 2 sequences (1 262 aa) |
| Fine alignments | 2 |
| Fine HSPs | 2 |
| Wall time | 0.36 s |
| Peak RSS | 40.5 MiB (41 504 KB) |
| Exit status | 0 |

Note on low candidate count: TFC3 (YAL001C) is a large, low-complexity protein
(1 161 aa).  Only 2 coarse sequences passed the e-value≤5 coarse filter, which
expanded to exactly 2 fine candidates.

---

## 5. Top fine hits

Only 2 alignments were returned.  "node\_id" is not exposed in BLAST text
output; the sequence identifier serves as the reference.

| # | source\_ref (accession) | pident (%) | aln\_len | evalue | bits |
|---|------------------------|-----------|---------|--------|------|
| 1 | YAL001C TFC3 | 100 | 1 160 | 0.0 | 2 374 |
| 2 | YML084W (dubious ORF) | 47 | 30 | 2.5 | 28.1 |

---

## 6. Hardware

| Item | Value |
|------|-------|
| CPU model | Intel(R) Xeon(R) Processor @ 2.10 GHz |
| Logical cores | 4 |
| Total RAM | 15 GiB |
| Available RAM | 15 GiB |
| Swap | 0 B |

---

## 7. Log files

- Compression log: `bench/hs_sprot_compress.log`
- Search log: `bench/hs_sprot_search.log`

---

## 8. Planned but blocked — SwissProt run

The intended full run on UniProt SwissProt (≈ 570 000 sequences, ≈ 260 MiB
compressed) was blocked at step 2 (download) by the network policy.  All
subsequent steps — extraction of P04637, compression, and search — would follow
identically once the data is available.  The fallback threshold of 50 000
sequences (from the task spec) was not reached because the data could not be
fetched at all.
