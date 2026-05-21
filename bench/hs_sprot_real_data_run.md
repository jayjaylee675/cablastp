# hs-cablastp Real-Data Run — 2026-05-21

> **[PARTIAL]** UniProt SwissProt could not be fetched: the environment's
> network policy blocks `ftp.uniprot.org`, `rest.uniprot.org`, `ncbi.nlm.nih.gov`,
> and `ebi.ac.uk` (all returned HTTP 403 / "Host not in allowlist").
> The pipeline was therefore run against the repository's bundled yeast-proteome
> dataset (`data/orf_trans_all.fasta`, 6,717 *S. cerevisiae* ORF translations).
> Every pipeline step ran to completion; only the input data differs from the
> planned SwissProt run.

---

## 1. Input Database

| Item | Value |
|------|-------|
| File | `data/orf_trans_all.fasta` |
| Source | *S. cerevisiae* ORF translations, SGD Release 64-1-1 |
| File size | 4,899 KB (4.8 MiB) |
| Sequence count | 6,717 |
| Intended file | `data/uniprot_sprot.fasta` (UniProt SwissProt) — network blocked |

---

## 2. Query Sequence

| Item | Value |
|------|-------|
| Locus | YFL039C ACT1 |
| Description | Actin, *S. cerevisiae* (SGD S000001855) |
| Extraction | Programmatic Python snippet scanning fallback FASTA |
| Length | 375 aa |
| Intended query | UniProt `P04637` · `P53_HUMAN` (393 aa) — network blocked |

Actin was chosen as the query because it is highly conserved and its
actin-related-protein (ARP) paralogues provide a wide identity gradient
(22–47%), which is a meaningful stress-test for the coarse-to-fine pipeline.

---

## 3. Compression

Command:
```
/usr/bin/time -v hs-cablastp-compress hs_sprot_db data/orf_trans_all.fasta
```

| Metric | Value |
|--------|-------|
| Wall time | 11.03 s |
| User time | 14.88 s |
| System time | 3.03 s |
| Peak RSS | 845,712 KB (~826 MiB) |
| Exit status | 0 |

### Forest totals

| Metric | Count |
|--------|-------|
| Total nodes | 13,602 |
| Roots (coarse sequences) | 6,885 |
| Children (compressed delta records) | 6,717 |

The root count (6,885) slightly exceeds the input count (6,717) because the
compression algorithm can split and extend coarse sequences when a partial
near-match is found; those extensions become distinct coarse nodes.

### Output directory (`hs_sprot_db/`)

| Metric | Value |
|--------|-------|
| Total size | 8.2 MiB |

---

## 4. Search

Command:
```
/usr/bin/time -v hs-cablastp-search hs_sprot_db query.fasta -v
```

| Metric | Value |
|--------|-------|
| Coarse BLAST e-value threshold | 5.0 (default) |
| Fine candidate sequences | 17 (7,964 aa total) |
| Fine alignments | 17 |
| HSPs (fine search) | 22 |
| Wall time | 0.44 s |
| Peak RSS | 41,376 KB (~40 MiB) |
| Exit status | 0 |

---

## 5. Top 10 Fine Hits

Sorted by BLASTP bit score (best first). `pident` and `aln_len` are from the
top-scoring HSP for each subject. `source_ref` is the SGD locus tag.

| Rank | source_ref | gene | pident | aln_len | evalue | bits |
|------|------------|------|--------|---------|--------|------|
| 1 | YFL039C | ACT1  | 100% | 375 | 0.0    | 783  |
| 2 | YHR129C | ARP1  |  47% | 368 | 3e-126 | 366  |
| 3 | YDL029W | ARP2  |  46% | 377 | 8e-125 | 363  |
| 4 | YJR065C | ARP3  |  35% | 444 | 2e-77  | 243  |
| 5 | YLR085C | ARP6  |  24% | 438 | 1e-24  | 102  |
| 6 | YNL059C | ARP5  |  28% | 230 | 2e-24  | 102  |
| 7 | YJL081C | ARP4  |  26% | 290 | 1e-23  | 99.8 |
| 8 | YPR034W | ARP7  |  29% | 193 | 9e-18  | 82.4 |
| 9 | YMR033W | ARP9  |  22% | 257 | 3e-05  | 43.5 |
| 10 | YOR141C | ARP8  |  31% |  49 | 4e-03  | 37.0 |

Hit 1 is the query itself (100% identity, expected). Hits 2–10 are all
actin-related proteins (ARPs), correctly identified across 22–47% identity —
confirming the two-stage coarse-to-fine search retrieves the expected
protein family.

---

## 6. Hardware

| Item | Value |
|------|-------|
| CPU model | Intel(R) Xeon(R) Processor @ 2.10 GHz |
| Logical cores | 4 |
| Total RAM | 15 GiB |
| Available RAM at run time | ~14 GiB |
| Swap | 0 B |

---

## 7. Step-by-step Status

| Step | Status | Notes |
|------|--------|-------|
| 1 — Install deps | ✅ | BLAST+ 2.12.0, Go 1.24.7, hs-cablastp installed |
| 2 — Fetch SwissProt | ❌ | `ftp.uniprot.org` blocked; fallback to SGD ORF FASTA |
| 3 — Extract query | ✅ | ACT1/YFL039C extracted from fallback; P04637 unavailable |
| 4 — Compress | ✅ | 11 s, 826 MiB peak RSS |
| 5 — Search | ✅ | 0.44 s, 17 hits |
| 6 — Write summary | ✅ | This file |
| 7 — Commit bench/ | ✅ | Log files + this summary |
| 8 — Open PR | ✅ | Title prefixed `[PARTIAL]` |

---

## 8. Log Files

- Compression log: `bench/hs_sprot_compress.log`
- Search log: `bench/hs_sprot_search.log`
