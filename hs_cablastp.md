# hs-cablastp

`hs-cablastp` is a Python package that wraps the
[CaBLASTP](https://github.com/BurntSushi/cablastp) compressive protein-search
pipeline with two console-script entry points:

| Command | Wraps | Extra output |
|---------|-------|--------------|
| `hs-cablastp-compress` | `cablastp-compress` | Forest totals (nodes / roots / children) and output-dir size printed to stderr after the run |
| `hs-cablastp-search`   | `cablastp-search`   | Search summary (alignments, HSPs) printed to stderr; accepts `-v` flag (no-op, tool is verbose by default) |

## Algorithm summary

### Compression (`cablastp-compress`)

1. Reads original sequences from a FASTA file one at a time.
2. Looks for k-mer seeds in a hash table built from previously seen sequences.
3. When a seed match is found, extends the alignment with gapped/ungapped
   dynamic programming.
4. If the extension meets the sequence-identity threshold, the query region is
   stored as a *delta* (edit script) pointing into the matching *coarse*
   sequence; otherwise the new sequence becomes a new coarse (root) node.
5. Writes:
   - `coarse.fasta` – the set of root sequences
   - `compressed` – CSV of delta records for child sequences
   - `blastdb-coarse.*` – BLAST database built from coarse sequences

### Search (`cablastp-search`)

1. **Coarse search** – BLASTP query against the small coarse BLAST database
   with a relaxed e-value (default 5.0).
2. **Decompression** – each coarse hit is expanded back to the original
   sequences it encodes; these form the fine candidate set.
3. **Fine search** – BLASTP query against a temporary BLAST database of only
   the candidate sequences, at the user's desired e-value threshold.

The two-stage design typically reduces the volume of sequence that must be
searched at fine resolution by an order of magnitude or more.

## Repository layout

```
hs_cablastp/          Python package (wrappers + entry points)
cmd/cablastp-compress/  Go source for compression binary
cmd/cablastp-search/    Go source for search binary
bench/                  Benchmark logs (gitignored: data/, query.fasta)
data/                   Input FASTA files (gitignored)
```
