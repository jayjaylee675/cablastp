# Running hs-cablastp

## Prerequisites

- BLAST+ (`blastp`, `makeblastdb`) – `sudo apt-get install ncbi-blast+`
- Go ≥ 1.18 – used to build the underlying `cablastp-compress` / `cablastp-search`
  binaries.
- Python ≥ 3.8 – for the `hs-cablastp-*` console-script wrappers.

## Installation

```bash
# 1. Build and install the Go binaries
go build -o /usr/local/bin/cablastp-compress ./cmd/cablastp-compress
go build -o /usr/local/bin/cablastp-search   ./cmd/cablastp-search

# 2. Install the Python wrappers
pip install -e .
```

Verify:
```bash
blastp -version
makeblastdb -version
hs-cablastp-compress --help
hs-cablastp-search   --help
```

## Quick start – SwissProt

```bash
# Download
mkdir -p data
curl -L -o data/uniprot_sprot.fasta.gz \
  https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
gunzip data/uniprot_sprot.fasta.gz

# Extract query (human p53, P04637)
python3 - <<'EOF'
import sys
out, capture = open("query.fasta","w"), False
for line in open("data/uniprot_sprot.fasta"):
    if line.startswith(">"):
        capture = "sp|P04637|" in line
    if capture:
        out.write(line)
EOF

# Compress
mkdir -p bench
/usr/bin/time -v hs-cablastp-compress hs_sprot_db data/uniprot_sprot.fasta \
  2>&1 | tee bench/hs_sprot_compress.log

# Search
/usr/bin/time -v hs-cablastp-search hs_sprot_db query.fasta -v \
  2>&1 | tee bench/hs_sprot_search.log
```

## Notes

- Compression is single-threaded by default (`-p 1`). For faster runs use
  `-p <cores>`, but RAM usage scales with parallelism.
- If peak RSS exceeds available memory during a full SwissProt run, compress
  only the first 50 000 sequences with `head -n $(awk …)` and note the fallback.
- The `hs-cablastp-compress` wrapper prints **Forest totals** and output-dir
  size to stderr after the run; these appear in the tee'd log.
- The `hs-cablastp-search` wrapper prints a **Search summary** hit count line
  to stderr after the run.
