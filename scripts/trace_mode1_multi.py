"""Find why the multi-query fine search drops the A0A0H3PV90 -> A0A5K1QXR9 HSP.

The single-query trace showed `node 2389` (a ROOT covering residues 413-561
of A0A5K1QXR9) carrying a valid HSP at bitscore 33.9, evalue 1.49e-05. The
50-query benchmark replay confirmed the multi-query search returns 52 hits
for this query but A0A5K1QXR9 is not among the subjects. Hypotheses:

  (A) Node 2389 isn't in the multi-query fine DB at all.
  (B) Node 2389 is in the fine DB but blastp's e-value for that HSP is
      higher than 1e-3 because the bigger fine DB inflates the e-value.
  (C) blastp returns the HSP but the SearchHit construction loses it.

Run hs_search with keep_temp=True, then inspect the fine DB directly.
"""

from __future__ import annotations

import subprocess
import tempfile
from pathlib import Path

from hs_cablastp.search import search as hs_search

DB_DIR = Path("bench/etrembl2k_hs")
QUERIES = Path("bench/etrembl2k_p03/queries.fasta")


def main() -> int:
    if not QUERIES.exists():
        print(f"missing {QUERIES}")
        return 1

    # Run the live multi-query search with keep_temp=True so we keep the
    # fine FASTA and the blastdb-fine.* files for inspection.
    result = hs_search(QUERIES, DB_DIR, fine_evalue=1e-3, keep_temp=True)
    fine_hits = result["fine_hits"]
    print(f"hs_search returned {len(fine_hits)} fine_hits")

    # search() doesn't expose its tmp_dir path in the result. Find the
    # most-recently-modified hs_fine_* in TEMP.
    import os
    temp_root = Path(tempfile.gettempdir())
    candidates = sorted(
        (p for p in temp_root.iterdir() if p.name.startswith("hs_fine_")),
        key=lambda p: p.stat().st_mtime, reverse=True,
    )
    if not candidates:
        print("no hs_fine_* tmp_dir found; can't inspect fine DB")
        return 1
    tmp_dir = candidates[0]
    print(f"using tmp_dir: {tmp_dir}")

    fine_fasta = tmp_dir / "fine.fasta"
    fine_db = tmp_dir / "blastdb-fine"
    assert fine_fasta.exists(), fine_fasta

    # Count fine FASTA entries; check if node_2389 is there.
    nodes_in_fine = set()
    with open(fine_fasta) as fh:
        for line in fh:
            if line.startswith(">"):
                head = line[1:].split(None, 1)[0]
                if head.startswith("node_"):
                    nodes_in_fine.add(int(head[5:]))
    print(f"fine DB has {len(nodes_in_fine)} node entries")
    print(f"  node_2389 (A0A5K1QXR9:413-561) in fine DB? {2389 in nodes_in_fine}")
    print(f"  node_2388 (A0A5K1QXR9 child) in fine DB? {2388 in nodes_in_fine}")
    print(f"  node_2387 (A0A5K1QXR9:0-382) in fine DB?  {2387 in nodes_in_fine}")

    # Build a single-query FASTA for A0A0H3PV90 and run blastp directly
    # against the fine DB, with a permissive e-value so we see the HSP
    # regardless of cutoff.
    q_only = tmp_dir / "q_only.fasta"
    with open(QUERIES) as fh, open(q_only, "w") as qfh:
        cur_header = None
        chunks = []
        keep = False
        for line in fh:
            if line.startswith(">"):
                if cur_header is not None and keep:
                    qfh.write(cur_header + "\n" + "".join(chunks))
                head = line[1:].split(None, 1)[0]
                keep = head == "tr|A0A0H3PV90|A0A0H3PV90_ECO5C"
                cur_header = line.rstrip("\n")
                chunks = []
            else:
                if keep:
                    chunks.append(line)
        if cur_header is not None and keep:
            qfh.write(cur_header + "\n" + "".join(chunks))

    proc = subprocess.run(
        [
            "blastp", "-query", str(q_only), "-db", str(fine_db),
            "-evalue", "100", "-outfmt", "6 qseqid sseqid pident length "
            "mismatch gapopen qstart qend sstart send evalue bitscore",
        ],
        capture_output=True, text=True,
    )
    print()
    print(f"manual blastp against multi-query fine DB (evalue<=100):")
    relevant = []
    for ln in proc.stdout.splitlines():
        if not ln.strip():
            continue
        fields = ln.split("\t")
        sseqid = fields[1]
        if sseqid in {"node_2387", "node_2388", "node_2389"}:
            relevant.append(fields)
    print(f"  rows referencing A0A5K1QXR9's nodes (2387/2388/2389): {len(relevant)}")
    for r in relevant:
        print(f"    sseqid={r[1]}  qstart={r[6]} qend={r[7]}  sstart={r[8]} send={r[9]}  "
              f"bitscore={r[11]}  evalue={r[10]}")

    # Same query against fine DB with the live cutoff
    proc = subprocess.run(
        [
            "blastp", "-query", str(q_only), "-db", str(fine_db),
            "-evalue", "1e-3", "-outfmt", "6 qseqid sseqid pident length "
            "mismatch gapopen qstart qend sstart send evalue bitscore",
        ],
        capture_output=True, text=True,
    )
    print()
    print(f"manual blastp against multi-query fine DB (evalue<=1e-3):")
    rows_strict = [ln.split("\t") for ln in proc.stdout.splitlines() if ln.strip()]
    print(f"  total rows: {len(rows_strict)}")
    relevant_strict = [r for r in rows_strict if r[1] in {"node_2387", "node_2388", "node_2389"}]
    print(f"  rows referencing A0A5K1QXR9's nodes: {len(relevant_strict)}")
    for r in relevant_strict:
        print(f"    sseqid={r[1]}  bitscore={r[11]}  evalue={r[10]}")

    print(f"\nfine FASTA: {fine_fasta}   bytes={fine_fasta.stat().st_size}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
