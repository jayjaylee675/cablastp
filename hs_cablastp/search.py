"""hs-cablastp-search: thin wrapper around cablastp-search.

Accepts an extra -v / --verbose flag (no-op, cablastp-search is verbose by
default).  Parses the BLAST output that reaches stdout to summarise hit
counts on stderr after the run.
"""

import os
import re
import sys
import subprocess
import tempfile


# Flags consumed by this wrapper; not forwarded to the Go binary.
_WRAPPER_FLAGS = {"-v", "--verbose"}


def _count_blast_hits(blast_output: bytes) -> dict:
    """Return coarse/fine hit counts from BLAST text output.

    For the coarse search the Go tool uses outfmt 5 (XML) internally; we
    never see it.  What reaches stdout is the final fine-search result in
    BLAST's default pairwise text format.  We count the "Sequences producing
    significant alignments" block lines and the HSP count.
    """
    counts = {"alignments": 0, "hsps": 0}
    in_align_block = False
    for raw in blast_output.splitlines():
        line = raw.decode("utf-8", errors="replace")
        if "Sequences producing significant alignments" in line:
            in_align_block = True
            continue
        if in_align_block:
            if line.strip() == "":
                if counts["alignments"] > 0:
                    in_align_block = False
                continue
            if re.match(r"^\S", line):
                counts["alignments"] += 1
        if line.startswith(" Score ="):
            counts["hsps"] += 1
    return counts


def main():
    raw_args = sys.argv[1:]
    args = [a for a in raw_args if a not in _WRAPPER_FLAGS]

    if not args or args[0] in ("-h", "--help"):
        subprocess.run(["cablastp-search"] + args)
        sys.exit(0)

    # Capture stdout so we can both pass it through and collect stats.
    proc = subprocess.run(["cablastp-search"] + args, stdout=subprocess.PIPE)

    blast_out = proc.stdout or b""
    if blast_out:
        sys.stdout.buffer.write(blast_out)
        sys.stdout.buffer.flush()

    if proc.returncode == 0 and blast_out:
        counts = _count_blast_hits(blast_out)
        print(
            f"\nSearch summary: {counts['alignments']} fine alignments, "
            f"{counts['hsps']} HSPs",
            file=sys.stderr,
        )

    sys.exit(proc.returncode)
