#!/usr/bin/env python3
"""Reservoir-sample N random records from a (possibly huge) FASTA file.

Streams the input so it works on files that don't fit comfortably in memory.

Usage:
    python scripts/subsample_fasta.py INPUT.fasta OUTPUT.fasta -n 1000 --seed 42
"""

from __future__ import annotations

import argparse
import random
from pathlib import Path

from cablastp.fasta import FastaReader


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("infasta")
    ap.add_argument("outfasta")
    ap.add_argument("-n", "--num", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    rng = random.Random(args.seed)
    reservoir: list[tuple[str, str]] = []
    with open(args.infasta, "rb") as fh:
        for i, rec in enumerate(FastaReader(fh)):
            item = (rec.name, rec.residues.decode("ascii").upper())
            if len(reservoir) < args.num:
                reservoir.append(item)
            else:
                j = rng.randint(0, i)
                if j < args.num:
                    reservoir[j] = item

    out = Path(args.outfasta)
    with open(out, "w") as oh:
        for hdr, s in reservoir:
            oh.write(f">{hdr}\n")
            for k in range(0, len(s), 60):
                oh.write(s[k:k + 60] + "\n")
    print(f"wrote {len(reservoir)} sequences to {out}")


if __name__ == "__main__":
    main()
