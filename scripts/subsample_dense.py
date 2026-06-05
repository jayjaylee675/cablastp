#!/usr/bin/env python3
"""Densely subsample a FASTA: concentrate redundancy instead of diluting it.

Uniform random sampling from a large, diverse corpus (e.g. all E. coli TrEMBL)
sprays across thousands of distinct gene families, so the sample has almost no
near-duplicates and the compressor finds little to delta-encode. This script
does the opposite: it groups records by gene name (the UniProt `GN=` header
field) and keeps the members of the most populous families. The same gene
across many strains is near-identical, so each kept family is a tight cluster
of redundant sequences.

A per-family cap spreads the budget over several families (multiple clusters)
rather than letting one huge family fill the whole sample.

Usage:
    python scripts/subsample_dense.py INPUT.fasta OUTPUT.fasta \
        -n 2000 --max-per-family 100 --seed 42
"""

from __future__ import annotations

import argparse
import random
import re
from pathlib import Path

from cablastp.fasta import FastaReader

_GN = re.compile(rb"GN=(\S+)")


def _name_bytes(name) -> bytes:
    return name if isinstance(name, bytes) else name.encode("ascii", "replace")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("infasta")
    ap.add_argument("outfasta")
    ap.add_argument("-n", "--num", type=int, default=2000)
    ap.add_argument("--max-per-family", type=int, default=100,
                    help="cap members taken from any one gene family (0 = unlimited)")
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    rng = random.Random(args.seed)

    # Pass 1: map gene family -> record indices (headers only; cheap).
    families: dict[bytes, list[int]] = {}
    with open(args.infasta, "rb") as fh:
        for i, rec in enumerate(FastaReader(fh)):
            m = _GN.search(_name_bytes(rec.name))
            if m:
                families.setdefault(m.group(1), []).append(i)

    # Rank families largest-first; take a (seeded) quota from each until full.
    cap = args.max_per_family or None
    ranked = sorted(families.items(), key=lambda kv: len(kv[1]), reverse=True)
    wanted: dict[int, int] = {}          # record index -> family rank (for output order)
    chosen_families: list[tuple[bytes, int]] = []
    for rank, (fam, idxs) in enumerate(ranked):
        if len(wanted) >= args.num:
            break
        quota = min(len(idxs), cap) if cap else len(idxs)
        quota = min(quota, args.num - len(wanted))
        picked = idxs if quota >= len(idxs) else sorted(rng.sample(idxs, quota))
        for idx in picked:
            wanted[idx] = rank
        chosen_families.append((fam, len(picked)))

    # Pass 2: pull the sequences for the chosen indices.
    records: dict[int, tuple[str, str]] = {}
    with open(args.infasta, "rb") as fh:
        for i, rec in enumerate(FastaReader(fh)):
            if i in wanted:
                records[i] = (rec.name if isinstance(rec.name, str)
                              else rec.name.decode("ascii", "replace"),
                              rec.residues.decode("ascii").upper())

    # Emit grouped by family (clusters stay contiguous).
    order = sorted(records, key=lambda i: (wanted[i], i))
    out = Path(args.outfasta)
    with open(out, "w") as oh:
        for i in order:
            hdr, s = records[i]
            oh.write(f">{hdr}\n")
            for k in range(0, len(s), 60):
                oh.write(s[k:k + 60] + "\n")

    print(f"wrote {len(order)} sequences from {len(chosen_families)} gene "
          f"families to {out}")
    print("top families used:")
    for fam, cnt in chosen_families[:10]:
        print(f"  {fam.decode('ascii','replace'):20s} {cnt}")


if __name__ == "__main__":
    main()
