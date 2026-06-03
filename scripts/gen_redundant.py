#!/usr/bin/env python3
"""Generate a small, controlled high-redundancy protein FASTA for clustering tests.

Takes real sequences from an input FASTA as "family founders", then emits each
founder plus several point-mutated variants at a target identity, shuffled. The
result has *known* redundancy: a healthy HS-CaBLASTP compressor should attach the
variants as child nodes under their founder, driving the child/total ratio far
above the 15% gate. If it does not, the clustering problem is algorithmic rather
than data-limited.

Usage:
    python scripts/gen_redundant.py data/redundant_small.fasta \
        --source data/orf_trans_all.fasta \
        --families 20 --variants 9 --identity 0.85 --min-len 80 --seed 42
"""

from __future__ import annotations

import argparse
import random
from pathlib import Path

from cablastp.fasta import FastaReader

AA = "ACDEFGHIKLMNPQRSTVWY"  # 20 standard amino acids


def load_founders(path: Path, n: int, min_len: int) -> list[tuple[str, str]]:
    founders: list[tuple[str, str]] = []
    with open(path, "rb") as fh:
        for rec in FastaReader(fh):
            seq = rec.residues.decode("ascii").upper().replace("*", "")
            # Keep only clean, long-enough sequences so variants can produce a
            # >= MIN_LENGTH match against the founder.
            if len(seq) >= min_len and set(seq) <= set(AA):
                founders.append((rec.name, seq))
            if len(founders) >= n:
                break
    return founders


def mutate(seq: str, identity: float, rng: random.Random) -> str:
    """Return a copy of `seq` with ~(1-identity) of residues substituted."""
    chars = list(seq)
    for i, c in enumerate(chars):
        if rng.random() > identity:
            repl = rng.choice(AA)
            while repl == c:
                repl = rng.choice(AA)
            chars[i] = repl
    return "".join(chars)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("out", nargs="?", default="data/redundant_small.fasta")
    ap.add_argument("--source", default="data/orf_trans_all.fasta")
    ap.add_argument("--families", type=int, default=20)
    ap.add_argument("--variants", type=int, default=9)
    ap.add_argument("--identity", type=float, default=0.85)
    ap.add_argument("--min-len", type=int, default=80)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    rng = random.Random(args.seed)
    founders = load_founders(Path(args.source), args.families, args.min_len)
    if len(founders) < args.families:
        raise SystemExit(
            f"only found {len(founders)} clean founders >= {args.min_len} aa "
            f"in {args.source} (wanted {args.families})"
        )

    records: list[tuple[str, str]] = []
    for fam_idx, (name, seq) in enumerate(founders):
        records.append((f"fam{fam_idx}_founder", seq))
        for v in range(args.variants):
            records.append((f"fam{fam_idx}_var{v}", mutate(seq, args.identity, rng)))

    # Shuffle so a family's members are not adjacent — exercises the seed table.
    rng.shuffle(records)

    out = Path(args.out)
    with open(out, "w") as fh:
        for hdr, s in records:
            fh.write(f">{hdr}\n")
            for i in range(0, len(s), 60):
                fh.write(s[i:i + 60] + "\n")

    ideal = args.families * args.variants  # children if every variant attaches
    print(
        f"wrote {len(records)} sequences "
        f"({args.families} families x {args.variants + 1}) "
        f"@ identity={args.identity} to {out}"
    )
    print(
        f"ground truth: {args.families} founders + {ideal} variants; "
        f"ideal child ratio ~ {ideal / len(records):.2f} "
        f"(every variant clustering under its founder)"
    )


if __name__ == "__main__":
    main()
