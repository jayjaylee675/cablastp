"""Deterministically sample N sequences from a FASTA file."""
import argparse
import random
import sys


def parse_fasta(path):
    records = []
    name = None
    chunks = []
    with open(path, "r", encoding="ascii", errors="replace", newline="") as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(chunks)))
                name = line
                chunks = []
            elif line:
                chunks.append(line)
    if name is not None:
        records.append((name, "".join(chunks)))
    return records


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("input")
    ap.add_argument("output")
    ap.add_argument("--n", type=int, required=True)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    records = parse_fasta(args.input)
    print(f"Loaded {len(records)} records from {args.input}", file=sys.stderr)

    rng = random.Random(args.seed)
    sampled = rng.sample(records, args.n)

    with open(args.output, "w", encoding="ascii") as fh:
        for header, seq in sampled:
            fh.write(header + "\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i+60] + "\n")

    print(f"Wrote {len(sampled)} records to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
