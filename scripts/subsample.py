"""Deterministically sample N sequences from a FASTA file."""
import random
import sys

def parse_fasta(path):
    records = []
    name = None
    chunks = []
    with open(path, "r", encoding="ascii", errors="replace") as fh:
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
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--n", type=int, required=True)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    records = parse_fasta(args.input)
    print(f"Total records: {len(records)}", file=sys.stderr)
    rng = random.Random(args.seed)
    sampled = rng.sample(records, args.n)
    total_res = 0
    with open(args.output, "w", encoding="ascii") as fh:
        for header, seq in sampled:
            fh.write(header + "\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i+60] + "\n")
            total_res += len(seq)
    n_records = sum(1 for h, s in sampled)
    print(f"Sampled {n_records} records, {total_res} residues -> {args.output}", file=sys.stderr)

if __name__ == "__main__":
    main()
