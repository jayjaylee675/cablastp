"""hs-cablastp-compress: thin wrapper around cablastp-compress that adds
forest-topology statistics after the compression run completes."""

import os
import sys
import subprocess


def _forest_totals(db_dir):
    coarse_fasta = os.path.join(db_dir, "coarse.fasta")
    comp_index   = os.path.join(db_dir, "compressed.index")

    roots = 0
    if os.path.exists(coarse_fasta):
        with open(coarse_fasta, "rb") as fh:
            for line in fh:
                if line.startswith(b">"):
                    roots += 1

    children = 0
    if os.path.exists(comp_index):
        # one 8-byte int64 offset per compressed sequence
        children = os.path.getsize(comp_index) // 8

    nodes = roots + children
    return nodes, roots, children


def main():
    args = sys.argv[1:]

    if not args or args[0] in ("-h", "--help"):
        subprocess.run(["cablastp-compress"] + args)
        sys.exit(0)

    proc = subprocess.run(["cablastp-compress"] + args)

    if proc.returncode == 0:
        db_dir = args[0]
        nodes, roots, children = _forest_totals(db_dir)
        print(
            f"\nForest totals: {nodes} nodes | {roots} roots (coarse) "
            f"| {children} children (compressed)",
            file=sys.stderr,
        )
        db_size = sum(
            os.path.getsize(os.path.join(db_dir, f))
            for f in os.listdir(db_dir)
            if os.path.isfile(os.path.join(db_dir, f))
        )
        print(
            f"Output dir size: {db_size / (1024**2):.1f} MiB",
            file=sys.stderr,
        )

    sys.exit(proc.returncode)
