"""CLI: hs-cablastp-compress."""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

from hs_cablastp.compress import (
    K_MER_SIZE, MAX_DEPTH, MIN_IDENTITY, MIN_LENGTH, compress_fasta,
)
from hs_cablastp.io import save_db


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="hs-cablastp-compress",
        description="Hierarchical Subsequence CaBLASTP compression.",
    )
    p.add_argument("db_dir", help="Output database directory")
    p.add_argument("fasta", nargs="+", help="One or more input FASTA files")
    p.add_argument("--k", type=int, default=K_MER_SIZE, help=f"Seed k-mer size (default {K_MER_SIZE})")
    p.add_argument("--min-identity", type=float, default=MIN_IDENTITY,
                   help=f"Minimum parent/child identity (default {MIN_IDENTITY})")
    p.add_argument("--min-length", type=int, default=MIN_LENGTH,
                   help=f"Minimum matched-region length (default {MIN_LENGTH})")
    p.add_argument("--max-depth", type=int, default=MAX_DEPTH,
                   help=f"Maximum tree depth (default {MAX_DEPTH})")
    p.add_argument("--makeblastdb", default="makeblastdb",
                   help="Path to the makeblastdb executable")
    p.add_argument("--overwrite", action="store_true",
                   help="Remove the database directory if it exists")
    return p


def main(argv=None) -> int:
    args = _build_parser().parse_args(argv)
    db_dir = Path(args.db_dir).resolve()
    if db_dir.exists():
        if not args.overwrite:
            print(f"error: {db_dir} exists (pass --overwrite)", file=sys.stderr)
            return 2
        import shutil
        shutil.rmtree(db_dir)

    def progress(n, final=False):
        sym = "done" if final else "..."
        print(f"  {sym} compressed {n} sequences", flush=True)

    t0 = time.perf_counter()
    comp = compress_fasta(
        args.fasta,
        k=args.k,
        min_identity=args.min_identity,
        min_length=args.min_length,
        max_depth=args.max_depth,
        progress=progress,
    )
    elapsed = time.perf_counter() - t0

    n_roots = sum(1 for n in comp.db.forest.values() if n.is_root)
    n_nodes = len(comp.db.forest)
    print(f"Forest: {n_nodes} total nodes, {n_roots} roots, {n_nodes - n_roots} children")
    print(f"Compression took {elapsed:.2f}s")

    params = {
        "k": args.k, "min_identity": args.min_identity,
        "min_length": args.min_length, "max_depth": args.max_depth,
        # Original (uncompressed) database size; used by search to set blastp's
        # -dbsize so fine-search e-values match a search of the full DB.
        "orig_db_seqs": comp.input_seqs,
        "orig_db_residues": comp.input_residues,
    }
    save_db(comp.db, db_dir, params, makeblastdb_path=args.makeblastdb)
    print(f"Database written to {db_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
