"""CLI: hs-cablastp-search."""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

from hs_cablastp.search import (
    COARSE_EVALUE, FINE_EVALUE, PRUNE_THRESHOLD, search,
)


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="hs-cablastp-search",
        description="Run an HS-CaBLASTP query against a compressed database.",
    )
    p.add_argument("db_dir", help="HS-CaBLASTP database directory")
    p.add_argument("query", help="Query FASTA file")
    p.add_argument("--blastp", default="blastp", help="Path to the blastp executable")
    p.add_argument("--makeblastdb", default="makeblastdb",
                   help="Path to the makeblastdb executable")
    p.add_argument("--coarse-evalue", type=float, default=COARSE_EVALUE,
                   help=f"Coarse-search E-value (default {COARSE_EVALUE})")
    p.add_argument("--fine-evalue", type=float, default=FINE_EVALUE,
                   help=f"Fine-search E-value (default {FINE_EVALUE})")
    p.add_argument("--prune-threshold", type=float, default=PRUNE_THRESHOLD,
                   help=f"Heuristic pruning threshold 0..1 (default {PRUNE_THRESHOLD})")
    p.add_argument("--keep-temp", action="store_true",
                   help="Keep the temporary fine BLAST database")
    p.add_argument("-v", "--verbose", action="store_true")
    return p


def main(argv=None) -> int:
    args = _build_parser().parse_args(argv)

    t0 = time.perf_counter()
    result = search(
        Path(args.query).resolve(),
        Path(args.db_dir).resolve(),
        blastp=args.blastp,
        makeblastdb=args.makeblastdb,
        coarse_evalue=args.coarse_evalue,
        fine_evalue=args.fine_evalue,
        prune_threshold=args.prune_threshold,
        keep_temp=args.keep_temp,
        verbose=args.verbose,
    )
    elapsed = time.perf_counter() - t0

    print(f"coarse HSPs:           {result['coarse_hits']}")
    print(f"matched root nodes:    {result['matched_roots']}")
    print(f"candidates after prune: {result['candidates']}")
    fine_hits = result["fine_hits"]
    print(f"fine HSPs:             {len(fine_hits)}")
    print(f"total time:            {elapsed:.2f}s")

    if fine_hits:
        print()
        print("# node_id  source_ref            pident  length  qstart  qend  evalue   bits")
        for h in fine_hits[:50]:
            print(
                f"node_{h.node_id:<6} {h.fasta_ref[:22]:22s} "
                f"{h.pident:6.1f} {h.length:6d} {h.qstart:6d} {h.qend:6d} "
                f"{h.evalue:.2e} {h.bitscore:6.1f}"
            )
        if len(fine_hits) > 50:
            print(f"... ({len(fine_hits) - 50} more)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
