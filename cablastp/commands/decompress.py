"""cablastp-decompress CLI (port of cmd/cablastp-decompress/main.go)."""

from __future__ import annotations

import argparse
import os
import sys
from typing import List, Optional

from cablastp import misc
from cablastp.db import new_read_db
from cablastp.fasta import FastaWriter


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="cablastp-decompress",
        usage="%(prog)s [flags] database-directory out-fasta-file",
    )
    parser.add_argument("-p", type=int, default=os.cpu_count() or 1,
                        dest="num_workers",
                        help="Maximum number of CPUs to use simultaneously.")
    parser.add_argument("--quiet", action="store_true", default=False,
                        help="Only emit errors to stderr.")
    parser.add_argument("--cpuprofile", default="",
                        help="(Ignored - no Go pprof equivalent.)")
    parser.add_argument("--memprofile", default="",
                        help="(Ignored - no Go pprof equivalent.)")
    parser.add_argument("paths", nargs=2,
                        help="database-directory out-fasta-file")
    return parser


def main(argv: Optional[List[str]] = None) -> int:
    if argv is None:
        argv = sys.argv[1:]

    parser = _build_parser()
    args = parser.parse_args(argv)

    if not args.quiet:
        misc.set_verbose(True)

    db_dir, out_fasta_path = args.paths

    out_file = open(out_fasta_path, "w")
    writer = FastaWriter(out_file, asterisk=True)

    try:
        db = new_read_db(db_dir)
    except Exception as err:
        sys.stderr.write("Could not open '%s' database: %s\n" % (db_dir, err))
        return 1
    misc.vprintln("")

    num_seqs = db.com_db.num_sequences()
    for org_seq_id in range(num_seqs):
        try:
            oseq = db.com_db.read_seq(db.coarse_db, org_seq_id)
        except Exception as err:
            sys.stderr.write(
                "Error reading seq id '%d': %s\n" % (org_seq_id, err)
            )
            return 1
        try:
            writer.write(oseq.fasta_seq())
        except Exception as err:
            misc.vprintf("Error writing seq '%s': %s\n", oseq.name, err)

    db.read_close()
    writer.flush()
    out_file.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
