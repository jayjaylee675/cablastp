"""cablastp-deltasearch CLI (port of cmd/cablastp-deltasearch/main.go)."""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from typing import List, Optional

from cablastp import misc
from cablastp.db import new_read_db, FILE_BLAST_COARSE, FILE_BLAST_FINE
from cablastp.commands._search_common import (
    expand_blast_hits,
    make_fine_blast_db,
    read_input_fasta,
    split_blast_args,
    write_fasta,
)


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="cablastp-deltasearch",
        usage=(
            "%(prog)s [flags] --rpspath rpspath database-directory "
            "query-fasta-file [--blast-args BLASTP_ARGUMENTS]"
        ),
    )
    parser.add_argument("--makeblastdb", default="makeblastdb",
                        help="Location of the 'makeblastdb' executable.")
    parser.add_argument("--deltablast", default="deltablast",
                        help="Location of the 'deltablast' executable.")
    parser.add_argument("--rpspath", default="",
                        help="Location of the 'rps' database (required).")
    parser.add_argument("--coarse-eval", type=float, default=5.0,
                        dest="coarse_eval",
                        help="E-value threshold for the coarse search.")
    parser.add_argument("--no-cleanup", action="store_true", default=False,
                        dest="no_cleanup",
                        help="Keep the temporary fine BLAST database.")
    parser.add_argument("-p", type=int, default=os.cpu_count() or 1,
                        dest="num_workers",
                        help="Maximum number of CPUs to use simultaneously.")
    parser.add_argument("--quiet", action="store_true", default=False,
                        help="Only emit errors to stderr.")
    parser.add_argument("--cpuprofile", default="", help="(Ignored.)")
    parser.add_argument("--memprofile", default="", help="(Ignored.)")
    # nargs=2 is intentionally stricter than Go (which silently ignores extra
    # positional args); we reject 3+ rather than quietly drop them.
    parser.add_argument("paths", nargs=2,
                        help="database-directory query-fasta-file")
    return parser


def _blast_coarse(args, db, query: bytes) -> bytes:
    cmd = [
        args.deltablast,
        "-db", os.path.join(db.path, FILE_BLAST_COARSE),
        "-rpsdb", args.rpspath,
        "-outfmt", "5",
        "-dbsize", str(db.BlastDBSize),
    ]
    proc = subprocess.Popen(
        cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    out, err = proc.communicate(input=query)
    if proc.returncode != 0:
        raise RuntimeError(
            "deltablast failed: %s" % err.decode("utf-8", errors="replace")
        )
    return out


def _blast_fine(args, db, blast_fine_dir: str, query: bytes,
                blast_args: List[str]) -> int:
    cmd = [
        args.deltablast,
        "-db", os.path.join(blast_fine_dir, FILE_BLAST_FINE),
        "-rpsdb", args.rpspath,
        "-dbsize", str(db.BlastDBSize),
        "-num_threads", str(args.num_workers),
    ]
    cmd.extend(blast_args)
    proc = subprocess.Popen(
        cmd, stdin=subprocess.PIPE, stdout=sys.stdout, stderr=sys.stderr,
    )
    proc.communicate(input=query)
    return proc.returncode


def main(argv: Optional[List[str]] = None) -> int:
    if argv is None:
        argv = sys.argv[1:]

    argv, blast_args = split_blast_args(argv)
    parser = _build_parser()
    args = parser.parse_args(argv)

    if not args.rpspath:
        sys.stderr.write("The '--rpspath' flag must be set.\n")
        parser.print_help(sys.stderr)
        return 1

    if not args.quiet:
        misc.set_verbose(True)

    db_dir, query_fasta_path = args.paths

    try:
        query_bytes = read_input_fasta(query_fasta_path)
    except Exception as err:
        sys.stderr.write("Could not read input fasta query: %s\n" % err)
        return 1

    try:
        db = new_read_db(db_dir)
    except Exception as err:
        sys.stderr.write("Could not open '%s' database: %s\n" % (db_dir, err))
        return 1

    misc.vprintln("\nBlasting query on coarse database...")
    try:
        coarse_xml = _blast_coarse(args, db, query_bytes)
    except Exception as err:
        sys.stderr.write("Error blasting coarse database: %s\n" % err)
        return 1

    misc.vprintln("Decompressing blast hits...")
    try:
        expanded = expand_blast_hits(db, coarse_xml, args.coarse_eval)
    except Exception as err:
        sys.stderr.write("%s\n" % err)
        return 1

    fasta_bytes = write_fasta(expanded)

    misc.vprintln("Building fine BLAST database...")
    try:
        tmp_dir = make_fine_blast_db(args.makeblastdb, fasta_bytes, FILE_BLAST_FINE)
    except Exception as err:
        sys.stderr.write("Could not create fine database to search on: %s\n" % err)
        return 1

    misc.vprintln("Blasting query on fine database...")
    rc = _blast_fine(args, db, tmp_dir, query_bytes, blast_args)
    if rc != 0:
        # Match Go's fatal "Error blasting fine database": deltablast's own stderr
        # already streamed through, but without this wrapper a non-zero exit was
        # returned silently with no context.
        sys.stderr.write("Error blasting fine database (deltablast exit %d)\n" % rc)

    if not args.no_cleanup:
        try:
            shutil.rmtree(tmp_dir)
        except Exception as err:
            sys.stderr.write("Could not delete fine BLAST database: %s\n" % err)

    db.read_close()
    return rc


if __name__ == "__main__":
    sys.exit(main())
