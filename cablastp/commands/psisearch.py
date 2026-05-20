"""cablastp-psisearch CLI (port of cmd/cablastp-psisearch/main.go)."""

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
        prog="cablastp-psisearch",
        usage=(
            "%(prog)s [flags] database-directory query-fasta-file "
            "[--blast-args BLASTP_ARGUMENTS]"
        ),
    )
    parser.add_argument("--makeblastdb", default="makeblastdb",
                        help="Location of the 'makeblastdb' executable.")
    parser.add_argument("--psiblast", default="psiblast",
                        help="Location of the 'psiblast' executable.")
    parser.add_argument("--coarse-eval", type=float, default=5.0,
                        dest="coarse_eval",
                        help="E-value threshold for the coarse search.")
    parser.add_argument("--no-cleanup", action="store_true", default=False,
                        dest="no_cleanup",
                        help="Keep the temporary fine BLAST database.")
    parser.add_argument("--num_iterations", type=int, default=1,
                        dest="iters",
                        help="Number of PSIBLAST iterations to perform.")
    parser.add_argument("-p", type=int, default=os.cpu_count() or 1,
                        dest="num_workers",
                        help="Maximum number of CPUs to use simultaneously.")
    parser.add_argument("--quiet", action="store_true", default=False,
                        help="Only emit errors to stderr.")
    parser.add_argument("--cpuprofile", default="", help="(Ignored.)")
    parser.add_argument("--memprofile", default="", help="(Ignored.)")
    parser.add_argument("paths", nargs=2,
                        help="database-directory query-fasta-file")
    return parser


def _blast_coarse(args, db, query: bytes) -> bytes:
    cmd = [
        args.psiblast,
        "-db", os.path.join(db.path, FILE_BLAST_COARSE),
        "-outfmt", "5",
        "-num_iterations", str(args.iters),
        "-dbsize", str(db.BlastDBSize),
    ]
    proc = subprocess.Popen(
        cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    out, err = proc.communicate(input=query)
    if proc.returncode != 0:
        raise RuntimeError(
            "psiblast failed: %s" % err.decode("utf-8", errors="replace")
        )
    return out


def _blast_fine(args, db, blast_fine_dir: str, query: bytes,
                blast_args: List[str]) -> int:
    cmd = [
        args.psiblast,
        "-db", os.path.join(blast_fine_dir, FILE_BLAST_FINE),
        "-num_iterations", str(args.iters),
        "-dbsize", str(db.BlastDBSize),
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

    if not args.no_cleanup:
        try:
            shutil.rmtree(tmp_dir)
        except Exception as err:
            sys.stderr.write("Could not delete fine BLAST database: %s\n" % err)

    db.read_close()
    return rc


if __name__ == "__main__":
    sys.exit(main())
