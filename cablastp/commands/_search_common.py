"""Shared helpers for the search/psisearch/deltasearch CLIs."""

from __future__ import annotations

import io
import os
import sys
import tempfile
from typing import List, TYPE_CHECKING

from cablastp.commands._blast_xml import BlastResults, parse_blast_xml
from cablastp.exec_util import run_command
from cablastp.misc import vprintf

if TYPE_CHECKING:
    from cablastp.db import DB
    from cablastp.seq import OriginalSeq


def split_blast_args(argv: List[str]):
    """Split a CLI argv on '--blast-args'.

    Returns (argv_before, blast_args_after).
    """
    if "--blast-args" not in argv:
        return argv, []
    idx = argv.index("--blast-args")
    return argv[:idx], argv[idx + 1:]


def read_input_fasta(path: str) -> bytes:
    with open(path, "rb") as f:
        return f.read()


def write_fasta(oseqs: List["OriginalSeq"]) -> bytes:
    out = io.BytesIO()
    for oseq in oseqs:
        out.write(b"> " + oseq.name.encode("utf-8") + b"\n")
        out.write(oseq.residues + b"\n")
    return out.getvalue()


def expand_blast_hits(
    db: "DB",
    blast_xml: bytes,
    coarse_eval: float,
) -> List["OriginalSeq"]:
    results: BlastResults = parse_blast_xml(blast_xml)
    used = set()
    oseqs: List["OriginalSeq"] = []
    for hit in results.hits:
        for hsp in hit.hsps:
            try:
                some = db.coarse_db.expand(
                    db.com_db, hit.accession, hsp.hit_from, hsp.hit_to,
                )
            except Exception as err:
                sys.stderr.write(
                    "Could not decompress coarse sequence %d (%d, %d): %s\n"
                    % (hit.accession, hsp.hit_from, hsp.hit_to, err)
                )
                continue
            if hsp.evalue > coarse_eval:
                continue
            for oseq in some:
                if oseq.id in used:
                    continue
                used.add(oseq.id)
                oseqs.append(oseq)
    return oseqs


def make_fine_blast_db(
    makeblastdb: str, fasta_bytes: bytes, fine_db_filename: str,
) -> str:
    tmp_dir = tempfile.mkdtemp(prefix="cablastp-fine-search-db-")
    args = [
        makeblastdb, "-dbtype", "prot",
        "-title", fine_db_filename,
        "-in", "-",
        "-out", os.path.join(tmp_dir, fine_db_filename),
    ]
    vprintf("Created temporary fine BLAST database in %s\n", tmp_dir)
    run_command(args, stdin=fasta_bytes)
    return tmp_dir
