"""Disk format for HS-CaBLASTP compressed databases.

Layout under <db_dir>/:
  forest.pkl    -- compact HSFOREST binary blob (see forest_codec.py)
                   NB: extension is legacy; pickle was replaced with a
                   hand-rolled codec that interns ref strings, varint-codes
                   ints, and packs each EditOp into ~3 bytes. Old pickle
                   forests must be rebuilt.
  meta.pkl      -- pickled params (k, min_identity, ...) and node-id counter
  coarse.fasta  -- root sequences, written for makeblastdb consumption
  blastdb-coarse.* -- output of makeblastdb -dbtype prot
"""

from __future__ import annotations

import pickle
import shutil
import subprocess
from pathlib import Path
from typing import Optional

from hs_cablastp.forest_codec import pack_forest, unpack_forest
from hs_cablastp.types import CompressedDB


COARSE_FASTA = "coarse.fasta"
BLAST_COARSE_DB = "blastdb-coarse"
FOREST_FILE = "forest.pkl"
META_FILE = "meta.pkl"


def save_db(
    db: CompressedDB,
    db_dir: Path,
    params: dict,
    makeblastdb_path: str = "makeblastdb",
) -> None:
    db_dir.mkdir(parents=True, exist_ok=True)

    coarse_path = db_dir / COARSE_FASTA
    with open(coarse_path, "w", encoding="ascii") as fh:
        for root in db.roots():
            fh.write(f">node_{root.node_id} {root.ref_original_seq}\n")
            seq = root.sequence or ""
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")

    with open(db_dir / FOREST_FILE, "wb") as fh:
        fh.write(pack_forest(db.forest))
    with open(db_dir / META_FILE, "wb") as fh:
        pickle.dump({"params": params, "next_id": db._next_id}, fh)

    # Build the BLAST coarse DB.
    if any(db.roots()):
        result = subprocess.run(
            [
                makeblastdb_path, "-dbtype", "prot",
                "-in", str(coarse_path),
                "-out", str(db_dir / BLAST_COARSE_DB),
            ],
            capture_output=True, text=True,
        )
        if result.returncode != 0:
            raise RuntimeError(
                f"makeblastdb failed (exit {result.returncode}):\n{result.stderr}"
            )


def _read_root_seqs_from_coarse(path: Path) -> dict[int, str]:
    """Parse coarse.fasta and return {node_id: sequence}.

    Headers in coarse.fasta are `>node_<id> <ref>`, so the node id is the
    first whitespace token after the leading `>`. We rely on that contract
    to round-trip root sequences without duplicating them in the forest blob.
    """
    out: dict[int, str] = {}
    with open(path, "r", encoding="ascii") as fh:
        current: int | None = None
        chunks: list[str] = []
        for line in fh:
            if line.startswith(">"):
                if current is not None:
                    out[current] = "".join(chunks)
                head = line[1:].split(None, 1)[0]
                current = int(head[5:]) if head.startswith("node_") else None
                chunks = []
            else:
                chunks.append(line.strip())
        if current is not None:
            out[current] = "".join(chunks)
    return out


def load_db(db_dir: Path) -> tuple[CompressedDB, dict]:
    with open(db_dir / FOREST_FILE, "rb") as fh:
        forest = unpack_forest(fh.read())
    # Reattach root residues from coarse.fasta (not stored in the forest blob).
    root_seqs = _read_root_seqs_from_coarse(db_dir / COARSE_FASTA)
    for nid, seq in root_seqs.items():
        node = forest.get(nid)
        if node is not None and node.is_root:
            node.sequence = seq
    with open(db_dir / META_FILE, "rb") as fh:
        meta = pickle.load(fh)
    db = CompressedDB()
    db.forest = forest
    db._next_id = meta["next_id"]
    return db, meta["params"]
