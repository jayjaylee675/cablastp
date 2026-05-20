"""Disk format for HS-CaBLASTP compressed databases.

Layout under <db_dir>/:
  forest.pkl    -- pickled CompressedDB.forest dict
  meta.pkl      -- params (k, min_identity, ...) and node-id counter
  coarse.fasta  -- root sequences, written for makeblastdb consumption
  blastdb-coarse.* -- output of makeblastdb -dbtype prot
"""

from __future__ import annotations

import pickle
import shutil
import subprocess
from pathlib import Path
from typing import Optional

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
        pickle.dump(db.forest, fh, protocol=pickle.HIGHEST_PROTOCOL)
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


def load_db(db_dir: Path) -> tuple[CompressedDB, dict]:
    with open(db_dir / FOREST_FILE, "rb") as fh:
        forest = pickle.load(fh)
    with open(db_dir / META_FILE, "rb") as fh:
        meta = pickle.load(fh)
    db = CompressedDB()
    db.forest = forest
    db._next_id = meta["next_id"]
    return db, meta["params"]
