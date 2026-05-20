"""Phase 2 (Coarse + Pruning) and Phase 3 (Reconstruction + Fine Search)."""

from __future__ import annotations

import csv
import os
import shutil
import subprocess
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

from hs_cablastp.alignment import apply_edit_script
from hs_cablastp.io import BLAST_COARSE_DB, load_db
from hs_cablastp.types import CompressedDB, EditOp, TreeNode


# Spec hyperparameters.
COARSE_EVALUE = 1e-3
FINE_EVALUE = 1e-10
PRUNE_THRESHOLD = 0.5  # Heuristic score required to keep a child node.


@dataclass
class _RootHit:
    root_id: int
    q_start: int
    q_end: int
    s_start: int  # On the root's sequence.
    s_end: int
    evalue: float
    bitscore: float


@dataclass
class _CandidateHit:
    """A node retained after pruning, with the implied query/parent ranges."""
    node_id: int
    root_id: int
    inherited_s_range: Tuple[int, int]  # Range on the *root* sequence still covered.


def _run_blastp(
    blastp: str, query: Path, db_prefix: Path, evalue: float,
) -> List[List[str]]:
    """Run blastp in tabular (outfmt 6) mode and return parsed rows."""
    cmd = [
        blastp,
        "-query", str(query),
        "-db", str(db_prefix),
        "-evalue", str(evalue),
        "-outfmt", "6 qseqid sseqid pident length mismatch gapopen "
                   "qstart qend sstart send evalue bitscore",
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"blastp failed (exit {proc.returncode}):\n{proc.stderr}")
    rows: List[List[str]] = []
    for line in proc.stdout.splitlines():
        if not line.strip():
            continue
        rows.append(line.split("\t"))
    return rows


def _heuristic_keep(
    child: TreeNode, q_range_on_parent: Tuple[int, int],
    threshold: float,
) -> Tuple[bool, Tuple[int, int]]:
    """Decide whether to keep a child node and return the updated coverage range.

    Logic: count diff_script ops that fall in the *intersection* of the child's
    parent_start..parent_end window with the q_range_on_parent. Many disruptions
    (deletes/substitutions) drop the score.
    """
    pstart, pend = child.parent_start, child.parent_end
    qa, qb = q_range_on_parent
    # Intersect with this child's coverage on the parent.
    isect_a = max(pstart, qa)
    isect_b = min(pend, qb)
    if isect_a >= isect_b:
        return False, (0, 0)

    # Count disruption ops in the intersected window. Positions in diff_script
    # are 0-based offsets within the parent's matched window [pstart, pend).
    window_len = isect_b - isect_a
    disruptions = 0
    for op in child.diff_script:
        abs_pos = pstart + op.position
        if isect_a <= abs_pos < isect_b:
            if op.op_type == "DELETE":
                disruptions += 2  # Deletions hurt more
            else:
                disruptions += 1
    score = 1.0 - (disruptions / max(window_len, 1))
    if score < threshold:
        return False, (0, 0)
    # The child's "new coordinate system" — translate the intersection range into
    # offsets relative to the child's matched segment, which is what its own
    # children will be indexed against.
    new_a = isect_a - pstart
    new_b = isect_b - pstart
    return True, (new_a, new_b)


def _prune_tree(
    db: CompressedDB, root_id: int, q_range_on_root: Tuple[int, int],
    threshold: float,
) -> List[_CandidateHit]:
    """DFS prune the subtree of root_id."""
    out: List[_CandidateHit] = []
    stack: List[Tuple[int, Tuple[int, int]]] = [(root_id, q_range_on_root)]
    while stack:
        node_id, q_range = stack.pop()
        node = db.forest[node_id]
        out.append(_CandidateHit(
            node_id=node_id, root_id=root_id, inherited_s_range=q_range,
        ))
        for child_id in node.children:
            child = db.forest[child_id]
            keep, new_range = _heuristic_keep(child, q_range, threshold)
            if keep:
                stack.append((child_id, new_range))
    return out


def reconstruct_sequence(db: CompressedDB, node_id: int) -> str:
    """Lossless reconstruction from root via diff scripts."""
    chain: List[TreeNode] = []
    cur: Optional[TreeNode] = db.forest[node_id]
    while cur is not None:
        chain.append(cur)
        cur = db.forest[cur.parent_id] if cur.parent_id is not None else None
    chain.reverse()
    root = chain[0]
    seq = root.sequence or ""
    for node in chain[1:]:
        segment = seq[node.parent_start:node.parent_end]
        seq = apply_edit_script(segment, node.diff_script)
    return seq


@dataclass
class SearchHit:
    node_id: int
    fasta_ref: str       # ref_original_seq for traceback to the source protein
    qseqid: str          # query identifier from the input FASTA
    pident: float
    length: int
    qstart: int
    qend: int
    sstart: int
    send: int
    evalue: float
    bitscore: float


def search(
    query_path: Path,
    db_dir: Path,
    *,
    blastp: str = "blastp",
    makeblastdb: str = "makeblastdb",
    coarse_evalue: float = COARSE_EVALUE,
    fine_evalue: float = FINE_EVALUE,
    prune_threshold: float = PRUNE_THRESHOLD,
    keep_temp: bool = False,
    verbose: bool = False,
) -> Dict[str, object]:
    """Run a full HS-CaBLASTP search; return a dict of results and stats."""
    db, _params = load_db(db_dir)

    # --- Phase 2 Step 1: coarse search ---
    coarse_rows = _run_blastp(
        blastp, query_path, db_dir / BLAST_COARSE_DB, coarse_evalue,
    )
    if verbose:
        print(f"coarse hits: {len(coarse_rows)}")

    # Group root hits by root node id; remember query/subject ranges.
    root_hits: List[_RootHit] = []
    matched_root_ids: Set[int] = set()
    for row in coarse_rows:
        sseqid = row[1]                       # e.g. "node_17"
        if not sseqid.startswith("node_"):
            continue
        root_id = int(sseqid.split("_", 1)[1])
        matched_root_ids.add(root_id)
        root_hits.append(_RootHit(
            root_id=root_id,
            q_start=int(row[6]) - 1, q_end=int(row[7]),
            s_start=int(row[8]) - 1, s_end=int(row[9]),
            evalue=float(row[10]), bitscore=float(row[11]),
        ))

    # --- Phase 2 Step 2: tree pruning ---
    candidates: Dict[int, _CandidateHit] = {}
    for h in root_hits:
        for cand in _prune_tree(db, h.root_id, (h.s_start, h.s_end), prune_threshold):
            # Keep the candidate with the widest covered range if duplicated.
            prev = candidates.get(cand.node_id)
            if (
                prev is None
                or (cand.inherited_s_range[1] - cand.inherited_s_range[0])
                > (prev.inherited_s_range[1] - prev.inherited_s_range[0])
            ):
                candidates[cand.node_id] = cand
    if verbose:
        print(f"candidates after pruning: {len(candidates)}")

    # --- Phase 3 Step 1: reconstruct candidate sequences ---
    tmp_dir = Path(tempfile.mkdtemp(prefix="hs_fine_"))
    try:
        fine_fasta = tmp_dir / "fine.fasta"
        node_ref: Dict[int, str] = {}
        with open(fine_fasta, "w", encoding="ascii") as fh:
            for cand in candidates.values():
                seq = reconstruct_sequence(db, cand.node_id)
                if not seq:
                    continue
                node = db.forest[cand.node_id]
                ref = node.ref_original_seq
                node_ref[cand.node_id] = ref
                fh.write(f">node_{cand.node_id} {ref}\n")
                for i in range(0, len(seq), 60):
                    fh.write(seq[i:i + 60] + "\n")

        if fine_fasta.stat().st_size == 0:
            return {
                "coarse_hits": len(coarse_rows),
                "candidates": 0,
                "fine_hits": [],
                "matched_roots": len(matched_root_ids),
            }

        fine_db_prefix = tmp_dir / "blastdb-fine"
        make = subprocess.run(
            [makeblastdb, "-dbtype", "prot",
             "-in", str(fine_fasta), "-out", str(fine_db_prefix)],
            capture_output=True, text=True,
        )
        if make.returncode != 0:
            raise RuntimeError(
                f"fine makeblastdb failed (exit {make.returncode}):\n{make.stderr}"
            )

        # --- Phase 3 Step 2: fine search ---
        fine_rows = _run_blastp(blastp, query_path, fine_db_prefix, fine_evalue)
        fine_hits: List[SearchHit] = []
        for row in fine_rows:
            sseqid = row[1]
            if not sseqid.startswith("node_"):
                continue
            node_id = int(sseqid.split("_", 1)[1])
            fine_hits.append(SearchHit(
                node_id=node_id,
                fasta_ref=node_ref.get(node_id, ""),
                qseqid=row[0],
                pident=float(row[2]),
                length=int(row[3]),
                qstart=int(row[6]), qend=int(row[7]),
                sstart=int(row[8]), send=int(row[9]),
                evalue=float(row[10]), bitscore=float(row[11]),
            ))

        return {
            "coarse_hits": len(coarse_rows),
            "matched_roots": len(matched_root_ids),
            "candidates": len(candidates),
            "fine_hits": fine_hits,
        }
    finally:
        if not keep_temp:
            shutil.rmtree(tmp_dir, ignore_errors=True)
