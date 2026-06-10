"""Phase 2 (Coarse + Pruning) and Phase 3 (Reconstruction + Fine Search)."""

from __future__ import annotations

import shutil
import subprocess
import sys
import tempfile
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

from hs_cablastp.alignment import apply_edit_script
from hs_cablastp.io import BLAST_COARSE_DB, load_db
from hs_cablastp.types import CompressedDB, EditOp, TreeNode


# Spec hyperparameters.
COARSE_EVALUE = 1e-3
FINE_EVALUE = 1e-10
PRUNE_THRESHOLD = 0.3  # Heuristic score required to keep a child node.


@dataclass
class _RootHit:
    root_id: int
    s_start: int  # On the root's sequence.
    s_end: int


@dataclass
class _CandidateHit:
    """A node retained after pruning, with the implied query/parent ranges."""
    node_id: int
    root_id: int
    inherited_s_range: Tuple[int, int]  # Range on the *root* sequence still covered.


def _run_blastp(
    blastp: str, query: Path, db_prefix: Path, evalue: float,
    dbsize: Optional[int] = None,
) -> List[List[str]]:
    """Run blastp in tabular (outfmt 6) mode and return parsed rows.

    When ``dbsize`` is given it is passed through as blastp's ``-dbsize`` so the
    reported e-values are computed against that effective database length rather
    than the actual (reconstructed-candidate) DB. The fine search uses this to
    pin e-values to the original uncompressed database size, making them match a
    search of the full DB and independent of how many queries are batched.
    """
    cmd = [
        blastp,
        "-query", str(query),
        "-db", str(db_prefix),
        "-evalue", str(evalue),
        "-outfmt", "6 qseqid sseqid pident length mismatch gapopen "
                   "qstart qend sstart send evalue bitscore",
    ]
    if dbsize is not None:
        cmd += ["-dbsize", str(dbsize)]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"blastp failed (exit {proc.returncode}):\n{proc.stderr}")
    rows: List[List[str]] = []
    for line in proc.stdout.splitlines():
        if not line.strip():
            continue
        rows.append(line.split("\t"))
    return rows


def _segment_to_recon_offset(diff_script: List[EditOp], seg_offset: int) -> int:
    """Map a 0-based offset within a child's parent-matched segment to the
    corresponding offset within that child's *reconstructed* sequence.

    A child's reconstruction is ``apply_edit_script(parent_segment, diff_script)``
    (see ``reconstruct_sequence``). INSERTs lengthen and DELETEs shorten the
    output relative to the parent segment, so the two coordinate frames drift
    apart by the net indel count. Pruning passes coverage ranges *down* the tree,
    where a grandchild's ``parent_start``/``parent_end`` are offsets into this
    child's reconstruction — so the inherited range must be expressed in that
    same frame, not the parent-segment frame, or depth>=2 intersections compare
    mismatched coordinates and silently drop legitimate candidates.
    """
    inserts: Dict[int, int] = {}
    deletes: Set[int] = set()
    for op in diff_script:
        if op.op_type == "INSERT":
            inserts[op.position] = inserts.get(op.position, 0) + 1
        elif op.op_type == "DELETE":
            deletes.add(op.position)
    out = 0
    for i in range(seg_offset):
        out += inserts.get(i, 0)
        if i not in deletes:
            out += 1
    out += inserts.get(seg_offset, 0)
    return out


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
    # Translate the intersection into the child's *reconstructed-sequence* frame —
    # the coordinate system its own children's parent_start/parent_end index
    # against. Subtracting pstart gives the offset within the child's matched
    # segment; the diff script then maps that through the child's indels into the
    # reconstruction frame (the two diverge once the child carries INSERTs/DELETEs,
    # which is exactly when the old parent-segment-frame return broke depth>=2).
    seg_a = isect_a - pstart
    seg_b = isect_b - pstart
    new_a = _segment_to_recon_offset(child.diff_script, seg_a)
    new_b = _segment_to_recon_offset(child.diff_script, seg_b)
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
    """Run a full HS-CaBLASTP search; return a dict of results and stats.

    The returned dict carries a timing split so callers can compare *algorithm*
    cost without the one-time DB-load overhead (which is a storage-format
    artifact, not an algorithmic difference): ``load_seconds`` is the time to
    deserialize the forest, ``search_seconds`` is everything after (coarse
    BLAST + pruning + reconstruction + fine BLAST).
    """
    _t0 = time.perf_counter()
    db, params = load_db(db_dir)
    _t_loaded = time.perf_counter()

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
        # Only the subject (root) range drives pruning; the query range and
        # alignment scores in cols 6/7/10/11 are unused here.
        root_hits.append(_RootHit(
            root_id=root_id,
            s_start=int(row[8]) - 1, s_end=int(row[9]),
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
    # Align-once / report-all. Many candidates reconstruct to the *same* residue
    # string: different strain proteins (distinct accessions) that are byte-
    # identical, or a parent and an empty-diff child. The forest already encodes
    # this (a child whose diff_script makes no net change over the hit region is,
    # there, identical to its parent), so re-aligning each copy is wasted work —
    # BLAST returns an identical HSP for an identical subject sequence. We
    # therefore put each *unique* reconstructed sequence into the fine DB once (a
    # "representative") and, after the fine search, copy every HSP on that
    # representative to every node that shares its sequence, each keeping its own
    # ref_original_seq. This shrinks the fine BLAST with zero recall loss: the
    # inherited hit is bit-for-bit what re-aligning the duplicate would have given.
    tmp_dir = Path(tempfile.mkdtemp(prefix="hs_fine_"))
    try:
        fine_fasta = tmp_dir / "fine.fasta"
        seq_to_rep: Dict[str, int] = {}                       # sequence -> representative node id
        rep_members: Dict[int, List[Tuple[int, str]]] = {}    # rep id -> [(node_id, ref), ...]
        with open(fine_fasta, "w", encoding="ascii") as fh:
            for cand in candidates.values():
                seq = reconstruct_sequence(db, cand.node_id)
                if not seq:
                    continue
                ref = db.forest[cand.node_id].ref_original_seq
                rep = seq_to_rep.get(seq)
                if rep is None:
                    # First node carrying this exact sequence: it is the copy that
                    # actually goes into the fine BLAST DB.
                    seq_to_rep[seq] = cand.node_id
                    rep_members[cand.node_id] = [(cand.node_id, ref)]
                    fh.write(f">node_{cand.node_id} {ref}\n")
                    for i in range(0, len(seq), 60):
                        fh.write(seq[i:i + 60] + "\n")
                else:
                    # Identical sequence already queued: this subject will inherit
                    # the representative's HSPs rather than be re-aligned.
                    rep_members[rep].append((cand.node_id, ref))

        if fine_fasta.stat().st_size == 0:
            return {
                "coarse_hits": len(coarse_rows),
                "candidates": len(candidates),
                "unique_candidates": len(seq_to_rep),
                "fine_hits": [],
                "matched_roots": len(matched_root_ids),
                "load_seconds": _t_loaded - _t0,
                "search_seconds": time.perf_counter() - _t_loaded,
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
        # Pin e-values to the original DB size so they match a full-DB search and
        # don't drift with the reconstructed-candidate count / query batching.
        orig_db_residues = params.get("orig_db_residues")
        if orig_db_residues is None:
            # Older DBs predate this key; without it e-values are computed against
            # the tiny candidate DB and are not comparable to a full-DB search.
            # Warn rather than silently change their meaning.
            print(
                "warning: this DB has no 'orig_db_residues'; fine e-values are "
                "computed against the candidate DB, not the original DB size "
                "(rebuild the DB to pin them).",
                file=sys.stderr,
            )
        fine_rows = _run_blastp(
            blastp, query_path, fine_db_prefix, fine_evalue,
            dbsize=orig_db_residues,
        )
        fine_hits: List[SearchHit] = []
        for row in fine_rows:
            sseqid = row[1]
            if not sseqid.startswith("node_"):
                continue
            rep_id = int(sseqid.split("_", 1)[1])
            # Fan this HSP out to every subject that shares the representative's
            # sequence. The alignment stats are identical by construction (same
            # subject residues), so they are copied verbatim; only node_id and the
            # source ref differ per member.
            members = rep_members.get(rep_id, [(rep_id, "")])
            pident = float(row[2]); length = int(row[3])
            qstart = int(row[6]); qend = int(row[7])
            sstart = int(row[8]); send = int(row[9])
            evalue = float(row[10]); bitscore = float(row[11])
            for node_id, ref in members:
                fine_hits.append(SearchHit(
                    node_id=node_id,
                    fasta_ref=ref,
                    qseqid=row[0],
                    pident=pident,
                    length=length,
                    qstart=qstart, qend=qend,
                    sstart=sstart, send=send,
                    evalue=evalue, bitscore=bitscore,
                ))

        return {
            "coarse_hits": len(coarse_rows),
            "matched_roots": len(matched_root_ids),
            "candidates": len(candidates),
            "unique_candidates": len(seq_to_rep),
            "fine_hits": fine_hits,
            "load_seconds": _t_loaded - _t0,
            "search_seconds": time.perf_counter() - _t_loaded,
        }
    finally:
        if not keep_temp:
            shutil.rmtree(tmp_dir, ignore_errors=True)
