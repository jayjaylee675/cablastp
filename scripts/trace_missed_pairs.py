"""Trace exactly where each missed (query, subject) pair gets dropped.

For each missed pair, replays the hs-cablastp search pipeline step by step:

  1. Coarse search: vanilla blastp(query, coarse FASTA db) at COARSE_EVALUE.
     Did any root in the subject's ancestry chain show up?
  2. Tree pruning: from each coarse-hit root, walk the subtree with the same
     _heuristic_keep logic the live search uses. Did any of the subject's
     forest nodes survive?
  3. Fine DB: reconstruct the surviving candidates' sequences (same code
     path as search.py phase 3) and check whether any of them represents
     the missed subject.
  4. Fine search: vanilla blastp(query, fine DB) at FINE_EVALUE. What is the
     actual reported e-value and bitscore for the (query, subject) pair?

Usage:
    python scripts/trace_missed_pairs.py
"""

from __future__ import annotations

import subprocess
import tempfile
from pathlib import Path
from typing import Dict, List, Set, Tuple

from hs_cablastp.io import BLAST_COARSE_DB, load_db
from hs_cablastp.search import (
    COARSE_EVALUE, PRUNE_THRESHOLD,
    _prune_tree, _run_blastp, reconstruct_sequence,
)

# Match the live benchmark, which passes pipeline_evalue=1e-3 as fine_evalue
# rather than the module default 1e-10. We want the gate as the pipeline
# actually saw it, not the search.py constant.
FINE_EVALUE = 1e-3
from hs_cablastp.types import TreeNode


FASTA = Path("data/ecoli_trembl_2k.fasta")
DB_DIR = Path("bench/etrembl2k_hs")

MISSED = [
    ("tr|A0A0H3PV90|A0A0H3PV90_ECO5C", "tr|A0A5K1QXR9|A0A5K1QXR9_ECOLX", 49.3),
    ("tr|A0A376W501|A0A376W501_ECOLX", "tr|A0A5R8R9U0|A0A5R8R9U0_ECO25", 62.4),
    ("tr|A0A3S4NPX6|A0A3S4NPX6_ECOLX", "tr|A0A5K1QXR9|A0A5K1QXR9_ECOLX", 50.8),
    ("tr|B7MXA2|B7MXA2_ECO81",         "tr|B7MFV9|B7MFV9_ECO45",         71.2),
]


def _first_token(s: str) -> str:
    return s.split()[0] if s else ""


def _load_fasta_record(path: Path, target: str) -> str:
    """Return the FASTA record (header + residues) for `target`'s first token."""
    with open(path) as fh:
        cur_header = None
        chunks: List[str] = []
        for line in fh:
            if line.startswith(">"):
                if cur_header is not None and _first_token(cur_header[1:]) == target:
                    return cur_header + "".join(chunks)
                cur_header = line.rstrip("\n")
                chunks = []
            else:
                chunks.append(line)
        if cur_header is not None and _first_token(cur_header[1:]) == target:
            return cur_header + "".join(chunks)
    raise KeyError(f"no record with first-token {target!r} in {path}")


def _orig_id_of(node: TreeNode) -> str:
    """Extract the orig accession from `ref_original_seq` (rpartition on `:`)."""
    ref = node.ref_original_seq or ""
    head, sep, tail = ref.rpartition(":")
    return _first_token(head) if head else _first_token(ref)


def _ancestor_root_of(db, nid: int) -> int:
    """Walk to the root of the tree containing `nid`."""
    cur = db.forest[nid]
    while cur.parent_id is not None:
        cur = db.forest[cur.parent_id]
    return cur.node_id


def main() -> int:
    db, _ = load_db(DB_DIR)
    forest = db.forest
    # orig_id -> list of node_ids that reference it
    by_orig: Dict[str, List[int]] = {}
    for nid, n in forest.items():
        by_orig.setdefault(_orig_id_of(n), []).append(nid)

    coarse_db = DB_DIR / BLAST_COARSE_DB

    print(f"forest: {len(forest)} nodes")
    print(f"coarse evalue cutoff: {COARSE_EVALUE}   fine evalue cutoff: {FINE_EVALUE}")
    print(f"prune threshold: {PRUNE_THRESHOLD}\n")

    for qid, sid, gt_score in MISSED:
        print("=" * 78)
        print(f"QUERY : {qid}")
        print(f"SUBJ  : {sid}    (vanilla blastp bitscore {gt_score})")

        # Subject's forest nodes + their ancestor roots
        subj_nodes = by_orig.get(sid, [])
        subj_roots = {_ancestor_root_of(db, n) for n in subj_nodes}
        print(f"  subject occupies {len(subj_nodes)} forest nodes, in {len(subj_roots)} root trees")
        for nid in subj_nodes:
            n = forest[nid]
            kind = "ROOT" if n.is_root else f"child(depth={n.depth}, parent={n.parent_id})"
            ref_tail = (n.ref_original_seq or "").rsplit(":", 1)[-1]
            print(f"    node {nid:>5}  {kind:>28}  coords={ref_tail}")

        # --- 1. Coarse search for THIS query ---
        with tempfile.TemporaryDirectory() as tmp:
            qpath = Path(tmp) / "q.fasta"
            qpath.write_text(_load_fasta_record(FASTA, qid), encoding="ascii")
            coarse_rows = _run_blastp("blastp", qpath, coarse_db, COARSE_EVALUE)
            coarse_hit_root_ids: Dict[int, Tuple[int, int, int, int, float, float]] = {}
            for row in coarse_rows:
                ss = row[1]
                if not ss.startswith("node_"):
                    continue
                rid = int(ss.split("_", 1)[1])
                # Keep best (lowest e-value) row per root
                ev = float(row[10])
                qs, qe = int(row[6]) - 1, int(row[7])
                sst, se = int(row[8]) - 1, int(row[9])
                bs = float(row[11])
                prev = coarse_hit_root_ids.get(rid)
                if prev is None or ev < prev[4]:
                    coarse_hit_root_ids[rid] = (qs, qe, sst, se, ev, bs)
            print(f"  coarse hits on {len(coarse_hit_root_ids)} roots total")

            overlap = subj_roots & set(coarse_hit_root_ids)
            print(f"  STEP 1 - subject's ancestor roots that were coarse-hit: "
                  f"{sorted(overlap) if overlap else 'NONE'}")
            if not overlap:
                print(f"  --> DROPPED AT COARSE STEP. fine search never sees this subject.")
                print()
                continue

            # --- 2. Pruning ---
            kept_subj_nodes: Set[int] = set()
            for rid in overlap:
                qs, qe, sst, se, ev, bs = coarse_hit_root_ids[rid]
                cands = _prune_tree(db, rid, (sst, se), PRUNE_THRESHOLD)
                cand_node_ids = {c.node_id for c in cands}
                kept = cand_node_ids & set(subj_nodes)
                kept_subj_nodes |= kept
                print(f"  STEP 2 - root {rid}: coarse hit s-range {sst}..{se}; pruning kept "
                      f"{len(cands)} cands; of subject's nodes: {sorted(kept) if kept else 'NONE'}")
            if not kept_subj_nodes:
                print(f"  --> DROPPED AT PRUNING STEP.")
                print()
                continue

            # --- 3. Build fine DB exactly as search.py would ---
            # Use the same candidate set the live pipeline would have built.
            all_cands: Dict[int, str] = {}  # node_id -> reconstructed sequence
            for rid in coarse_hit_root_ids:
                qs, qe, sst, se, ev, bs = coarse_hit_root_ids[rid]
                cands = _prune_tree(db, rid, (sst, se), PRUNE_THRESHOLD)
                for c in cands:
                    if c.node_id not in all_cands:
                        seq = reconstruct_sequence(db, c.node_id)
                        if seq:
                            all_cands[c.node_id] = seq
            present = {n for n in subj_nodes if n in all_cands}
            print(f"  STEP 3 - fine DB will have {len(all_cands)} candidate sequences total")
            print(f"         - of subject's nodes, reconstructed and present: "
                  f"{sorted(present) if present else 'NONE'}")
            if not present:
                print(f"  --> DROPPED BEFORE FINE DB.")
                print()
                continue

            # --- 4. Run fine blastp ---
            fine_fasta = Path(tmp) / "fine.fasta"
            with open(fine_fasta, "w", encoding="ascii") as fh:
                for nid, seq in all_cands.items():
                    fh.write(f">node_{nid} {forest[nid].ref_original_seq}\n")
                    for i in range(0, len(seq), 60):
                        fh.write(seq[i:i + 60] + "\n")
            fine_db_prefix = Path(tmp) / "fine_db"
            mk = subprocess.run(
                ["makeblastdb", "-dbtype", "prot", "-in", str(fine_fasta),
                 "-out", str(fine_db_prefix)],
                capture_output=True, text=True,
            )
            if mk.returncode != 0:
                print(f"  makeblastdb failed: {mk.stderr[-300:]}")
                print()
                continue

            # Run with a very permissive e-value so we see the row even when
            # the production cutoff would drop it. We'll print both the actual
            # row and whether it would have passed FINE_EVALUE.
            fine_rows = _run_blastp("blastp", qpath, fine_db_prefix, 100.0)
            target_hits = []
            for row in fine_rows:
                ss = row[1]
                if not ss.startswith("node_"):
                    continue
                nid = int(ss.split("_", 1)[1])
                if nid in subj_nodes:
                    target_hits.append((nid, float(row[10]), float(row[11]), int(row[6]), int(row[7]), int(row[8]), int(row[9])))
            print(f"  STEP 4 - fine blastp (permissive) returned {len(target_hits)} HSPs on subject's nodes:")
            for nid, ev, bs, qstart, qend, sstart, send in target_hits:
                pass_fine = "PASS" if ev <= FINE_EVALUE else "FAIL"
                print(f"           node {nid}  q[{qstart}..{qend}] -> s[{sstart}..{send}]  "
                      f"bitscore={bs:.1f}  evalue={ev:.2e}  vs FINE_EVALUE={FINE_EVALUE}: {pass_fine}")
            if not target_hits:
                print(f"  --> fine blastp found NO HSP for query -> subject's nodes "
                      f"(even at evalue<=100). Fragment-vs-full-sequence statistical "
                      f"effect: the subject's residue range in any single retained "
                      f"node is too short to score a detectable local alignment.")
            else:
                passed = [t for t in target_hits if t[1] <= FINE_EVALUE]
                if not passed:
                    print(f"  --> DROPPED AT FINE EVALUE GATE. Best HSP e-value "
                          f"{min(t[1] for t in target_hits):.2e} > FINE_EVALUE {FINE_EVALUE}.")
                else:
                    print(f"  --> UNEXPECTED: {len(passed)} HSP(s) pass FINE_EVALUE; "
                          f"investigate why the live search didn't return it.")
        print()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
