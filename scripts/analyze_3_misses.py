"""Deeper analysis of the 3 remaining missed pairs at PRUNE_THRESHOLD=0.3.

For the 2 Mode-2 (pruning) misses: compute the exact `_heuristic_keep` score
on each subject child node against every coarse hit window, so we can see
how far below 0.3 they actually sit.

For the 1 Mode-1 (HSP-exists-but-not-reported) miss: invoke the live
`hs_cablastp.search.search()` end-to-end and inspect every fine hit it
returns for the subject's nodes — that's the only way to tell whether the
HSP is being dropped before or after `SearchHit` construction.
"""

from __future__ import annotations

import tempfile
from pathlib import Path
from typing import Dict, List

from hs_cablastp.io import BLAST_COARSE_DB, load_db
from hs_cablastp.search import (
    COARSE_EVALUE, PRUNE_THRESHOLD, _run_blastp, search as hs_search,
)


FASTA = Path("data/ecoli_trembl_2k.fasta")
DB_DIR = Path("bench/etrembl2k_hs")

MISSED = [
    # (query_full_id, subject_full_id, gt_bitscore, kind)
    ("tr|A0A0H3PV90|A0A0H3PV90_ECO5C", "tr|A0A5K1QXR9|A0A5K1QXR9_ECOLX", 49.3, "mode1"),
    ("tr|A0A376W501|A0A376W501_ECOLX", "tr|A0A5R8R9U0|A0A5R8R9U0_ECO25", 62.4, "mode2"),
    ("tr|A0A3S4NPX6|A0A3S4NPX6_ECOLX", "tr|A0A5K1QXR9|A0A5K1QXR9_ECOLX", 50.8, "mode2"),
]


def _first_token(s: str) -> str:
    return s.split()[0] if s else ""


def _orig_id_of(node) -> str:
    ref = node.ref_original_seq or ""
    head, sep, tail = ref.rpartition(":")
    return _first_token(head) if head else _first_token(ref)


def _heuristic_score(child, q_range_on_parent):
    """Re-derive the score exactly as `_heuristic_keep` computes it."""
    pstart, pend = child.parent_start, child.parent_end
    qa, qb = q_range_on_parent
    isect_a = max(pstart, qa)
    isect_b = min(pend, qb)
    if isect_a >= isect_b:
        return None, 0, 0
    window_len = isect_b - isect_a
    disruptions = 0
    for op in child.diff_script:
        abs_pos = pstart + op.position
        if isect_a <= abs_pos < isect_b:
            disruptions += 2 if op.op_type == "DELETE" else 1
    score = 1.0 - (disruptions / max(window_len, 1))
    return score, disruptions, window_len


def _load_query_record(target: str) -> str:
    with open(FASTA) as fh:
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
    raise KeyError(target)


def _instrument_multi_query_benchmark(qids: list[str]) -> None:
    """Run the live benchmark's hs_search with the exact 50-query FASTA the
    ROC benchmark sampled, then for each requested subject id, dump every
    fine_hit that mentions it."""
    print("\n" + "#" * 78)
    print("# MULTI-QUERY benchmark replay")
    print("#" * 78)
    queries_path = Path("bench/etrembl2k_p03/queries.fasta")
    if not queries_path.exists():
        print(f"  skipped: {queries_path} doesn't exist")
        return
    result = hs_search(queries_path, DB_DIR, fine_evalue=1e-3)
    fine_hits = result["fine_hits"]
    print(f"  total fine_hits across all 50 queries: {len(fine_hits)}")
    for qid in qids:
        acc_q = qid.split("|")[1] if qid.startswith(("sp|", "tr|")) else qid.split()[0]
        for_q = [h for h in fine_hits if h.qseqid == qid]
        print(f"  query {qid} ({acc_q}):  {len(for_q)} fine_hits returned")
        # Show all distinct fasta_refs and their best bitscores for this query
        per_acc: Dict[str, float] = {}
        for h in for_q:
            ref = (h.fasta_ref or "").split(":", 1)[0]
            acc_s = ref.split("|")[1] if ref.startswith(("sp|", "tr|")) else ref
            if acc_s not in per_acc or h.bitscore > per_acc[acc_s]:
                per_acc[acc_s] = h.bitscore
        print(f"    distinct subject accessions hit: {len(per_acc)}")
        # show specific subjects
        for sid in ["A0A5K1QXR9", "A0A5R8R9U0"]:
            if sid in per_acc:
                print(f"      {sid}: best bitscore {per_acc[sid]:.1f}  <-- PRESENT in hs_scores")
            else:
                print(f"      {sid}: NOT in hs_scores")


def main() -> int:
    db, _ = load_db(DB_DIR)
    print(f"PRUNE_THRESHOLD currently in force: {PRUNE_THRESHOLD}\n")

    # orig_id -> list of subject nodes
    by_orig: Dict[str, List[int]] = {}
    for nid, n in db.forest.items():
        by_orig.setdefault(_orig_id_of(n), []).append(nid)

    coarse_db = DB_DIR / BLAST_COARSE_DB

    for qid, sid, gt_score, kind in MISSED:
        print("=" * 78)
        print(f"[{kind}] {qid}  ->  {sid}  (vanilla bitscore {gt_score})")

        # Run a fresh coarse blast for this single query so we can read
        # all HSPs against the right roots.
        with tempfile.TemporaryDirectory() as tmp:
            qpath = Path(tmp) / "q.fasta"
            qpath.write_text(_load_query_record(qid), encoding="ascii")
            coarse_rows = _run_blastp("blastp", qpath, coarse_db, COARSE_EVALUE)

        subj_nodes = by_orig.get(sid, [])
        subj_child_nodes = [
            db.forest[n] for n in subj_nodes if not db.forest[n].is_root
        ]

        # For each subject child node, find every coarse-hit row that lands
        # on its ancestor root and compute the exact pruning score.
        if kind == "mode2":
            print(f"  Subject's CHILD nodes (the ones at risk of pruning):")
            for child in subj_child_nodes:
                # ancestor root
                cur = child
                while cur.parent_id is not None:
                    cur = db.forest[cur.parent_id]
                root_id = cur.node_id
                # all coarse hits on this root
                row_for_root = [r for r in coarse_rows
                                if r[1] == f"node_{root_id}"]
                print(f"    child node {child.node_id}  parent={child.parent_id}  "
                      f"parent_window={child.parent_start}..{child.parent_end}  "
                      f"ops={len(child.diff_script)} "
                      f"(SUB={sum(1 for o in child.diff_script if o.op_type=='SUBSTITUTE')}, "
                      f"INS={sum(1 for o in child.diff_script if o.op_type=='INSERT')}, "
                      f"DEL={sum(1 for o in child.diff_script if o.op_type=='DELETE')})")
                if not row_for_root:
                    print(f"      no coarse hit on ancestor root {root_id}")
                    continue
                for r in row_for_root:
                    sst, se = int(r[8]) - 1, int(r[9])
                    score, disrupt, window = _heuristic_score(child, (sst, se))
                    if score is None:
                        print(f"      coarse hit on root {root_id} s-range {sst}..{se}: "
                              f"NO INTERSECTION with child window {child.parent_start}..{child.parent_end}")
                        continue
                    pass_ = "PASS" if score >= PRUNE_THRESHOLD else "DROP"
                    print(f"      coarse hit on root {root_id} s-range {sst}..{se}: "
                          f"window={window} disruptions={disrupt} score={score:.3f}  {pass_}")
            # Highest pruning score we'd need to clear to keep any subject child
            best_score = None
            for child in subj_child_nodes:
                cur = child
                while cur.parent_id is not None:
                    cur = db.forest[cur.parent_id]
                root_id = cur.node_id
                for r in coarse_rows:
                    if r[1] != f"node_{root_id}":
                        continue
                    sst, se = int(r[8]) - 1, int(r[9])
                    score, _d, _w = _heuristic_score(child, (sst, se))
                    if score is not None and (best_score is None or score > best_score):
                        best_score = score
            if best_score is not None:
                print(f"  highest pruning score across subject's children: {best_score:.3f}")
                print(f"  -> would need PRUNE_THRESHOLD <= {best_score:.3f} to keep any of them")
            else:
                print(f"  no coarse-hit window intersects any subject child's parent window")
                print(f"  -> pruning is not the gate here; the coarse step doesn't land on the right region")

        # Mode-1: run live search and check what fine_hits has
        if kind == "mode1":
            with tempfile.TemporaryDirectory() as tmp:
                qpath = Path(tmp) / "q.fasta"
                qpath.write_text(_load_query_record(qid), encoding="ascii")
                result = hs_search(qpath, DB_DIR, fine_evalue=1e-3)
            fine_hits = result["fine_hits"]
            print(f"  live search returned {len(fine_hits)} total fine_hits")
            subj_fine_hits = [
                h for h in fine_hits if h.node_id in subj_nodes
            ]
            print(f"  of those, on subject's nodes: {len(subj_fine_hits)}")
            for h in subj_fine_hits:
                print(f"    node_id={h.node_id}  qseqid={h.qseqid!r}  "
                      f"fasta_ref={h.fasta_ref!r}  bitscore={h.bitscore}  evalue={h.evalue}")
            if not subj_fine_hits:
                # Check by orig_id: maybe node id mapping is off
                print(f"  searching all fine_hits for any with fasta_ref pointing at {sid}...")
                from_orig = [h for h in fine_hits if sid.split()[0] in (h.fasta_ref or "")]
                print(f"    {len(from_orig)} hits mention {sid.split()[0]}")
                for h in from_orig[:5]:
                    print(f"      node_id={h.node_id}  fasta_ref={h.fasta_ref!r}  "
                          f"bitscore={h.bitscore}  evalue={h.evalue}")
        print()
    # Replay the multi-query benchmark with all 50 queries so we can see
    # whether the live aggregation drops the surviving fine_hits.
    _instrument_multi_query_benchmark([m[0] for m in MISSED])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
