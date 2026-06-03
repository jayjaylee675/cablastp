"""Phase 1: Hierarchical compression and tree building."""

from __future__ import annotations

from collections import OrderedDict
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

from hs_cablastp.alignment import (
    alignment_identity,
    apply_edit_script,
    make_edit_script,
    needleman_wunsch,
    ungapped_extend,
)
from hs_cablastp.types import CompressedDB, EditOp, SeedTable, TreeNode


# Defaults from the spec.
K_MER_SIZE = 4
MIN_IDENTITY = 0.7
MIN_LENGTH = 30           # Step 3: relaxed from 40 to admit shorter homologies.
MAX_DEPTH = 5

# Step 3: depth-dependent identity gating. Matching against a root (depth 0) is
# allowed to be looser because roots are the broadest cluster representatives;
# attaching deeper requires tighter agreement to keep diff scripts small.
MIN_IDENTITY_ROOT = 0.5   # target.depth == 0
MIN_IDENTITY_CHILD = 0.6  # target.depth >= 1

# Implementation knobs (not in the spec).
MAX_SEED_HITS = 32        # cap candidates examined per seed
GAPPED_BAND = 10          # NW band width
EXTEND_DROP = 12          # X-drop for ungapped extension
MISS_STRIDE = 2           # sparse seeding: advance by this many on a failed seed
RECON_CACHE_MAX = 1000    # bounded LRU size for reconstructed-sequence cache


@dataclass
class _Match:
    target_id: int           # The node we attach to (root or descendant)
    target_seq: str          # Reconstructed sequence of the target node
    target_start: int        # Matched range on target_seq
    target_end: int
    q_start: int             # Matched range on the query (input) sequence
    q_end: int
    identity: float
    parent_aln: str
    child_aln: str

    @property
    def length(self) -> int:
        return self.q_end - self.q_start


class _Compressor:
    def __init__(
        self,
        k: int = K_MER_SIZE,
        min_identity: float = MIN_IDENTITY,
        min_length: int = MIN_LENGTH,
        max_depth: int = MAX_DEPTH,
        min_identity_root: float = MIN_IDENTITY_ROOT,
        min_identity_child: float = MIN_IDENTITY_CHILD,
    ):
        self.k = k
        self.min_identity = min_identity
        self.min_length = min_length
        self.max_depth = max_depth
        self.min_identity_root = min_identity_root
        self.min_identity_child = min_identity_child
        self.db = CompressedDB()
        self.seeds = SeedTable(k)
        # Bounded LRU: node_id -> reconstructed string. The original unbounded
        # dict grew with every node ever reconstructed; under a deep/wide forest
        # that is an unbounded heap leak. An OrderedDict capped at RECON_CACHE_MAX
        # keeps memory flat while still caching hot ancestors during a descent.
        self._seq_cache: "OrderedDict[int, str]" = OrderedDict()

    def _min_identity_for_depth(self, depth: int) -> float:
        """Step 3: looser gate for root targets, tighter for deeper nodes."""
        return self.min_identity_root if depth == 0 else self.min_identity_child

    def reconstruct(self, node_id: int) -> str:
        cached = self._seq_cache.get(node_id)
        if cached is not None:
            self._seq_cache.move_to_end(node_id)  # mark most-recently-used
            return cached
        node = self.db.forest[node_id]
        if node.is_root:
            seq = node.sequence or ""
        else:
            parent_seq = self.reconstruct(node.parent_id)  # type: ignore[arg-type]
            parent_segment = parent_seq[node.parent_start:node.parent_end]
            seq = apply_edit_script(parent_segment, node.diff_script)
        self._seq_cache[node_id] = seq
        self._seq_cache.move_to_end(node_id)
        if len(self._seq_cache) > RECON_CACHE_MAX:
            self._seq_cache.popitem(last=False)  # evict least-recently-used
        return seq

    def _invalidate_cache(self, node_id: int) -> None:
        # Children inherit from parent; updating a parent's sequence (only possible
        # for roots which we never mutate) would force cache eviction.
        self._seq_cache.pop(node_id, None)

    # --- Phase 1 core ---

    def compress_sequence(self, fasta_id: str, sequence: str) -> None:
        seq = sequence
        N = len(seq)
        if N == 0:
            return
        # Too short to seed-match; preserve the whole sequence as a root.
        if N < self.k:
            self._make_root(fasta_id, seq, 0, N)
            return
        pos = 0
        pending_root_start = 0
        while pos + self.k <= N:
            kmer = seq[pos:pos + self.k]
            hits = self.seeds.lookup(kmer)
            match: Optional[_Match] = None
            if hits:
                match = self._best_match(seq, pos, hits)

            if match is not None:
                # Flush any preceding unmatched region as a new root.
                # MUST preserve every leftover residue (debug.md "CRITICAL").
                if match.q_start > pending_root_start:
                    self._make_root(
                        fasta_id, seq, pending_root_start, match.q_start
                    )
                # Attach child to the best target (possibly deeper than root).
                target_node = self.db.forest[match.target_id]
                if target_node.depth < self.max_depth:
                    diff_script = make_edit_script(match.parent_aln, match.child_aln)
                    child = self.db.new_node(
                        is_root=False,
                        sequence=None,
                        parent_id=match.target_id,
                        diff_script=diff_script,
                        depth=target_node.depth + 1,
                        ref_original_seq=f"{fasta_id}:{match.q_start}-{match.q_end}",
                        parent_start=match.target_start,
                        parent_end=match.target_end,
                    )
                    target_node.children.append(child.node_id)
                    # Step 2: index the child's own k-mers so later sequences
                    # can seed-match it directly (deeper clustering) instead of
                    # only ever attaching to roots.
                    child_seq = seq[match.q_start:match.q_end]
                    self.seeds.add_node_sequence(child.node_id, child_seq)
                else:
                    # Depth cap reached: emit this region as a new root instead.
                    self._make_root(
                        fasta_id, seq, match.q_start, match.q_end
                    )
                pos = match.q_end
                pending_root_start = pos
            else:
                # Sparse seeding: when the window at `pos` yields no usable
                # match, neighbouring windows (which share k-1 residues) almost
                # always fail too, so re-running extension/NW at every offset is
                # wasteful. Stride forward to skip the redundant work. Data is
                # still preserved: the skipped residues stay inside the pending
                # root region and are flushed by _make_root below.
                pos += MISS_STRIDE if hits else 1

        # Flush trailing region. MUST preserve every leftover residue.
        if N > pending_root_start:
            self._make_root(fasta_id, seq, pending_root_start, N)

    def _make_root(self, fasta_id: str, seq: str, start: int, end: int) -> TreeNode:
        sub = seq[start:end]
        node = self.db.new_node(
            is_root=True,
            sequence=sub,
            parent_id=None,
            depth=0,
            ref_original_seq=f"{fasta_id}:{start}-{end}",
        )
        self.seeds.add_root_sequence(node.node_id, sub)
        return node

    def _best_match(
        self, seq: str, pos: int, hits: List[Tuple[int, int]]
    ) -> Optional[_Match]:
        best: Optional[_Match] = None
        for root_id, root_off in hits[:MAX_SEED_HITS]:
            # Targets may now be child nodes (Step 2), whose sequence is not
            # stored literally — reconstruct it (cached LRU) rather than reading
            # node.sequence, which is None for children.
            root_seq = self.reconstruct(root_id)
            if not root_seq:
                continue
            min_identity = self._min_identity_for_depth(
                self.db.forest[root_id].depth
            )
            sa, ea, sb, eb = ungapped_extend(
                seq, pos, root_seq, root_off, self.k, max_drop=EXTEND_DROP,
            )
            q_sub = seq[sa:ea]
            t_sub = root_seq[sb:eb]
            if len(q_sub) >= self.min_length:
                # ungapped_extend advances both pointers in lockstep, so q_sub
                # and t_sub are equal-length and gap-free: the alignment is the
                # 1:1 pairing itself. No Needleman-Wunsch needed (this was the
                # dominant cost in the profile). Identity is a direct compare.
                parent_aln, child_aln = t_sub, q_sub
            else:
                # Too short ungapped: fall back to constrained gapped alignment.
                q_sub, t_sub, sa, ea, sb, eb = self._gapped_extend(
                    seq, root_seq, sa, ea, sb, eb
                )
                if len(q_sub) < self.min_length:
                    continue
                # NB: an older revision gated this call with a cheap
                # best_diagonal_identity prefilter to skip NW on candidates
                # unlikely to reach the threshold. Profiling after parasail
                # was wired in (~0.05 ms/NW vs ~1 ms/prefilter scan) showed
                # the prefilter cost more than the NW it avoided, *and* its
                # band=4 gap-free assumption silently rejected real homologs
                # with wider indels. With parasail on, run NW unconditionally.
                parent_aln, child_aln = needleman_wunsch(
                    t_sub, q_sub, band=GAPPED_BAND
                )
            identity = alignment_identity(parent_aln, child_aln)
            if identity < min_identity:
                continue

            cand = _Match(
                target_id=root_id,
                target_seq=root_seq,
                target_start=sb,
                target_end=eb,
                q_start=sa,
                q_end=ea,
                identity=identity,
                parent_aln=parent_aln,
                child_aln=child_aln,
            )
            # Try to refine: walk into existing children to see if any aligns better.
            cand = self._descend_to_best(cand)
            if (
                best is None
                or cand.identity * cand.length > best.identity * best.length
            ):
                best = cand
        return best

    def _gapped_extend(
        self, seq: str, root_seq: str, sa: int, ea: int, sb: int, eb: int
    ) -> Tuple[str, str, int, int, int, int]:
        # Expand the window symmetrically until min_length residues are covered.
        target_len = self.min_length
        cur_len = max(ea - sa, eb - sb)
        if cur_len < target_len:
            extra = target_len - cur_len
            sa = max(0, sa - extra // 2)
            sb = max(0, sb - extra // 2)
            ea = min(len(seq), sa + target_len)
            eb = min(len(root_seq), sb + target_len)
        return seq[sa:ea], root_seq[sb:eb], sa, ea, sb, eb

    def _descend_to_best(self, match: _Match) -> _Match:
        """Walk into the target node's subtree to find a descendant whose
        reconstructed segment aligns at least as well as the current target."""
        best = match
        current_id = match.target_id
        while True:
            current = self.db.forest[current_id]
            if not current.children or current.depth >= self.max_depth - 1:
                break
            improved = False
            for child_id in current.children:
                child = self.db.forest[child_id]
                child_seq = self.reconstruct(child_id)
                if not child_seq:
                    continue
                # Align the query subsequence against child_seq.
                q_sub = match_q_sub(best)
                if len(child_seq) == len(q_sub):
                    # Equal length => a gap-free 1:1 pairing is itself a valid
                    # alignment, so the banded NW (the profile's dominant cost,
                    # reached almost entirely through this descent loop) cannot
                    # change the residue correspondence. Skip it. This mirrors
                    # the same gap-free fast path in _best_match and fires on
                    # the common case here: indel-free point-mutation variants.
                    parent_aln, child_aln = child_seq, q_sub
                else:
                    parent_aln, child_aln = needleman_wunsch(
                        child_seq, q_sub, band=GAPPED_BAND
                    )
                ident = alignment_identity(parent_aln, child_aln)
                child_min_identity = self._min_identity_for_depth(child.depth)
                if ident >= child_min_identity and ident > best.identity:
                    best = _Match(
                        target_id=child_id,
                        target_seq=child_seq,
                        target_start=0,
                        target_end=len(child_seq),
                        q_start=match.q_start,
                        q_end=match.q_end,
                        identity=ident,
                        parent_aln=parent_aln,
                        child_aln=child_aln,
                    )
                    current_id = child_id
                    improved = True
                    break
            if not improved:
                break
        return best


def match_q_sub(m: "_Match") -> str:
    # Helper for descent step (avoids storing the query string on _Match).
    # The query substring lives at m.q_start..m.q_end on the original sequence,
    # but we already have its aligned-to-parent form (gaps stripped).
    return m.child_aln.replace("-", "")


def compress_fasta(
    input_paths: List[str],
    *,
    k: int = K_MER_SIZE,
    min_identity: float = MIN_IDENTITY,
    min_length: int = MIN_LENGTH,
    max_depth: int = MAX_DEPTH,
    progress=None,
) -> _Compressor:
    """Compress one or more FASTA files into an HS-CaBLASTP forest."""
    from cablastp.fasta import FastaReader

    comp = _Compressor(
        k=k, min_identity=min_identity, min_length=min_length, max_depth=max_depth,
    )
    total = 0
    for path in input_paths:
        with open(path, "rb") as fh:
            reader = FastaReader(fh)
            for rec in reader:
                seq = rec.residues.decode("ascii").upper()
                comp.compress_sequence(rec.name, seq)
                total += 1
                if progress and total % 100 == 0:
                    progress(total)
    if progress:
        progress(total, final=True)
    return comp
