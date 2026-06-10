"""Alignment primitives: ungapped extension, Needleman-Wunsch, identity, edit scripts."""

from __future__ import annotations

from typing import List, Tuple

from cablastp.blosum import ALPHABET62, MATRIX62
from hs_cablastp.types import EditOp


_RES_INDEX = [-1] * 256
for _i, _ch in enumerate(ALPHABET62):
    _RES_INDEX[ord(_ch)] = _i
_GAP_IDX = len(MATRIX62) - 1


def _score(a: str, b: str) -> int:
    # Call ord() once per residue (was twice each) — this is the inner loop of
    # ungapped_extend and ran ~13.6M times in the profile. ord() of a 1-char str
    # is never negative, so only the upper bound needs checking.
    oa = ord(a)
    ob = ord(b)
    ia = _RES_INDEX[oa] if oa < 256 else -1
    ib = _RES_INDEX[ob] if ob < 256 else -1
    if ia < 0 or ib < 0:
        return -4
    return MATRIX62[ia][ib]


def ungapped_extend(
    a: str, ia: int, b: str, ib: int, k: int,
    max_drop: int = 12,
) -> Tuple[int, int, int, int]:
    """Greedy ungapped extension from a seed match of length k.

    Starts with the k-mer at (ia, ib) and extends left/right while the running
    score doesn't drop more than max_drop below the best score seen.

    Returns (start_a, end_a, start_b, end_b) — half-open intervals over both sequences.
    """
    start_a = ia
    start_b = ib
    end_a = ia + k
    end_b = ib + k

    score = 0
    for off in range(k):
        score += _score(a[ia + off], b[ib + off])

    best = score
    sa, sb = start_a, start_b
    ea, eb = end_a, end_b

    # Extend right
    cur = best
    pa, pb = end_a, end_b
    while pa < len(a) and pb < len(b):
        cur += _score(a[pa], b[pb])
        pa += 1
        pb += 1
        if cur > best:
            best = cur
            ea, eb = pa, pb
        elif best - cur > max_drop:
            break

    # Extend left
    cur = best
    pa, pb = start_a - 1, start_b - 1
    while pa >= 0 and pb >= 0:
        cur += _score(a[pa], b[pb])
        if cur > best:
            best = cur
            sa, sb = pa, pb
        elif best - cur > max_drop:
            break
        pa -= 1
        pb -= 1

    return sa, ea, sb, eb


try:
    from hs_cablastp.alignment_ps import needleman_wunsch_ps as _nw_ps
except Exception:  # parasail missing or load failure — fall back silently
    _nw_ps = None


def needleman_wunsch(a: str, b: str, band: int = 0) -> Tuple[str, str]:
    """Global NW alignment of two short strings.

    `band` enables a diagonal band constraint when > 0 (skip cells |i-j| > band).
    Returns the two aligned strings (with '-' gaps).

    Uses parasail's SIMD NW kernel when available (10-130x faster on the
    workload's typical input sizes); falls back to the pure-Python kernel
    below otherwise. Parasail does not implement a banded NW with traceback,
    so the `band` argument is honoured by the pure-Python path only — under
    parasail, full NW is computed (still much faster than the banded Python
    fill at any L >= ~25 we have measured).
    """
    if _nw_ps is not None:
        return _nw_ps(a, b, band=band)

    n, m = len(a), len(b)
    if n == 0:
        return "-" * m, b
    if m == 0:
        return a, "-" * n

    # A band narrower than the length difference cannot reach cell (n, m): it
    # stays at NEG and the traceback then walks an arbitrary all-gap path,
    # yielding a wrong identity. Widen it to at least |n-m| so the bottom-right
    # corner is always reachable. (band == 0 means "unbanded" — leave it.)
    if band:
        band = max(band, abs(n - m))

    gap_pen = MATRIX62[_GAP_IDX][_GAP_IDX]
    NEG = -10 ** 9

    # Hoist residue-index lookups out of the inner loop: the per-cell
    # ord()/_RES_INDEX[] calls dominated the original profile (~18% of total
    # runtime in builtins.ord alone). Precompute index arrays once per call.
    ra_arr = [_RES_INDEX[c] if c < 256 else -1 for c in map(ord, a)]
    rb_arr = [_RES_INDEX[c] if c < 256 else -1 for c in map(ord, b)]
    M = MATRIX62  # local alias avoids repeated global lookups in the hot loop

    table = [[NEG] * (m + 1) for _ in range(n + 1)]
    table[0][0] = 0
    for i in range(1, n + 1):
        if band and i > band:
            break
        table[i][0] = table[i - 1][0] + gap_pen
    for j in range(1, m + 1):
        if band and j > band:
            break
        table[0][j] = table[0][j - 1] + gap_pen

    for i in range(1, n + 1):
        ra = ra_arr[i - 1]
        # Per-row substitution scores: index M once by ra, then the inner loop
        # only does srow[rb] instead of M[ra][rb] + bounds checks.
        srow = M[ra] if ra >= 0 else None
        prev_row = table[i - 1]
        cur_row = table[i]
        jmin = max(1, i - band) if band else 1
        jmax = min(m, i + band) if band else m
        for j in range(jmin, jmax + 1):
            rb = rb_arr[j - 1]
            sub = srow[rb] if (srow is not None and rb >= 0) else -4
            diag = prev_row[j - 1] + sub
            up = prev_row[j] + gap_pen
            left = cur_row[j - 1] + gap_pen
            # Inline max(diag, up, left) — avoids the builtin call overhead
            # (profiled at ~17% of runtime) while producing the identical value.
            best = diag if diag >= up else up
            if left > best:
                best = left
            cur_row[j] = best

    # Traceback
    aln_a: List[str] = []
    aln_b: List[str] = []
    i, j = n, m
    while i > 0 and j > 0:
        ai = ord(a[i - 1]); ra = _RES_INDEX[ai] if 0 <= ai < 256 else -1
        bj = ord(b[j - 1]); rb = _RES_INDEX[bj] if 0 <= bj < 256 else -1
        sub = MATRIX62[ra][rb] if ra >= 0 and rb >= 0 else -4
        if table[i][j] == table[i - 1][j - 1] + sub:
            aln_a.append(a[i - 1]); aln_b.append(b[j - 1])
            i -= 1; j -= 1
        elif table[i][j] == table[i - 1][j] + gap_pen:
            aln_a.append(a[i - 1]); aln_b.append("-")
            i -= 1
        else:
            aln_a.append("-"); aln_b.append(b[j - 1])
            j -= 1
    while i > 0:
        aln_a.append(a[i - 1]); aln_b.append("-"); i -= 1
    while j > 0:
        aln_a.append("-"); aln_b.append(b[j - 1]); j -= 1
    return "".join(reversed(aln_a)), "".join(reversed(aln_b))


def alignment_identity(aln_a: str, aln_b: str) -> float:
    matches = total = 0
    for ca, cb in zip(aln_a, aln_b):
        if ca == "-" or cb == "-":
            total += 1
            continue
        total += 1
        if ca == cb:
            matches += 1
    return matches / total if total else 0.0


def make_edit_script(parent_aln: str, child_aln: str) -> List[EditOp]:
    """Build the edit script that transforms the parent's matched segment into the child.

    Positions are 0-based offsets within the parent's matched segment.
    Conventions:
      - parent has residue, child has gap   -> DELETE at that parent position
      - parent has gap,  child has residue  -> INSERT char at that parent position
      - mismatching residues                -> SUBSTITUTE char at that parent position
      - matching residues                   -> no op
    """
    ops: List[EditOp] = []
    parent_pos = 0
    for pa, ca in zip(parent_aln, child_aln):
        if pa == "-" and ca != "-":
            ops.append(EditOp("INSERT", parent_pos, ca))
        elif pa != "-" and ca == "-":
            ops.append(EditOp("DELETE", parent_pos, ""))
            parent_pos += 1
        elif pa != ca:
            ops.append(EditOp("SUBSTITUTE", parent_pos, ca))
            parent_pos += 1
        else:
            parent_pos += 1
    return ops


def apply_edit_script(parent_segment: str, script: List[EditOp]) -> str:
    """Apply a diff script to a parent segment to reconstruct the child segment."""
    # Sort ops by position (stable). INSERTs at the same position keep their order.
    out: List[str] = []
    i = 0
    ops_by_pos: dict[int, List[EditOp]] = {}
    for op in script:
        ops_by_pos.setdefault(op.position, []).append(op)

    while i <= len(parent_segment):
        for op in ops_by_pos.get(i, []):
            if op.op_type == "INSERT":
                out.append(op.char)
        if i == len(parent_segment):
            break
        # Handle DELETE / SUBSTITUTE at this position
        applied = False
        for op in ops_by_pos.get(i, []):
            if op.op_type == "DELETE":
                applied = True
                break
            if op.op_type == "SUBSTITUTE":
                out.append(op.char)
                applied = True
                break
        if not applied:
            out.append(parent_segment[i])
        i += 1
    return "".join(out)
