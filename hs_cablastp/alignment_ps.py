"""Parasail-backed Needleman-Wunsch — prototype.

Drop-in for hs_cablastp.alignment.needleman_wunsch using parasail's SIMD NW.
Same signature, same scoring (BLOSUM62 + the existing kernel's effective
zero gap penalty — see KNOWN_ISSUES.md / discussion: MATRIX62[gap][gap]=0).

Caveats vs the pure-Python kernel:
- With gap_open = gap_extend = 0, many NW alignments tie on score. Parasail
  resolves ties differently than the pure-Python kernel, so the returned
  gap-padded strings can differ even when the score (and downstream identity)
  is identical.
- Parasail's banded NW does not provide trace, so we always use full NW with
  trace and accept the (n*m) work. For the candidate-window sizes used in
  hs-cablastp this is still well within the SIMD kernel's sweet spot.
"""

from __future__ import annotations

from typing import Tuple

import parasail

import cablastp.blosum as _bm

# Build a parasail Matrix whose substitution scores mirror cablastp's
# MATRIX62 exactly. We index only the 23 residue columns; parasail handles
# gaps via its own open/extend parameters.
_ALPHA = _bm.ALPHABET62          # 'ABCDEFGHIKLMNPQRSTVWXYZ'
_M62 = _bm.MATRIX62
_PS_MATRIX = parasail.matrix_create(_ALPHA, 0, 0)
for _i in range(len(_ALPHA)):
    for _j in range(len(_ALPHA)):
        _PS_MATRIX.set_value(_i, _j, int(_M62[_i][_j]))

# Match the pure-Python kernel: it uses MATRIX62[gap][gap] = 0 as the per-cell
# gap penalty (effectively no gap penalty). Replicate so the score function
# matches and the downstream identity computation behaves the same.
_GAP_OPEN = -int(_M62[-1][-1])   # 0 today; positive when the bug is fixed
_GAP_EXTEND = _GAP_OPEN

# Restrict input residues to the parasail matrix alphabet. Unknown chars
# crash the kernel; map them to 'X' (index 20 in BLOSUM62 alphabet).
_VALID = set(_ALPHA)


def _clean(s: str) -> str:
    if all(c in _VALID for c in s):
        return s
    return "".join(c if c in _VALID else "X" for c in s)


def needleman_wunsch_ps(a: str, b: str, band: int = 0) -> Tuple[str, str]:
    n, m = len(a), len(b)
    if n == 0:
        return "-" * m, b
    if m == 0:
        return a, "-" * n
    sa = _clean(a)
    sb = _clean(b)
    r = parasail.nw_trace_scan_16(sa, sb, _GAP_OPEN, _GAP_EXTEND, _PS_MATRIX)
    # parasail's `query` is the first arg's aligned form, `ref` is the second's.
    # Both are equal-length strings padded with '-' for gaps — exactly the
    # shape our downstream make_edit_script expects.
    return r.traceback.query, r.traceback.ref
