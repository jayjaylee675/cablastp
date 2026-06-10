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

from typing import List, Tuple

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

# Gap penalty. The pure-Python kernel charges matrix[residue][gap] for every gap
# cell — the residue-vs-gap score (last MATRIX62 column = -4), NOT
# MATRIX62[gap][gap]=0 (the gap-vs-gap corner the DP never visits — scoring gaps
# as free there is the bug). parasail charges `open + (len-1)*extend` for a gap
# of `len`, so open == extend == 4 reproduces the reference's linear -4 gap.
_GAP_OPEN = -int(_M62[0][-1])    # residue-vs-gap penalty = 4
_GAP_EXTEND = _GAP_OPEN

# Restrict input residues to the parasail matrix alphabet. Unknown chars
# crash the kernel; map them to 'X' (index 20 in BLOSUM62 alphabet).
_VALID = set(_ALPHA)


def _clean(s: str) -> str:
    if all(c in _VALID for c in s):
        return s
    return "".join(c if c in _VALID else "X" for c in s)


def _restore(aligned: str, original: str) -> str:
    """Re-thread the ORIGINAL residues into parasail's gap structure.

    We align '_clean'-ed copies so parasail's matrix lookup can't crash on
    residues outside ALPHABET62 (U, O, J, '*'), but those substituted 'X's must
    never leak into the returned alignment: make_edit_script stores the aligned
    residues verbatim, so an 'X' would make the forest's reconstruction lossy.
    `aligned` is `original` with '-' gaps inserted, so its non-gap characters
    map 1:1, in order, back onto `original`.
    """
    out: List[str] = []
    k = 0
    for ch in aligned:
        if ch == "-":
            out.append("-")
        else:
            out.append(original[k])
            k += 1
    return "".join(out)


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
    # shape our downstream make_edit_script expects. When unknown residues had
    # to be cleansed to 'X', thread the originals back so the edit script (and
    # thus reconstruction) stays lossless.
    qa, rb = r.traceback.query, r.traceback.ref
    if sa != a:
        qa = _restore(qa, a)
    if sb != b:
        rb = _restore(rb, b)
    return qa, rb
