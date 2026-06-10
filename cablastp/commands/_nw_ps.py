"""Parasail-backed Needleman-Wunsch for cablastp.

Same signature and contract as cablastp.commands._nw.nw_align:
  - input: bytes (ASCII residues)
  - output: (ref_aln, org_aln) bytes, equal length, '-' (ord 45) for gaps

The `mem` arena from the pure-Python kernel is unused: parasail allocates its
own scratch internally. We accept it to keep the call site identical.

Approximation: the pure-Python/Go kernel constrains long alignments to a
diagonal band of width max(len)/4, but parasail's banded NW provides no
traceback, so this backend always runs the full (n*m) NW. The optimal global
alignment can therefore differ from the reference when the best path strays
outside that band — rare for the near-diagonal homologies cablastp aligns.
"""

from __future__ import annotations

from typing import Tuple

import parasail

import cablastp.blosum as _bm

# Build a parasail Matrix whose substitution scores mirror cablastp's MATRIX62
# (the BLOSUM62 used by the existing kernel) one-for-one. We index only the
# 23-residue alphabet — parasail handles gaps via its own open/extend args.
_ALPHA_STR = _bm.ALPHABET62                  # 'ABCDEFGHIKLMNPQRSTVWXYZ'
_ALPHA_BYTES = _ALPHA_STR.encode("ascii")
_M62 = _bm.MATRIX62
_PS_MATRIX = parasail.matrix_create(_ALPHA_STR, 0, 0)
for _i in range(len(_ALPHA_STR)):
    for _j in range(len(_ALPHA_STR)):
        _PS_MATRIX.set_value(_i, _j, int(_M62[_i][_j]))

# Gap penalty. The pure-Python/Go kernel charges matrix[residue][gap] for every
# gap cell — the residue-vs-gap score in the last MATRIX62 column, which is -4.
# (NOT MATRIX62[gap][gap]=0, the gap-vs-gap *corner* the DP never visits — that
# bug made parasail score gaps as free and diverge from the reference.) parasail
# uses positive open/extend penalties and charges `open + (len-1)*extend` for a
# gap of `len`, so open == extend == 4 reproduces the reference's linear -4 gap.
_GAP_OPEN = -int(_M62[0][-1])    # residue-vs-gap penalty = 4
_GAP_EXTEND = _GAP_OPEN

# Map any unknown residue byte to 'A' (BLOSUM62 index 0) so parasail's lookup
# can't crash AND scores it exactly as the Go/pure-Python kernel does: both fold
# out-of-alphabet residues onto index 0 (_RES_TRANS in _nw.py defaults to 0), not
# 'X'. An unknown vs 'A' then scores MATRIX62[0][...] in every backend instead of
# diverging on the 'X' row. The ORIGINAL bytes are threaded back via _restore, so
# this only affects scoring, never the residues emitted into the diff. Build the
# 256-byte table once.
_UNKNOWN = _ALPHA_BYTES[0]   # ord('A'), the index-0 residue
_TRANS = bytearray(256)
for _b in range(256):
    _TRANS[_b] = _b if _b in _ALPHA_BYTES else _UNKNOWN
_TRANS = bytes(_TRANS)


def _restore(aligned: str, original: bytes) -> bytes:
    """Re-thread the ORIGINAL residue bytes into parasail's gap structure.

    We align copies whose out-of-alphabet residues were folded to 'A' (so
    parasail's matrix lookup can't crash on residues outside ALPHABET62), but
    those substituted residues must never leak into the returned alignment:
    cablastp stores the aligned residues verbatim in the diff/edit script, so a
    folded residue here would make decompression lossy. The aligned string is the
    input with '-' gaps inserted, so its non-gap characters map 1:1, in order,
    back onto `original`.
    """
    out = bytearray()
    k = 0
    gap = ord("-")
    for ch in aligned.encode("ascii"):
        if ch == gap:
            out.append(gap)
        else:
            out.append(original[k])
            k += 1
    return bytes(out)


def nw_align_ps(rseq: bytes, oseq: bytes, mem) -> Tuple[bytes, bytes]:
    if not rseq:
        return b"-" * len(oseq), bytes(oseq)
    if not oseq:
        return bytes(rseq), b"-" * len(rseq)
    # Align copies whose out-of-alphabet residues are folded to 'A' (so
    # parasail's lookup can't crash) but thread the ORIGINAL bytes back in.
    rs = rseq.translate(_TRANS).decode("ascii")
    os = oseq.translate(_TRANS).decode("ascii")
    r = parasail.nw_trace_scan_16(rs, os, _GAP_OPEN, _GAP_EXTEND, _PS_MATRIX)
    # parasail returns the aligned forms of the inputs in order:
    #   .query == aligned form of the first arg (= rseq)
    #   .ref   == aligned form of the second arg (= oseq)
    return _restore(r.traceback.query, rseq), _restore(r.traceback.ref, oseq)
