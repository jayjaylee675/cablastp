"""Needleman-Wunsch sequence alignment (port of nw.go)."""

from __future__ import annotations

from typing import Tuple

from cablastp.blosum import ALPHABET62, MATRIX62
from cablastp.commands._memory import DYNAMIC_TABLE_SIZE, Memory

# Translate ASCII residue characters to BLOSUM62 indices.
_RES_TRANS = [0] * 256
for _idx, _ch in enumerate(ALPHABET62):
    _RES_TRANS[ord(_ch)] = _idx


def nw_align(rseq: bytes, oseq: bytes, mem: Memory) -> Tuple[bytes, bytes]:
    """Needleman-Wunsch global alignment, constrained for long sequences."""
    gap = len(MATRIX62) - 1
    r = len(rseq) + 1
    c = len(oseq) + 1

    constrained = True
    constraint = max(r, c) // 4
    if r <= 11 or c <= 11:
        constrained = False

    if r * c > DYNAMIC_TABLE_SIZE:
        table = [0] * (r * c)
    else:
        table = mem.table
        for i in range(r * c):
            table[i] = 0

    matrix = MATRIX62
    res_trans = _RES_TRANS
    gap_char = ord("-")

    for i in range(1, r):
        i2 = (i - 1) * c
        i3 = i * c
        rval = res_trans[rseq[i - 1]]
        for j in range(1, c):
            if constrained and (abs(i - j) > constraint):
                continue
            oval = res_trans[oseq[j - 1]]

            off = i2 + (j - 1)
            sdiag = table[off] + matrix[rval][oval]
            sup = table[off + 1] + matrix[rval][gap]
            sleft = table[off + c] + matrix[gap][oval]
            if sdiag > sup and sdiag > sleft:
                table[i3 + j] = sdiag
            elif sup > sleft:
                table[i3 + j] = sup
            else:
                table[i3 + j] = sleft

    ref_aln = bytearray()
    org_aln = bytearray()

    i = r - 1
    j = c - 1
    while i > 0 and j > 0:
        rval = res_trans[rseq[i - 1]]
        oval = res_trans[oseq[j - 1]]
        sdiag = table[(i - 1) * c + (j - 1)] + matrix[rval][oval]
        sup = table[(i - 1) * c + j] + matrix[gap][oval]
        sleft = table[i * c + (j - 1)] + matrix[rval][gap]
        if sdiag > sup and sdiag > sleft:
            i -= 1
            j -= 1
            ref_aln.append(rseq[i])
            org_aln.append(oseq[j])
        elif sup > sleft:
            i -= 1
            ref_aln.append(rseq[i])
            org_aln.append(gap_char)
        else:
            j -= 1
            ref_aln.append(gap_char)
            org_aln.append(oseq[j])

    while i > 0:
        i -= 1
        ref_aln.append(rseq[i])
        org_aln.append(gap_char)
    while j > 0:
        j -= 1
        ref_aln.append(gap_char)
        org_aln.append(oseq[j])

    ref_aln.reverse()
    org_aln.reverse()
    return bytes(ref_aln), bytes(org_aln)
