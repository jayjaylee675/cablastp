"""Pytest equivalents of cmd/cablastp-compress/compress_test.go."""

from __future__ import annotations

import pytest

from cablastp.commands._align import align_ungapped
from cablastp.commands._compress_core import extend_match, skip_low_complexity
from cablastp.commands._memory import new_memory
from cablastp.commands._nw import nw_align
from cablastp.blosum import ALPHABET62, MATRIX62


_RES_IDX = {c: i for i, c in enumerate(ALPHABET62)}


def _sub_score(x: int, y: int) -> int:
    """BLOSUM62 score for two residue bytes; unknown residues fall back to the
    kernels' -4 (the matrix gap row/col), so this matches their scoring."""
    ix = _RES_IDX.get(chr(x), -1)
    iy = _RES_IDX.get(chr(y), -1)
    if ix < 0 or iy < 0:
        return -4
    return MATRIX62[ix][iy]


def _aln_score(a: bytes, b: bytes) -> int:
    """Score a gap-padded alignment under BLOSUM62 with a linear -4 gap."""
    gap = ord("-")
    s = 0
    for x, y in zip(a, b):
        s += -4 if (x == gap or y == gap) else _sub_score(x, y)
    return s


def _nw_opt_score(a: bytes, b: bytes) -> int:
    """Independent full (unbanded) Needleman-Wunsch optimum, linear -4 gap.

    A backend-agnostic ground truth for the optimality assertion below.
    """
    n, m = len(a), len(b)
    dp = [[0] * (m + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
        dp[i][0] = -4 * i
    for j in range(1, m + 1):
        dp[0][j] = -4 * j
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            dp[i][j] = max(
                dp[i - 1][j - 1] + _sub_score(a[i - 1], b[j - 1]),
                dp[i - 1][j] - 4,
                dp[i][j - 1] - 4,
            )
    return dp[n][m]


@pytest.mark.parametrize(
    "seq,window_size,region_size,expected_skipped,expected_leftover",
    [
        (b"ABCDDDDDDDDDDDDDDDDDDXYZ", 10, 5, 21, b"XYZ"),
        (b"DDDDDDABCDEF", 10, 5, 6, b"ABCDEF"),
        (b"DDDDDDABCDEFFFFFFFFFFXYZXYZ", 10, 5, 6, b"ABCDEFFFFFFFFFFXYZXYZ"),
        (b"ABCDEFFFFFFFFFFFFFFFXYZXYZ", 10, 5, 20, b"XYZXYZ"),
    ],
)
def test_skip_low_complexity(
    seq, window_size, region_size, expected_skipped, expected_leftover,
):
    skipped = skip_low_complexity(seq, window_size, region_size)
    assert skipped == expected_skipped
    assert seq[skipped:] == expected_leftover


@pytest.mark.parametrize(
    "seq1,seq2",
    [
        (b"ABCD", b"ABCD"),
        (b"PPPGHIKLMNPQR", b"GAAAHIKLMN"),
        (b"GHIKLMNPQRSTVW", b"GAAAHIKLMNPQRSTVW"),
        (b"X" * 62, b"X" * 62),
        (b"NNNNNNNN", b"NNNNNNNN"),
        (b"N" * 37, b"N" * 37),
    ],
)
def test_needleman_wunsch(seq1, seq2):
    # nw_align routes through parasail when installed and the pure-Python
    # fallback otherwise. Rather than pin backend-specific gap-padded strings
    # (which made this test RED whenever parasail was present and never
    # exercised the production backend), assert the properties any correct
    # global alignment must satisfy.
    mem = new_memory()
    a, b = nw_align(seq1, seq2, mem)
    gap = ord("-")

    # Equal-length, gap-padded output.
    assert len(a) == len(b)
    # Lossless: dropping the gaps recovers each input verbatim and in order.
    # (This is exactly what the parasail 'X'-rewrite bug used to violate.)
    assert a.replace(b"-", b"") == seq1
    assert b.replace(b"-", b"") == seq2
    # No column aligns a gap against a gap.
    assert all(not (x == gap and y == gap) for x, y in zip(a, b))
    # Identical inputs must align on the diagonal with no gaps.
    if seq1 == seq2:
        assert a == seq1 and b == seq2
    # The alignment is optimal under BLOSUM62 with the kernels' linear -4 gap.
    assert _aln_score(a, b) == _nw_opt_score(seq1, seq2)


@pytest.mark.parametrize(
    "rseq,oseq,rmseq,omseq",
    [
        (
            b"ABCDEFGHIKLMNPQR",
            b"ABCDEFGHIKLMNPQR",
            b"ABCDEFGHIKLMNPQR",
            b"ABCDEFGHIKLMNPQR",
        ),
        (
            b"ABCDEFGHIKLMNPQRSTVW",
            b"ABCDEFGAAAHIKLMNPQRSTVW",
            b"ABCDEFGHIKLMNPQRSTVW",
            b"ABCDEFGAAAHIKLMNPQRSTVW",
        ),
        (
            b"ABCDEFGHIKLMNPQRSTVW",
            b"ABCDEFGAAAHIKLMNPQRSTBBBBBBBBBBBBBBBBBBBVW",
            b"ABCDEF",
            b"ABCDEF",
        ),
    ],
)
def test_extend_match(rseq, oseq, rmseq, omseq):
    flag_match_kmer_size = 3
    flag_ungapped_window_size = 10
    flag_ext_seq_id_threshold = 50
    flag_gapped_window_size = 25

    mem = new_memory()
    cor_match, org_match = extend_match(
        rseq, oseq,
        flag_gapped_window_size, flag_ungapped_window_size,
        flag_match_kmer_size, flag_ext_seq_id_threshold,
        mem,
    )
    assert cor_match == rmseq
    assert org_match == omseq


@pytest.mark.parametrize(
    "rseq,oseq,answer",
    [
        (b"A", b"A", 0),
        (b"AB", b"AB", 0),
        (b"ABC", b"ABC", 3),
        (b"ABCD", b"ABCD", 3),
        (b"ABCYEFG", b"ABCZEFG", 3),
        (b"ABCYEFGH", b"ABCZEFGH", 8),
        (b"ABCDEFGHIJKLMNOP", b"ABCDEFGHIJKLMNOP", 15),
        (b"ABCDEF", b"ABC", 3),
        (b"ABC", b"ABCDEF", 3),
        (b"ABCDEFGHIKLMNPQR", b"ABCDEFGHIKLMNPQR", 15),
    ],
)
def test_ungapped_extension(rseq, oseq, answer):
    flag_match_kmer_size = 3
    flag_ungapped_window_size = 10
    flag_ext_seq_id_threshold = 50

    result = align_ungapped(
        rseq, oseq,
        flag_ungapped_window_size, flag_match_kmer_size,
        flag_ext_seq_id_threshold,
    )
    assert result == answer
