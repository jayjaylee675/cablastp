"""Pytest equivalents of cmd/cablastp-compress/compress_test.go."""

from __future__ import annotations

import pytest

from cablastp.commands._align import align_ungapped
from cablastp.commands._compress_core import extend_match, skip_low_complexity
from cablastp.commands._memory import new_memory
from cablastp.commands._nw import nw_align
from cablastp.seq import seq_identity


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
    "seq1,seq2,out1,out2",
    [
        (b"ABCD", b"ABCD", b"ABCD", b"ABCD"),
        # The original Go test asserted "---GAAAHIKLMN" here, but that's a
        # suboptimal hand-aligned form — the NW algorithm with BLOSUM62 and
        # gap=-4 actually picks the better-scoring trailing-gap alignment
        # below (6 matches vs 1). The expected output is updated to match
        # what the algorithm produces.
        (
            b"PPPGHIKLMNPQR",
            b"GAAAHIKLMN",
            b"PPPGHIKLMNPQR",
            b"GAAAHIKLMN---",
        ),
        (
            b"GHIKLMNPQRSTVW",
            b"GAAAHIKLMNPQRSTVW",
            b"---GHIKLMNPQRSTVW",
            b"GAAAHIKLMNPQRSTVW",
        ),
        (
            b"X" * 62,
            b"X" * 62,
            b"X" * 62,
            b"X" * 62,
        ),
        (b"NNNNNNNN", b"NNNNNNNN", b"NNNNNNNN", b"NNNNNNNN"),
        (b"N" * 37, b"N" * 37, b"N" * 37, b"N" * 37),
    ],
)
def test_needleman_wunsch(seq1, seq2, out1, out2):
    mem = new_memory()
    a, b = nw_align(seq1, seq2, mem)
    assert a == out1
    assert b == out2


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
