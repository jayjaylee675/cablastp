"""Alignment helpers used during compression (port of align.go)."""

from __future__ import annotations

from cablastp.seq import seq_identity


def align_len(seq: bytes) -> int:
    """Number of non-gap residues in `seq`."""
    return sum(1 for r in seq if r != ord("-"))


def align_ungapped(
    rseq: bytes,
    oseq: bytes,
    window_size: int,
    kmer_size: int,
    id_threshold: int,
) -> int:
    """Greedy ungapped extension using exact K-mer matches.

    Mirrors `alignUngapped` from align.go: returns the number of residues
    consumed by successful K-mer matches, or 0 if none.
    """
    length = 0
    scanned = 0
    successive = 0
    try_next_window = True
    while try_next_window:
        try_next_window = False
        for _ in range(window_size):
            if scanned >= len(rseq) or scanned >= len(oseq):
                break

            if rseq[scanned] == oseq[scanned]:
                successive += 1
            else:
                successive = 0
            scanned += 1

            if successive == kmer_size:
                if (scanned - kmer_size) - length > 0:
                    id = seq_identity(
                        rseq[length:scanned - kmer_size],
                        oseq[length:scanned - kmer_size],
                    )
                    if id < id_threshold:
                        successive -= 1
                        continue
                length = scanned
                successive = 0
                try_next_window = True
                break
    return length
