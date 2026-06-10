"""Sequence types: Sequence, CoarseSeq, OriginalSeq (port of seq.go)."""

from __future__ import annotations

import threading
from typing import Optional, TYPE_CHECKING

if TYPE_CHECKING:
    from cablastp.fasta import FastaSequence
    from cablastp.links import LinkToCompressed


def seq_identity(seq1: bytes, seq2: bytes) -> int:
    """Return the percentage identity (0..100) of two equal-length sequences."""
    if len(seq1) != len(seq2):
        raise ValueError(
            "Sequence identity requires that len(seq1) == len(seq2), "
            "but %d != %d." % (len(seq1), len(seq2))
        )
    if len(seq1) == 0 and len(seq2) == 0:
        return 0
    same = sum(1 for a, b in zip(seq1, seq2) if a == b)
    return (same * 100) // len(seq1)


def is_low_complexity(residues: bytes, offset: int, window: int) -> bool:
    """True if residues at `offset` lie inside a window of identical residues."""
    repeats = 1
    last = 0
    start = max(0, offset - window)
    end = min(len(residues), offset + window)
    for i in range(start, end):
        if residues[i] == last:
            repeats += 1
            if repeats >= window:
                return True
            continue
        last = residues[i]
        repeats = 1
    return False


def _repetitive(bs: bytes) -> bool:
    if len(bs) <= 1:
        return False
    first = bs[0]
    return all(b == first for b in bs[1:])


class Sequence:
    """Base sequence type, embedded by CoarseSeq and OriginalSeq."""

    __slots__ = ("name", "residues", "offset", "id")

    def __init__(self, id: int, name: str, residues: bytes):
        residues_str = residues.decode("ascii", errors="replace").upper().replace("*", "")
        self.name = name
        self.residues: bytes = residues_str.encode("ascii")
        self.offset = 0
        self.id = id

    def _new_sub_sequence(self, start: int, end: int) -> "Sequence":
        if start < 0 or start >= end or end > len(self.residues):
            raise ValueError(
                "Invalid sub sequence (%d, %d) for sequence with length %d."
                % (start, end, len(self.residues))
            )
        sub = Sequence(self.id, self.name, self.residues[start:end])
        # Mirror Go's newSubSequence, which sets the offset to `start` (relative
        # to this sequence) rather than accumulating self.offset. Current callers
        # only subsequence offset-0 sequences so the values coincide, but the Go
        # contract is a plain reset.
        sub.offset = start
        return sub

    def fasta_seq(self) -> "FastaSequence":
        from cablastp.fasta import FastaSequence
        return FastaSequence(name=self.name, residues=self.residues)

    def __len__(self) -> int:
        return len(self.residues)

    # `Len` exists in the Go API; expose a Python-friendly alias.
    def Len(self) -> int:  # noqa: N802
        return len(self.residues)

    def __str__(self) -> str:
        body = self.residues.decode("ascii")
        if self.offset == 0:
            return "> %s (%d)\n%s" % (self.name, self.id, body)
        return "> %s (%d) (%d, %d)\n%s" % (
            self.name, self.id, self.offset, len(self), body,
        )


class CoarseSeq(Sequence):
    """A reference sequence in the coarse database (with link list)."""

    __slots__ = ("links", "_link_lock")

    def __init__(self, id: int, name: str, residues: bytes):
        super().__init__(id, name, residues)
        self.links: Optional["LinkToCompressed"] = None
        self._link_lock = threading.RLock()

    def new_sub_sequence(self, start: int, end: int) -> "CoarseSeq":
        sub = self._new_sub_sequence(start, end)
        copy = CoarseSeq(self.id, self.name, sub.residues)
        copy.offset = sub.offset
        copy.links = None
        return copy

    def add_link(self, link: "LinkToCompressed") -> None:
        with self._link_lock:
            self._add_link_locked(link)

    def _add_link_locked(self, link: "LinkToCompressed") -> None:
        if self.links is None:
            self.links = link
            return
        node = self.links
        while node.next is not None:
            node = node.next
        node.next = link


class OriginalSeq(Sequence):
    """An input sequence read from a FASTA file."""

    def new_sub_sequence(self, start: int, end: int) -> "OriginalSeq":
        sub = self._new_sub_sequence(start, end)
        copy = OriginalSeq(self.id, self.name, sub.residues)
        copy.offset = sub.offset
        return copy


def new_fasta_seq(id: int, fasta_sequence: "FastaSequence") -> OriginalSeq:
    """Create an OriginalSeq from a FastaSequence."""
    return OriginalSeq(id, fasta_sequence.name, fasta_sequence.residues)


def new_fasta_coarse_seq(id: int, fasta_sequence: "FastaSequence") -> CoarseSeq:
    return CoarseSeq(id, fasta_sequence.name, fasta_sequence.residues)
