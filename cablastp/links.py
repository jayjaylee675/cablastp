"""Bridges between coarse and compressed sequences (port of link_to_*.go)."""

from __future__ import annotations

from typing import Optional, Tuple

from cablastp.seqdiff import EditScript, new_edit_script


class LinkToCoarse:
    """A piece of a compressed original sequence pointing into the coarse DB.

    `diff` is the textual edit script that, when applied to the indicated coarse
    sub-sequence, reproduces this slice of the original sequence exactly.
    """

    __slots__ = ("diff", "coarse_seq_id", "coarse_start", "coarse_end")

    def __init__(
        self,
        diff: str,
        coarse_seq_id: int,
        coarse_start: int,
        coarse_end: int,
    ):
        self.diff = diff
        self.coarse_seq_id = int(coarse_seq_id)
        self.coarse_start = int(coarse_start) & 0xFFFF
        self.coarse_end = int(coarse_end) & 0xFFFF

    def __str__(self) -> str:
        return (
            "reference sequence id: %d, reference range: (%d, %d)\n%s"
            % (self.coarse_seq_id, self.coarse_start, self.coarse_end, self.diff)
        )


def new_link_to_coarse(
    coarse_seq_id: int,
    coarse_start: int,
    coarse_end: int,
    alignment: Tuple[bytes, bytes],
) -> LinkToCoarse:
    return LinkToCoarse(
        diff=str(new_edit_script(alignment)),
        coarse_seq_id=coarse_seq_id,
        coarse_start=coarse_start,
        coarse_end=coarse_end,
    )


def new_link_to_coarse_no_diff(
    coarse_seq_id: int,
    coarse_start: int,
    coarse_end: int,
) -> LinkToCoarse:
    return LinkToCoarse(
        diff="",
        coarse_seq_id=coarse_seq_id,
        coarse_start=coarse_start,
        coarse_end=coarse_end,
    )


class LinkToCompressed:
    """A node in the linked list of links from a coarse sequence to compressed
    original sequences."""

    __slots__ = ("org_seq_id", "coarse_start", "coarse_end", "next")

    def __init__(self, org_seq_id: int, coarse_start: int, coarse_end: int):
        self.org_seq_id = int(org_seq_id) & 0xFFFFFFFF
        self.coarse_start = int(coarse_start) & 0xFFFF
        self.coarse_end = int(coarse_end) & 0xFFFF
        self.next: Optional["LinkToCompressed"] = None

    def __str__(self) -> str:
        return "original sequence id: %d, coarse range: (%d, %d)" % (
            self.org_seq_id, self.coarse_start, self.coarse_end,
        )
