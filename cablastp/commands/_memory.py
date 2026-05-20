"""Per-worker memory arena (port of cmd/cablastp-compress/memory.go)."""

from __future__ import annotations

from typing import List, Tuple

MEM_SEQ_SIZE = 10000
DYNAMIC_TABLE_SIZE = MEM_SEQ_SIZE * MEM_SEQ_SIZE
NUM_SEEDS = 100


class Memory:
    """Reusable buffers used by alignment / seed lookup hot paths."""

    __slots__ = ("table", "ref", "org", "seeds")

    def __init__(self):
        self.table: List[int] = [0] * (MEM_SEQ_SIZE * MEM_SEQ_SIZE)
        self.ref: bytearray = bytearray(MEM_SEQ_SIZE)
        self.org: bytearray = bytearray(MEM_SEQ_SIZE)
        self.seeds: List[Tuple[int, int]] = []


def new_memory() -> Memory:
    return Memory()
