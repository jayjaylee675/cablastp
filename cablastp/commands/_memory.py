"""Per-worker memory arena (port of cmd/cablastp-compress/memory.go)."""

from __future__ import annotations

from typing import List, Optional, Tuple

MEM_SEQ_SIZE = 10000
DYNAMIC_TABLE_SIZE = MEM_SEQ_SIZE * MEM_SEQ_SIZE
NUM_SEEDS = 100


class Memory:
    """Reusable buffers used by alignment / seed lookup hot paths."""

    __slots__ = ("_table", "ref", "org", "seeds")

    def __init__(self):
        # The NW DP table is ~800 MB and is touched ONLY by the pure-Python NW
        # fallback. Under the default parasail backend it is never used, so
        # allocate it lazily on first access instead of eagerly per worker ×
        # os.cpu_count() at startup (which cost ~800 MB/worker for nothing).
        self._table: Optional[List[int]] = None
        self.ref: bytearray = bytearray(MEM_SEQ_SIZE)
        self.org: bytearray = bytearray(MEM_SEQ_SIZE)
        self.seeds: List[Tuple[int, int]] = []

    @property
    def table(self) -> List[int]:
        if self._table is None:
            self._table = [0] * DYNAMIC_TABLE_SIZE
        return self._table


def new_memory() -> Memory:
    return Memory()
