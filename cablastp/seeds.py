"""K-mer seeds table (port of seeds.go)."""

from __future__ import annotations

import threading
from typing import List, Optional, Tuple

from cablastp.blosum import ALPHABET62
from cablastp.seq import is_low_complexity, CoarseSeq

SEED_ALPHA_SIZE: int = len(ALPHABET62)
SEED_ALPHA_NUMS: List[int] = [-1] * 26
REVERSE_SEED_ALPHA_NUMS: bytearray = bytearray(26)


def _populate_alpha_tables() -> None:
    amino_val = 0
    for i in range(26):
        letter = ord("A") + i
        if chr(letter) in ALPHABET62:
            SEED_ALPHA_NUMS[i] = amino_val
            REVERSE_SEED_ALPHA_NUMS[amino_val] = letter
            amino_val += 1
        else:
            SEED_ALPHA_NUMS[i] = -1


_populate_alpha_tables()


def amino_value(letter: int) -> int:
    """Return the base-N index of `letter` (ASCII value of an upper-case residue)."""
    val = SEED_ALPHA_NUMS[letter - ord("A")]
    if val == -1:
        raise ValueError("Invalid amino acid letter: %s" % chr(letter))
    return val


class SeedLoc:
    """A node in the linked list of (sequence index, residue index) seed locations."""

    __slots__ = ("seq_ind", "res_ind", "next")

    def __init__(self, seq_ind: int, res_ind: int):
        self.seq_ind = int(seq_ind) & 0xFFFFFFFF
        self.res_ind = int(res_ind) & 0xFFFF
        self.next: Optional["SeedLoc"] = None


class Seeds:
    """A table of seed-location lists indexed by hashed K-mer."""

    def __init__(self, seed_size: int, low_complexity_window: int):
        self.seed_size = seed_size
        self.low_complexity_window = low_complexity_window
        self.lock = threading.RLock()
        self.powers: List[int] = []
        p = 1
        for _ in range(seed_size + 1):
            self.powers.append(p)
            p *= SEED_ALPHA_SIZE
        self.locs: List[Optional[SeedLoc]] = [None] * self.powers[seed_size]
        self.num_seeds: int = 0

    def NumSeeds(self) -> int:  # noqa: N802
        with self.lock:
            return self.num_seeds

    def num_seeds_count(self) -> int:
        return self.NumSeeds()

    def MaybeWipe(self, seed_table_size_gb: float) -> None:  # noqa: N802
        seed_loc_size = 16
        max_bytes = seed_table_size_gb * 1024.0 * 1024.0 * 1024.0
        max_seeds = int(max_bytes) // seed_loc_size
        if self.NumSeeds() >= max_seeds:
            print("Blowing away seeds table...")
            with self.lock:
                for i in range(len(self.locs)):
                    self.locs[i] = None
                self.num_seeds = 0

    maybe_wipe = MaybeWipe

    def Add(self, coarse_seq_index: int, cor_seq: CoarseSeq) -> None:  # noqa: N802
        with self.lock:
            residues = cor_seq.residues
            limit = len(residues) - self.seed_size
            for i in range(0, limit):
                if is_low_complexity(residues, i, self.low_complexity_window):
                    continue
                kmer = residues[i:i + self.seed_size]
                kmer_index = self.hash_kmer(kmer)
                loc = SeedLoc(coarse_seq_index, i)
                self.num_seeds += 1
                if self.locs[kmer_index] is None:
                    self.locs[kmer_index] = loc
                else:
                    node = self.locs[kmer_index]
                    while node.next is not None:
                        node = node.next
                    node.next = loc

    add = Add

    def Lookup(self, kmer: bytes, mem: List[Tuple[int, int]]) -> List[Tuple[int, int]]:  # noqa: N802
        """Return all (seq_ind, res_ind) seed locations for `kmer`.

        `mem` is reused as backing storage to avoid allocations.
        """
        with self.lock:
            head = self.locs[self.hash_kmer(kmer)]
            if head is None:
                return []
            mem.clear()
            node: Optional[SeedLoc] = head
            while node is not None:
                mem.append((node.seq_ind, node.res_ind))
                node = node.next
            return mem

    lookup = Lookup

    def hash_kmer(self, kmer: bytes) -> int:
        h = 0
        last_pow = len(kmer) - 1
        for i, b in enumerate(kmer):
            h += SEED_ALPHA_NUMS[b - ord("A")] * self.powers[last_pow - i]
        return h

    hashKmer = hash_kmer

    def unhash_kmer(self, hashed: int) -> bytes:
        residues = bytearray()
        base = SEED_ALPHA_SIZE
        h = hashed
        for _ in range(self.seed_size):
            ones_zeroed = (h // base) * base
            digit = h - ones_zeroed
            residues.append(REVERSE_SEED_ALPHA_NUMS[digit])
            h //= base
        residues.reverse()
        return bytes(residues)

    unhashKmer = unhash_kmer
