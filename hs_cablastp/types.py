"""Core data structures for HS-CaBLASTP."""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple


@dataclass
class EditOp:
    op_type: str   # 'INSERT', 'DELETE', or 'SUBSTITUTE'
    position: int  # Position in the parent's matched segment (0-based)
    char: str      # Amino acid character (residue for INSERT/SUBSTITUTE; '' for DELETE)


@dataclass
class TreeNode:
    node_id: int
    is_root: bool
    sequence: Optional[str]               # Set only when is_root is True
    parent_id: Optional[int]
    children: List[int] = field(default_factory=list)
    diff_script: List[EditOp] = field(default_factory=list)
    depth: int = 0
    ref_original_seq: str = ""             # "<fasta_id>:<start>-<end>" for traceability
    parent_start: int = 0                  # Start offset within the parent's sequence (children only)
    parent_end: int = 0                    # End offset within the parent's sequence (children only)


# Murphy-10 reduced amino-acid alphabet (Murphy et al. 2000, "Simplified amino
# acid alphabets for protein fold recognition and implications for folding").
# Groups residues by physicochemical similarity into 10 classes; lookups hash
# similar-but-non-identical k-mers to the same bucket, so e.g. an `L` -> `I`
# substitution still triggers a seed hit. The downstream ungapped_extend and
# identity check operate on the *original* residues, so reduced-alphabet
# seeding adds sensitivity without admitting false matches.
_REDUCE_MAP = {
    "L": "1", "V": "1", "I": "1", "M": "1",   # aliphatic hydrophobic
    "C": "2",
    "A": "3",
    "G": "4",
    "S": "5", "T": "5",                       # polar / small
    "P": "6",
    "F": "7", "Y": "7", "W": "7",             # aromatic
    "E": "8", "D": "8", "N": "8", "Q": "8",   # acidic + amide
    "B": "8", "Z": "8",                       # ambiguous: fold in with their groups
    "K": "9", "R": "9",                       # basic
    "H": "0",                                 # basic / aromatic
}


def _reduce_kmer(kmer: str) -> Optional[str]:
    """Return the Murphy-10 reduced form of `kmer`, or None if it contains an
    unmappable residue (gap, stop, or ambiguity char not in _REDUCE_MAP)."""
    out = []
    for c in kmer:
        m = _REDUCE_MAP.get(c)
        if m is None:
            return None
        out.append(m)
    return "".join(out)


class SeedTable:
    """Maps Murphy-10 reduced k-mer string -> list of (node_id, offset).

    Indexes k-mers from both root and child nodes. The offset is relative to
    the start of that node's own (reconstructed) sequence, so a hit on a child
    points directly at the child rather than forcing a descent from the root.

    Keys are Murphy-10 reduced (10 chars vs the full 20-letter alphabet), so
    chemically similar k-mers collide and all become candidates from a single
    lookup. The downstream extension/identity check operates on the original
    residues, so reduced-alphabet seeding adds sensitivity without admitting
    false matches.
    """

    def __init__(self, k: int):
        self.k = k
        self.table: Dict[str, List[Tuple[int, int]]] = defaultdict(list)

    def add_sequence(self, node_id: int, sequence: str) -> None:
        """Index every clean k-mer of `sequence` under `node_id`."""
        k = self.k
        for i in range(len(sequence) - k + 1):
            kmer = sequence[i:i + k]
            reduced = _reduce_kmer(kmer)
            if reduced is None:
                continue
            self.table[reduced].append((node_id, i))

    # Back-compat alias; roots and children index the same way.
    add_root_sequence = add_sequence
    add_node_sequence = add_sequence

    def lookup(self, kmer: str) -> List[Tuple[int, int]]:
        reduced = _reduce_kmer(kmer)
        if reduced is None:
            return []
        return self.table.get(reduced, [])


class CompressedDB:
    """The hierarchical forest plus the flat coarse FASTA of roots."""

    def __init__(self):
        self.forest: Dict[int, TreeNode] = {}
        self._next_id: int = 0

    def new_node(self, **kwargs) -> TreeNode:
        node = TreeNode(node_id=self._next_id, **kwargs)
        self.forest[self._next_id] = node
        self._next_id += 1
        return node

    def roots(self):
        for node in self.forest.values():
            if node.is_root:
                yield node
