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


class SeedTable:
    """Maps k-mer string -> list of (node_id, offset_in_node).

    Indexes k-mers from both root and child nodes. The offset is relative to
    the start of that node's own (reconstructed) sequence, so a hit on a child
    points directly at the child rather than forcing a descent from the root.
    """

    def __init__(self, k: int):
        self.k = k
        self.table: Dict[str, List[Tuple[int, int]]] = defaultdict(list)

    def add_sequence(self, node_id: int, sequence: str) -> None:
        """Index every clean k-mer of `sequence` under `node_id`."""
        k = self.k
        for i in range(len(sequence) - k + 1):
            kmer = sequence[i:i + k]
            if "*" in kmer or "-" in kmer:
                continue
            self.table[kmer].append((node_id, i))

    # Back-compat alias; roots and children index the same way.
    add_root_sequence = add_sequence
    add_node_sequence = add_sequence

    def lookup(self, kmer: str) -> List[Tuple[int, int]]:
        return self.table.get(kmer, [])


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
