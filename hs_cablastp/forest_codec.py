"""Compact binary serializer for the hs-cablastp forest.

Replaces pickle's per-object framing (50+ bytes of type tags, attribute
names, and refcount-friendly markers per TreeNode and EditOp) with a
hand-rolled layout. The wins, in order of magnitude:

1. Intern `ref_original_seq` strings into a single table. Every fragment of
   the same input shares the same orig id, so a 23k-input compression that
   produces 30k nodes only needs ~23k unique refs instead of 30k copies.
2. Pack each `EditOp` into 1 byte tag + varint position + 1 byte char (≤6
   bytes, usually 3) instead of pickle's ~30 bytes per dataclass.
3. Varint integers for ids, children, coords — short ints become 1–2 bytes
   instead of pickle's 5-byte BININT2/BININT framing.

Layout (no padding, all little-endian where it matters):

  magic        b"HSFOREST"
  version      uint8                          (=2)
  num_refs     varuint
  refs         [varuint len, bytes utf-8]*
  num_nodes    varuint
  nodes        node*

  node:
    node_id        varuint
    flags          uint8   (bit0 is_root, bit2 has_parent; bit1 unused in v2)
    parent_block   [varuint parent_id, varuint parent_start, varuint parent_end]
                                                  iff has_parent
    depth          uint8
    ref_idx        varuint
    num_children   varuint
    children       varuint*
    num_ops        varuint
    ops            op*

Root residues are NOT stored in this blob. They live in coarse.fasta (where
makeblastdb consumes them anyway) and `load_db` reattaches them to roots by
node id at load time. Save 50%+ of the v1 blob this way on Swiss-Prot-like
data where roots dominate the residue total.

  op:
    tag        uint8     (0=INSERT, 1=DELETE, 2=SUBSTITUTE)
    position   varuint
    char       uint8     (residue ord, or 0 for DELETE / empty)

Old pickle-format forest.pkl files won't load with this codec; rebuild the
DB. The magic-byte check makes the failure mode loud rather than silent.
"""

from __future__ import annotations

from typing import Dict, List, Tuple

from hs_cablastp.types import EditOp, TreeNode

_MAGIC = b"HSFOREST"
_VERSION = 2

_OP_TAG = {"INSERT": 0, "DELETE": 1, "SUBSTITUTE": 2}
_OP_NAME = ["INSERT", "DELETE", "SUBSTITUTE"]


def _write_varuint(out: bytearray, n: int) -> None:
    """Unsigned varint, LEB128-style. Negative inputs would crash — we never
    have negative integers in the forest (node ids, positions, depths are all
    >= 0)."""
    while n >= 128:
        out.append((n & 0x7F) | 0x80)
        n >>= 7
    out.append(n)


def _read_varuint(buf: bytes, pos: int) -> Tuple[int, int]:
    n = 0
    shift = 0
    while True:
        b = buf[pos]
        pos += 1
        n |= (b & 0x7F) << shift
        if not (b & 0x80):
            return n, pos
        shift += 7


def pack_forest(forest: Dict[int, TreeNode]) -> bytes:
    # Pass 1: intern `ref_original_seq` strings.
    refs: Dict[str, int] = {}
    for node in forest.values():
        r = node.ref_original_seq or ""
        if r not in refs:
            refs[r] = len(refs)

    out = bytearray()
    out.extend(_MAGIC)
    out.append(_VERSION)

    _write_varuint(out, len(refs))
    for r in refs:
        b = r.encode("utf-8")
        _write_varuint(out, len(b))
        out.extend(b)

    _write_varuint(out, len(forest))
    # Sort by node_id for deterministic output; saves nothing on size but
    # makes byte-identical reproducibility checks possible.
    for nid in sorted(forest):
        node = forest[nid]
        _write_varuint(out, nid)
        # bit 1 (has_sequence) is unused in v2: root residues live in
        # coarse.fasta only and load_db reattaches them by node id.
        flags = 0
        if node.is_root:
            flags |= 1
        if node.parent_id is not None:
            flags |= 4
        out.append(flags)

        if node.parent_id is not None:
            _write_varuint(out, node.parent_id)
            _write_varuint(out, node.parent_start)
            _write_varuint(out, node.parent_end)

        out.append(node.depth)
        _write_varuint(out, refs[node.ref_original_seq or ""])

        _write_varuint(out, len(node.children))
        for c in node.children:
            _write_varuint(out, c)

        _write_varuint(out, len(node.diff_script))
        for op in node.diff_script:
            out.append(_OP_TAG[op.op_type])
            _write_varuint(out, op.position)
            out.append(ord(op.char) if op.char else 0)

    return bytes(out)


def unpack_forest(buf: bytes) -> Dict[int, TreeNode]:
    if not buf.startswith(_MAGIC):
        raise ValueError(
            "forest blob is not in HSFOREST binary format; old pickle DBs "
            "must be rebuilt with the current codec"
        )
    pos = len(_MAGIC)
    ver = buf[pos]
    pos += 1
    if ver != _VERSION:
        raise ValueError(f"unsupported HSFOREST version {ver}")

    n_refs, pos = _read_varuint(buf, pos)
    refs: List[str] = []
    for _ in range(n_refs):
        sz, pos = _read_varuint(buf, pos)
        refs.append(buf[pos:pos + sz].decode("utf-8"))
        pos += sz

    n_nodes, pos = _read_varuint(buf, pos)
    forest: Dict[int, TreeNode] = {}
    for _ in range(n_nodes):
        nid, pos = _read_varuint(buf, pos)
        flags = buf[pos]
        pos += 1
        is_root = bool(flags & 1)
        has_parent = bool(flags & 4)

        # Root residues come from coarse.fasta in v2; caller (load_db)
        # reattaches them after this returns. Leave sequence=None here.
        sequence: str | None = None

        parent_id = None
        parent_start = 0
        parent_end = 0
        if has_parent:
            parent_id, pos = _read_varuint(buf, pos)
            parent_start, pos = _read_varuint(buf, pos)
            parent_end, pos = _read_varuint(buf, pos)

        depth = buf[pos]
        pos += 1
        ref_idx, pos = _read_varuint(buf, pos)

        n_children, pos = _read_varuint(buf, pos)
        children: List[int] = []
        for _ in range(n_children):
            c, pos = _read_varuint(buf, pos)
            children.append(c)

        n_ops, pos = _read_varuint(buf, pos)
        diff_script: List[EditOp] = []
        for _ in range(n_ops):
            tag = buf[pos]
            pos += 1
            position, pos = _read_varuint(buf, pos)
            ch_byte = buf[pos]
            pos += 1
            ch = chr(ch_byte) if ch_byte else ""
            diff_script.append(EditOp(
                op_type=_OP_NAME[tag], position=position, char=ch,
            ))

        forest[nid] = TreeNode(
            node_id=nid,
            is_root=is_root,
            sequence=sequence,
            parent_id=parent_id,
            children=children,
            diff_script=diff_script,
            depth=depth,
            ref_original_seq=refs[ref_idx],
            parent_start=parent_start,
            parent_end=parent_end,
        )
    return forest
