"""EditScript: substitution/insertion/deletion modifications (port of seqdiff.go)."""

from __future__ import annotations

from typing import List, Optional, Tuple

MOD_SUBSTITUTION = 0
MOD_DELETION = 1
MOD_INSERTION = 2


def _byte_to_mod_kind(ch: str) -> int:
    if ch == "s":
        return MOD_SUBSTITUTION
    if ch == "i":
        return MOD_INSERTION
    if ch == "d":
        return MOD_DELETION
    raise ValueError("Invalid mod kind: %r" % ch)


class _Mod:
    __slots__ = ("kind", "start", "end", "residues")

    def __init__(self, kind: int):
        self.kind = kind
        self.start = 0
        self.end = 0
        self.residues = bytearray()

    def add_residue(self, residue: int) -> None:
        if self.kind != MOD_DELETION:
            self.residues.append(residue)
        if self.kind != MOD_INSERTION:
            self.end += 1

    def __str__(self) -> str:
        body = self.residues.decode("ascii", errors="replace")
        if self.kind == MOD_SUBSTITUTION:
            return "s(%d,%d)%s" % (self.start, self.end, body)
        if self.kind == MOD_INSERTION:
            return "i(%d,%d)%s" % (self.start, self.end, body)
        if self.kind == MOD_DELETION:
            return "d(%d,%d)%s" % (self.start, self.end, body)
        raise ValueError("Invalid kind '%d' for an EditScript modification." % self.kind)


class EditScript:
    """A list of substitutions/insertions/deletions between two sequences."""

    def __init__(self, mods: Optional[List[_Mod]] = None):
        self.mods: List[_Mod] = mods if mods is not None else []

    def apply(self, from_seq: bytes) -> bytes:
        out = bytearray()
        last_end = 0
        for mod in self.mods:
            out.extend(from_seq[last_end:mod.start])
            out.extend(mod.residues)
            last_end = mod.end
        out.extend(from_seq[last_end:])
        return bytes(out)

    # `Apply` for parity with the Go API.
    def Apply(self, from_seq: bytes) -> bytes:  # noqa: N802
        return self.apply(from_seq)

    def __str__(self) -> str:
        parts: List[str] = []
        last_dist = 0
        for m in self.mods:
            dist = m.start - last_dist
            last_dist = m.start
            body = m.residues.decode("ascii", errors="replace").upper()
            if m.kind == MOD_SUBSTITUTION:
                parts.append("s%d%s" % (dist, body))
            elif m.kind == MOD_INSERTION:
                parts.append("i%d%s" % (dist, body))
            elif m.kind == MOD_DELETION:
                parts.append("d%d%s" % (dist, "-" * (m.end - m.start)))
            else:
                raise ValueError(
                    "Invalid kind '%d' for an EditScript modification." % m.kind
                )
        return "".join(parts)


def new_edit_script(alignment: Tuple[bytes, bytes]) -> EditScript:
    """Build an EditScript from a pair of aligned sequences."""
    return _new_edit_script(alignment[0], alignment[1])


def _new_edit_script(from_seq: bytes, to_seq: bytes) -> EditScript:
    if len(from_seq) != len(to_seq):
        raise ValueError(
            "A new edit script can only be generated with two sequences "
            "of equal length. Lengths of %d and %d were provided."
            % (len(from_seq), len(to_seq))
        )

    mods: List[_Mod] = []
    cur: Optional[_Mod] = None
    from_index = 0
    dash = ord("-")

    for i in range(len(from_seq)):
        a, b = from_seq[i], to_seq[i]
        if a == b:
            new_kind = -1
        elif a == dash:
            new_kind = MOD_INSERTION
        elif b == dash:
            new_kind = MOD_DELETION
        else:
            new_kind = MOD_SUBSTITUTION

        if cur is not None:
            if cur.kind == new_kind:
                cur.add_residue(b)
            else:
                mods.append(cur)
                if new_kind == -1:
                    cur = None
                else:
                    cur = _Mod(new_kind)
                    cur.start = from_index
                    cur.end = from_index
                    cur.add_residue(b)
        elif new_kind != -1:
            cur = _Mod(new_kind)
            cur.start = from_index
            cur.end = from_index
            cur.add_residue(b)

        if a != dash:
            from_index += 1

    if cur is not None:
        mods.append(cur)

    return EditScript(mods)


def parse_edit_script(script: str) -> EditScript:
    """Parse the textual edit-script form back into an EditScript."""
    mods: List[_Mod] = []
    cur: Optional[_Mod] = None
    i = 0
    n = len(script)
    while i < n:
        b = script[i]
        if b in ("s", "i", "d"):
            if cur is not None:
                mods.append(cur)
            new_mod = _Mod(_byte_to_mod_kind(b))
            digits: List[str] = []
            j = i + 1
            while j < n and script[j].isdigit():
                digits.append(script[j])
                j += 1
            if not digits:
                raise ValueError(
                    "Expected an offset number after '%s' in column %d of '%s'."
                    % (b, i, script)
                )
            try:
                num = int("".join(digits))
            except ValueError as err:
                raise ValueError(
                    "Expected an offset number after '%s' in column %d of '%s', "
                    "but got '%s' instead." % (b, i, script, "".join(digits))
                ) from err
            new_mod.start = num
            if cur is not None:
                new_mod.start += cur.start
            cur = new_mod
            cur.end = cur.start
            i = j
            continue

        if b.isdigit():
            raise ValueError(
                "Expected a residue at column %d in '%s', but got a number "
                "'%s' instead." % (i, script, b)
            )
        if cur is None:
            raise ValueError(
                "Expected 's', 'i' or 'd' but got '%s' at column %d in '%s'."
                % (b, i, script)
            )
        cur.add_residue(ord(b))
        i += 1

    if cur is not None:
        mods.append(cur)
    return EditScript(mods)
