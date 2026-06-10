"""FASTA reader/writer and original-sequence iteration (port of fasta.go).

The Go code relied on `github.com/TuftsBCB/io/fasta` for parsing; this module
provides a minimal equivalent so cablastp can stand on its own.
"""

from __future__ import annotations

import gzip
import io
import threading
from dataclasses import dataclass
from queue import Queue
from typing import IO, Iterable, Iterator, List, Optional, Sequence, Union


@dataclass
class FastaSequence:
    """A simple FASTA record."""
    name: str
    residues: bytes

    def bytes(self) -> bytes:
        return self.residues


class FastaReader:
    """Iterate over (name, residues) records in a FASTA stream."""

    def __init__(self, stream: IO):
        self._stream = stream
        # We need text reading; if the stream is binary, wrap it.
        if isinstance(stream, (io.RawIOBase, io.BufferedIOBase, gzip.GzipFile)):
            self._text = io.TextIOWrapper(stream, encoding="utf-8", newline="")
        else:
            self._text = stream
        self._pending: Optional[str] = None

    def read(self) -> Optional[FastaSequence]:
        """Return the next FASTA record, or None at EOF."""
        # Skip blank lines until we find a header.
        line = self._pending
        self._pending = None
        if line is None:
            line = self._text.readline()
            while line and not line.strip():
                line = self._text.readline()
        if not line:
            return None
        if not line.startswith(">"):
            raise ValueError("FASTA record must begin with '>': %r" % line)
        name = line[1:].strip()

        residues = bytearray()
        while True:
            nxt = self._text.readline()
            if not nxt:
                break
            if nxt.startswith(">"):
                self._pending = nxt
                break
            residues.extend(nxt.strip().encode("ascii"))

        return FastaSequence(name=name, residues=bytes(residues))

    def __iter__(self) -> Iterator[FastaSequence]:
        while True:
            seq = self.read()
            if seq is None:
                return
            yield seq

    def detach(self) -> None:
        """Release an internally-wrapped binary stream without closing it.

        When constructed over a binary handle we wrap it in a TextIOWrapper; if
        that wrapper is later garbage-collected it closes the underlying handle.
        A caller that must keep the handle open (e.g. CoarseDB reads the coarse
        FASTA and then appends to the same handle) calls this after reading to
        EOF to take back ownership of the stream. No-op when no wrapper was
        created (the stream was already text). The reader must not be used after.
        """
        if isinstance(self._text, io.TextIOWrapper) and self._text is not self._stream:
            self._text.detach()
            self._text = None


class FastaWriter:
    """Write FASTA records.

    `asterisk=True` emits a trailing '*' after each sequence body, matching
    the behaviour of TuftsBCB/io/fasta's Asterisk option.
    """

    def __init__(self, stream: IO, asterisk: bool = False, line_width: int = 60):
        self._stream = stream
        self.asterisk = asterisk
        self.line_width = line_width

    def _write(self, s: str) -> None:
        # Fall back to a write-as-bytes attempt only if the stream rejects strings.
        try:
            self._stream.write(s)
        except TypeError:
            self._stream.write(s.encode("utf-8"))

    def write(self, seq: FastaSequence) -> None:
        self._write(">" + seq.name + "\n")
        residues = seq.residues if isinstance(seq.residues, (bytes, bytearray)) else bytes(seq.residues)
        text = residues.decode("ascii")
        if self.line_width and self.line_width > 0:
            for i in range(0, len(text), self.line_width):
                self._write(text[i:i + self.line_width] + "\n")
        else:
            self._write(text + "\n")
        if self.asterisk:
            self._write("*\n")

    def flush(self) -> None:
        if hasattr(self._stream, "flush"):
            self._stream.flush()


@dataclass
class ReadOriginalSeq:
    """Value yielded while iterating a FASTA file as OriginalSeq."""
    seq: Optional["OriginalSeq"] = None
    err: Optional[BaseException] = None


def read_original_seqs(
    file_name: str,
    ignore: Sequence[int] = (),
) -> "Queue[ReadOriginalSeq]":
    """Read a FASTA file and stream OriginalSeq records onto a queue.

    `ignore` is a sequence of byte values (e.g., b'JOU') that should be replaced
    with 'X' in the residue strings.

    Returns a Queue. A `ReadOriginalSeq` with both `seq` and `err` equal to
    None signals end-of-stream.
    """
    # Local import to avoid circularity at module import time.
    from cablastp.seq import OriginalSeq, new_fasta_seq

    if file_name.endswith(".gz"):
        raw = gzip.open(file_name, "rb")
    else:
        raw = open(file_name, "rb")

    reader = FastaReader(raw)
    queue: "Queue[ReadOriginalSeq]" = Queue(maxsize=200)

    ignore_set = {int(b) for b in ignore}

    def producer() -> None:
        try:
            i = 0
            while True:
                try:
                    sequence = reader.read()
                except Exception as err:
                    queue.put(ReadOriginalSeq(seq=None, err=err))
                    queue.put(ReadOriginalSeq(seq=None, err=None))
                    return
                if sequence is None:
                    queue.put(ReadOriginalSeq(seq=None, err=None))
                    return

                if ignore_set:
                    residues = bytearray(sequence.residues)
                    for idx, residue in enumerate(residues):
                        if residue in ignore_set:
                            residues[idx] = ord("X")
                    sequence = FastaSequence(sequence.name, bytes(residues))

                queue.put(
                    ReadOriginalSeq(seq=new_fasta_seq(i, sequence), err=None)
                )
                i += 1
        finally:
            try:
                raw.close()
            except Exception:
                pass

    threading.Thread(target=producer, daemon=True).start()
    return queue


def iter_original_seqs(
    file_name: str,
    ignore: Sequence[int] = (),
) -> Iterator["OriginalSeq"]:
    """Convenience iterator over OriginalSeq records from a FASTA file."""
    queue = read_original_seqs(file_name, ignore)
    while True:
        item = queue.get()
        if item.err is not None:
            raise item.err
        if item.seq is None:
            return
        yield item.seq
