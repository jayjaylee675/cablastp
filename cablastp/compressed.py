"""CompressedDB: per-original-sequence compressed records (port of compressed.go)."""

from __future__ import annotations

import csv
import io
import os
import queue
import struct
import threading
from typing import IO, List, Optional, TYPE_CHECKING

from cablastp.links import LinkToCoarse, new_link_to_coarse_no_diff
from cablastp.misc import vprintln, vprintf

if TYPE_CHECKING:
    from cablastp.coarse import CoarseDB
    from cablastp.db import DB
    from cablastp.seq import OriginalSeq

FILE_COMPRESSED = "compressed"
FILE_INDEX = "compressed.index"


class CompressedSeq:
    """An original sequence represented as an ordered list of LinkToCoarse."""

    __slots__ = ("id", "name", "links")

    def __init__(self, id: int, name: str, links: Optional[List[LinkToCoarse]] = None):
        self.id = id
        self.name = name
        self.links: List[LinkToCoarse] = links if links is not None else []

    def add(self, link: LinkToCoarse) -> None:
        self.links.append(link)

    Add = add

    def __str__(self) -> str:
        lines = []
        for link in self.links:
            lines.append(
                "coarse id: %d, start: %d, end: %d\n%s"
                % (link.coarse_seq_id, link.coarse_start, link.coarse_end, link.diff)
            )
        return "\n".join(lines)

    def decompress(self, coarse: "CoarseDB") -> "OriginalSeq":
        """Recover the original sequence by following each link."""
        from cablastp.seq import OriginalSeq
        from cablastp.seqdiff import parse_edit_script

        residues = bytearray()
        for lk in self.links:
            if lk.coarse_seq_id < 0 or lk.coarse_seq_id >= coarse.num_sequences():
                raise ValueError(
                    "Cannot decompress compressed sequence (id: %d), "
                    "because a link refers to an invalid coarse sequence id: %d."
                    % (self.id, lk.coarse_seq_id)
                )
            edit = parse_edit_script(lk.diff)
            coarse_seq = coarse.read_coarse_seq(int(lk.coarse_seq_id))
            sub = coarse_seq.residues[lk.coarse_start:lk.coarse_end]
            residues.extend(edit.apply(sub))
        return OriginalSeq(self.id, self.name, bytes(residues))

    Decompress = decompress


def new_compressed_seq(id: int, name: str) -> CompressedSeq:
    return CompressedSeq(id=id, name=name, links=[])


# Sentinel for the writer queue.
_WRITER_DONE = object()


class CompressedDB:
    """The compressed-side of a cablastp database (writes/reads on disk)."""

    def __init__(self):
        self.file: Optional[IO] = None
        self.index: Optional[IO] = None
        self._index_size: int = 0
        self._writer_queue: Optional["queue.Queue"] = None
        self._writer_thread: Optional[threading.Thread] = None
        self._writer_exc: Optional[BaseException] = None
        self._seq_cache: Optional[dict] = None
        self._closed = False

    # ------------------------------------------------------------------ open

    @classmethod
    def open_for_write(cls, append: bool, db: "DB") -> "CompressedDB":
        vprintln("\tOpening compressed database...")

        cdb = cls()
        cdb._writer_queue = queue.Queue(maxsize=500)

        compressed_path = db.file_path(FILE_COMPRESSED)
        index_path = db.file_path(FILE_INDEX)

        if append:
            cdb.file = open(compressed_path, "r+b")
            cdb.file.seek(0, os.SEEK_END)
            cdb.index = open(index_path, "r+b")
            cdb.index.seek(0, os.SEEK_END)
        else:
            cdb.file = open(compressed_path, "w+b")
            cdb.index = open(index_path, "w+b")

        cdb._index_size = os.fstat(cdb.index.fileno()).st_size

        cdb._writer_thread = threading.Thread(target=cdb._writer, daemon=True)
        cdb._writer_thread.start()

        vprintln("\tDone opening compressed database.")
        return cdb

    @classmethod
    def open_for_read(cls, db: "DB") -> "CompressedDB":
        vprintln("\tOpening compressed database...")
        cdb = cls()
        cdb._seq_cache = {}
        cdb.file = open(db.file_path(FILE_COMPRESSED), "rb")
        cdb.index = open(db.file_path(FILE_INDEX), "rb")
        cdb._index_size = os.fstat(cdb.index.fileno()).st_size
        vprintln("\tDone opening compressed database.")
        return cdb

    # ------------------------------------------------------------------ ops

    def num_sequences(self) -> int:
        return self._index_size // 8

    NumSequences = num_sequences

    def seq_get(self, coarsedb: "CoarseDB", org_seq_id: int) -> "OriginalSeq":
        if self._writer_queue is not None:
            raise RuntimeError(
                "A compressed database cannot be read while it is also "
                "being modified."
            )
        # Never `self._seq_cache or {}`: an empty instance dict is falsy, so that
        # rebinds `cache` to a throwaway local and the populated entries are lost
        # (the `is None` write-back never fires once _seq_cache is a dict). Seed
        # the real dict on first use, then operate on it directly so the cache
        # actually persists across calls.
        if self._seq_cache is None:
            self._seq_cache = {}
        cache = self._seq_cache
        if org_seq_id not in cache:
            cache[org_seq_id] = self.read_seq(coarsedb, org_seq_id)
        return cache[org_seq_id]

    SeqGet = seq_get

    def write(self, cseq: CompressedSeq) -> None:
        assert self._writer_queue is not None
        self._writer_queue.put(cseq)

    Write = write

    def read_close(self) -> None:
        if self.file:
            self.file.close()
        if self.index:
            self.index.close()

    def write_close(self) -> None:
        if self._closed:
            return
        self._closed = True
        if self._writer_queue is not None:
            self._writer_queue.put(_WRITER_DONE)
        if self._writer_thread is not None:
            self._writer_thread.join()
        # Surface any failure the writer thread hit. Without this, an exception
        # in the daemon writer would just kill the thread: join() returns
        # normally and the caller reports a truncated DB as success.
        if self._writer_exc is not None:
            exc = self._writer_exc
            self._writer_exc = None
            raise exc

    # ------------------------------------------------------------------ I/O

    def _org_seq_offset(self, id: int) -> int:
        try_off = id * 8
        self.index.seek(try_off, os.SEEK_SET)
        data = self.index.read(8)
        if len(data) != 8:
            raise IOError("Unexpected EOF reading compressed index for id %d" % id)
        return struct.unpack(">q", data)[0]

    def read_seq(self, coarsedb: "CoarseDB", org_seq_id: int) -> "OriginalSeq":
        off = self._org_seq_offset(org_seq_id)
        self.file.seek(off, os.SEEK_SET)
        return self.read_next_seq(coarsedb, org_seq_id)

    ReadSeq = read_seq

    def read_next_seq(self, coarsedb: "CoarseDB", org_seq_id: int) -> "OriginalSeq":
        # The original code uses a streaming CSV reader. We read one logical
        # line (which may contain newlines inside a CSV-quoted diff field).
        record = _read_csv_record(self.file)
        if record is None:
            raise IOError("[csv reader]: id out of range")
        cseq = _parse_compressed_record(org_seq_id, record)
        return cseq.decompress(coarsedb)

    ReadNextSeq = read_next_seq

    # ----------------------------------------------------------------- writer

    def _writer(self) -> None:
        # Run the write loop in this worker thread. Any exception here would
        # otherwise silently kill the daemon thread — join() returns normally and
        # the caller reports a truncated DB as success — so capture it for
        # write_close() to re-raise on the main thread, and always close the
        # handles so a partial DB is at least flushed and consistent on disk.
        try:
            self._writer_loop()
        except BaseException as exc:  # noqa: BLE001 — re-raised in write_close()
            self._writer_exc = exc
        finally:
            if self.index is not None:
                self.index.close()
            if self.file is not None:
                self.file.close()

    def _writer_loop(self) -> None:
        # Maintains the original-sequence ordering of incoming compressed
        # sequences before flushing to disk, matching the Go writer.
        saved: List[CompressedSeq] = []
        next_index = self.num_sequences()

        if self._index_size > 0:
            byte_offset = os.fstat(self.file.fileno()).st_size
        else:
            byte_offset = 0

        text_buf = io.StringIO()
        writer = csv.writer(text_buf, delimiter=",", lineterminator="\n",
                            quoting=csv.QUOTE_MINIMAL)

        while True:
            possible = self._writer_queue.get()
            if possible is _WRITER_DONE:
                break
            if not isinstance(possible, CompressedSeq):
                raise TypeError(
                    "BUG: writer received non-CompressedSeq value: %r" % possible
                )
            if possible.id < next_index:
                raise RuntimeError(
                    "BUG: Next sequence expected is '%d', but we have an "
                    "earlier sequence: %d" % (next_index, possible.id)
                )
            saved.append(possible)

            cseq, saved = _next_seq_to_write(next_index, saved)
            while cseq is not None:
                text_buf.seek(0)
                text_buf.truncate(0)

                record: List[str] = [cseq.name]
                for link in cseq.links:
                    record.extend([
                        str(link.coarse_seq_id),
                        str(link.coarse_start),
                        str(link.coarse_end),
                        link.diff,
                    ])

                writer.writerow(record)
                encoded = text_buf.getvalue().encode("utf-8")
                self.file.write(encoded)
                self.index.write(struct.pack(">q", byte_offset))

                byte_offset += len(encoded)
                next_index += 1
                cseq, saved = _next_seq_to_write(next_index, saved)


def _next_seq_to_write(next_index: int, saved: List[CompressedSeq]):
    for i, cseq in enumerate(saved):
        if cseq.id < next_index:
            raise RuntimeError(
                "Cannot keep sequences (%d) earlier than the next index (%d)"
                % (cseq.id, next_index)
            )
        if cseq.id == next_index:
            chosen = cseq
            return chosen, saved[:i] + saved[i + 1:]
    return None, saved


def _read_csv_record(stream: IO) -> Optional[List[str]]:
    """Read a single CSV record from a binary stream.

    Handles quoted fields containing newlines.
    """
    chunks = bytearray()
    in_quotes = False
    saw_data = False

    while True:
        b = stream.read(1)
        if not b:
            break
        saw_data = True
        c = b[0]
        chunks.append(c)
        if c == ord('"'):
            in_quotes = not in_quotes
        elif c == ord("\n") and not in_quotes:
            break

    if not saw_data:
        return None

    text = chunks.decode("utf-8")
    reader = csv.reader([text], delimiter=",")
    for row in reader:
        return row
    return None


def _parse_compressed_record(id: int, record: List[str]) -> CompressedSeq:
    cseq = CompressedSeq(id=id, name=record[0], links=[])
    i = 1
    while i + 3 < len(record):
        try:
            coarse_seq_id = int(record[i])
            coarse_start = int(record[i + 1])
            coarse_end = int(record[i + 2])
        except ValueError:
            return cseq
        link = new_link_to_coarse_no_diff(coarse_seq_id, coarse_start, coarse_end)
        link.diff = record[i + 3]
        cseq.add(link)
        i += 4
    return cseq
