"""CoarseDB: the deduplicated reference sequences (port of coarse.go + io.go)."""

from __future__ import annotations

import csv
import gzip
import io
import os
import struct
import threading
import time
from typing import IO, Dict, List, Optional, TYPE_CHECKING

from cablastp.fasta import FastaReader
from cablastp.links import LinkToCompressed
from cablastp.misc import vprintf, vprintln
from cablastp.seeds import Seeds
from cablastp.seq import CoarseSeq, OriginalSeq, new_fasta_coarse_seq

if TYPE_CHECKING:
    from cablastp.compressed import CompressedDB
    from cablastp.db import DB


# Hard-coded file names that make up a coarse database on disk.
FILE_COARSE_FASTA = "coarse.fasta"
FILE_COARSE_FASTA_INDEX = "coarse.fasta.index"
FILE_COARSE_LINKS = "coarse.links"
FILE_COARSE_PLAIN_LINKS = "coarse.links.plain"
FILE_COARSE_LINKS_INDEX = "coarse.links.index"
FILE_COARSE_SEEDS = "coarse.seeds"
FILE_COARSE_PLAIN_SEEDS = "coarse.seeds.plain"


class CoarseDB:
    """Set of unique reference sequences plus the seeds index."""

    def __init__(self):
        self.seqs: List[CoarseSeq] = []
        self.seeds: Optional[Seeds] = None
        self.fasta_cache: Dict[int, CoarseSeq] = {}
        self._fasta_index_size: int = 0

        self.file_fasta: Optional[IO] = None
        self.file_fasta_index: Optional[IO] = None
        self.file_seeds: Optional[IO] = None
        self.file_links: Optional[IO] = None
        self.file_links_index: Optional[IO] = None

        self._seq_lock: Optional[threading.RLock] = None
        self._read_only: bool = False
        self._seqs_read: int = 0

        self._plain: bool = False
        self._plain_links: Optional[IO] = None
        self._plain_seeds: Optional[IO] = None

    # ---------------------------------------------------------- factories

    @classmethod
    def open_for_write(cls, append: bool, db: "DB") -> "CoarseDB":
        vprintln("\tOpening coarse database...")

        cdb = cls()
        cdb.seeds = Seeds(db.MapSeedSize, db.SeedLowComplexity)
        cdb._seq_lock = threading.RLock()
        cdb._read_only = db.ReadOnly
        cdb._plain = db.SavePlain

        cdb.file_fasta = _open_for_write(append, db.file_path(FILE_COARSE_FASTA))
        cdb.file_fasta_index = _open_for_write(append, db.file_path(FILE_COARSE_FASTA_INDEX))
        cdb.file_seeds = _open_for_write(append, db.file_path(FILE_COARSE_SEEDS))
        cdb.file_links = _open_for_write(append, db.file_path(FILE_COARSE_LINKS))
        cdb.file_links_index = _open_for_write(append, db.file_path(FILE_COARSE_LINKS_INDEX))

        cdb._fasta_index_size = os.fstat(cdb.file_fasta_index.fileno()).st_size

        if cdb._plain:
            cdb._plain_links = _open_for_write(append, db.file_path(FILE_COARSE_PLAIN_LINKS))
            cdb._plain_seeds = _open_for_write(append, db.file_path(FILE_COARSE_PLAIN_SEEDS))

        if append:
            cdb._load()
            for f in (cdb.file_seeds, cdb.file_links, cdb.file_links_index):
                f.truncate(0)
                f.seek(0, os.SEEK_SET)
            if cdb._plain:
                cdb._plain_seeds.truncate(0)
                cdb._plain_seeds.seek(0, os.SEEK_SET)
                cdb._plain_links.truncate(0)
                cdb._plain_links.seek(0, os.SEEK_SET)

        vprintln("\tDone opening coarse database.")
        return cdb

    @classmethod
    def open_for_read(cls, db: "DB") -> "CoarseDB":
        vprintln("\tOpening coarse database...")
        cdb = cls()
        cdb.seeds = Seeds(db.MapSeedSize, db.SeedLowComplexity)
        cdb._read_only = False
        cdb._plain = db.SavePlain

        cdb.file_fasta = open(db.file_path(FILE_COARSE_FASTA), "rb")
        cdb.file_fasta_index = open(db.file_path(FILE_COARSE_FASTA_INDEX), "rb")
        cdb.file_links = open(db.file_path(FILE_COARSE_LINKS), "rb")
        cdb.file_links_index = open(db.file_path(FILE_COARSE_LINKS_INDEX), "rb")

        cdb._fasta_index_size = os.fstat(cdb.file_fasta_index.fileno()).st_size

        vprintln("\tDone opening coarse database.")
        return cdb

    # ---------------------------------------------------------- mutation

    def add(self, oseq: bytes) -> "tuple[int, CoarseSeq]":
        """Add `oseq` as a new coarse sequence and seed it. Returns (id, seq)."""
        with self._seq_lock:
            id = len(self.seqs)
            cor_seq = CoarseSeq(id, "", oseq)
            self.seqs.append(cor_seq)
        self.seeds.add(id, cor_seq)
        return id, cor_seq

    Add = add

    def coarse_seq_get(self, i: int) -> CoarseSeq:
        with self._seq_lock:
            return self.seqs[i]

    CoarseSeqGet = coarse_seq_get

    # ----------------------------------------------------------- reading

    def num_sequences(self) -> int:
        return self._fasta_index_size // 8

    NumSequences = num_sequences

    def read_coarse_seq(self, id: int) -> CoarseSeq:
        if id in self.fasta_cache:
            return self.fasta_cache[id]

        off = self._coarse_offset(id)
        self.file_fasta.seek(off, os.SEEK_SET)

        header = self.file_fasta.readline()
        if not header.startswith(b">"):
            raise IOError(
                "Could not scan coarse sequence %d: missing '>' header." % id
            )
        try:
            cor_seq_id = int(header[1:].strip())
        except ValueError as err:
            raise IOError(
                "Could not scan coarse sequence %d: %s" % (id, err)
            ) from err
        if cor_seq_id != id:
            raise IOError(
                "Expected to read coarse sequence %d but read coarse "
                "sequence %d instead." % (id, cor_seq_id)
            )
        residues = self.file_fasta.readline().strip()

        coarse_seq = CoarseSeq(id, "", residues)
        self.fasta_cache[id] = coarse_seq
        return coarse_seq

    ReadCoarseSeq = read_coarse_seq

    def expand(
        self,
        comdb: "CompressedDB",
        id: int,
        start: int,
        end: int,
    ) -> List[OriginalSeq]:
        off = self._link_offset(id)
        self.file_links.seek(off, os.SEEK_SET)
        num_links = struct.unpack(">I", _read_exact(self.file_links, 4))[0]

        s = start & 0xFFFF
        e = end & 0xFFFF

        seen = set()
        oseqs: List[OriginalSeq] = []
        for _ in range(num_links):
            link = self._read_link()
            if e < link.coarse_start or s > link.coarse_end:
                continue
            if link.org_seq_id in seen:
                continue
            seen.add(link.org_seq_id)
            oseq = comdb.read_seq(self, int(link.org_seq_id))
            oseqs.append(oseq)
        return oseqs

    Expand = expand

    # ---------------------------------------------------------- offsets

    def _coarse_offset(self, id: int) -> int:
        try_off = id * 8
        self.file_fasta_index.seek(try_off, os.SEEK_SET)
        return struct.unpack(">q", _read_exact(self.file_fasta_index, 8))[0]

    def _link_offset(self, id: int) -> int:
        try_off = id * 8
        self.file_links_index.seek(try_off, os.SEEK_SET)
        return struct.unpack(">q", _read_exact(self.file_links_index, 8))[0]

    def _read_link(self) -> LinkToCompressed:
        org_seq_id = struct.unpack(">I", _read_exact(self.file_links, 4))[0]
        coarse_start = struct.unpack(">H", _read_exact(self.file_links, 2))[0]
        coarse_end = struct.unpack(">H", _read_exact(self.file_links, 2))[0]
        return LinkToCompressed(org_seq_id, coarse_start, coarse_end)

    # ---------------------------------------------------------- close/save

    def read_close(self) -> None:
        for f in (self.file_fasta, self.file_fasta_index,
                  self.file_links, self.file_links_index):
            if f:
                f.close()

    def write_close(self) -> None:
        for f in (self.file_fasta, self.file_fasta_index,
                  self.file_seeds, self.file_links, self.file_links_index):
            if f:
                f.close()
        if self._plain:
            if self._plain_links:
                self._plain_links.close()
            if self._plain_seeds:
                self._plain_seeds.close()

    def save(self) -> None:
        with self._seq_lock:
            errors: List[BaseException] = []
            threads: List[threading.Thread] = []

            def runner(fn):
                try:
                    fn()
                except BaseException as err:
                    errors.append(err)

            jobs = [self._save_fasta, self._save_links]
            if not self._read_only:
                jobs.append(self._save_seeds)
            if self._plain:
                jobs.append(self._save_links_plain)
                if not self._read_only:
                    jobs.append(self._save_seeds_plain)

            for job in jobs:
                t = threading.Thread(target=runner, args=(job,))
                threads.append(t)
                t.start()
            for t in threads:
                t.join()

            if errors:
                raise errors[0]

    # ------------------------------------------------------- file format

    def _load(self) -> None:
        self._read_fasta()
        self._read_links()
        if self.file_seeds is not None:
            self._read_seeds()

    def _read_fasta(self) -> None:
        vprintf("\t\tReading %s...\n", FILE_COARSE_FASTA)
        start = time.time()

        self.file_fasta.seek(0, os.SEEK_SET)
        reader = FastaReader(self.file_fasta)
        i = 0
        while True:
            seq = reader.read()
            if seq is None:
                break
            self.seqs.append(new_fasta_coarse_seq(i, seq))
            i += 1
        self._seqs_read = len(self.seqs)
        # FastaReader wrapped our binary handle in a TextIOWrapper; releasing it
        # (we've consumed to EOF, so the handle is positioned at end) keeps
        # self.file_fasta open for the subsequent --append _save_fasta, which
        # fstat()s and writes to this same handle. Without this the wrapper would
        # close the handle on GC and _save_fasta would raise on a closed file.
        reader.detach()

        vprintf("\t\tDone reading %s (%.3fs).\n", FILE_COARSE_FASTA, time.time() - start)

    def _save_fasta(self) -> None:
        vprintf("Writing %s...\n", FILE_COARSE_FASTA)
        vprintf("Writing %s...\n", FILE_COARSE_FASTA_INDEX)
        start = time.time()

        byte_off = 0
        if self._fasta_index_size > 0:
            byte_off = os.fstat(self.file_fasta.fileno()).st_size

        for i in range(self._seqs_read, len(self.seqs)):
            body = b"> %d\n%s\n" % (i, bytes(self.seqs[i].residues))
            self.file_fasta.write(body)
            self.file_fasta_index.write(struct.pack(">q", byte_off))
            byte_off += len(body)

        vprintf("Done writing %s (%.3fs).\n", FILE_COARSE_FASTA, time.time() - start)
        vprintf("Done writing %s (%.3fs).\n", FILE_COARSE_FASTA_INDEX, time.time() - start)

    def _read_seeds(self) -> None:
        from cablastp.seeds import SeedLoc

        vprintf("\t\tReading %s... (this could take a while)\n", FILE_COARSE_SEEDS)
        start = time.time()

        self.file_seeds.seek(0, os.SEEK_SET)
        with gzip.GzipFile(fileobj=self.file_seeds, mode="rb") as gr:
            while True:
                hdr = gr.read(4)
                if len(hdr) < 4:
                    break
                hash_val = struct.unpack(">I", hdr)[0]
                cnt = struct.unpack(">I", _read_exact(gr, 4))[0]
                for _ in range(cnt):
                    seq_ind = struct.unpack(">I", _read_exact(gr, 4))[0]
                    res_ind = struct.unpack(">H", _read_exact(gr, 2))[0]
                    sl = SeedLoc(seq_ind, res_ind)
                    if self.seeds.locs[hash_val] is None:
                        self.seeds.locs[hash_val] = sl
                    else:
                        node = self.seeds.locs[hash_val]
                        while node.next is not None:
                            node = node.next
                        node.next = sl

        vprintf("\t\tDone reading %s (%.3fs).\n", FILE_COARSE_SEEDS, time.time() - start)

    def _save_seeds(self) -> None:
        vprintf("Writing %s... (this could take a while)\n", FILE_COARSE_SEEDS)
        start = time.time()

        with gzip.GzipFile(fileobj=self.file_seeds, mode="wb", compresslevel=1) as gw:
            for i in range(self.seeds.powers[self.seeds.seed_size]):
                head = self.seeds.locs[i]
                if head is None:
                    continue
                gw.write(struct.pack(">i", i))
                cnt = 0
                node = head
                while node is not None:
                    cnt += 1
                    node = node.next
                gw.write(struct.pack(">i", cnt))
                node = head
                while node is not None:
                    gw.write(struct.pack(">I", node.seq_ind))
                    gw.write(struct.pack(">H", node.res_ind))
                    node = node.next

        vprintf("Done writing %s (%.3fs).\n", FILE_COARSE_SEEDS, time.time() - start)

    def _save_seeds_plain(self) -> None:
        vprintf("Writing %s...\n", FILE_COARSE_PLAIN_SEEDS)
        start = time.time()

        text_buf = io.StringIO()
        writer = csv.writer(text_buf, delimiter=",", lineterminator="\n")
        for i in range(self.seeds.powers[self.seeds.seed_size]):
            head = self.seeds.locs[i]
            if head is None:
                continue
            record: List[str] = [self.seeds.unhash_kmer(i).decode("ascii")]
            node = head
            while node is not None:
                record.append(str(node.seq_ind))
                record.append(str(node.res_ind))
                node = node.next
            writer.writerow(record)
        self._plain_seeds.write(text_buf.getvalue().encode("utf-8"))

        vprintf("Done writing %s (%.3fs).\n", FILE_COARSE_PLAIN_SEEDS, time.time() - start)

    def _read_links(self) -> None:
        vprintf("\t\tReading %s...\n", FILE_COARSE_LINKS)
        start = time.time()

        self.file_links.seek(0, os.SEEK_SET)
        coarse_seq_id = 0
        while True:
            buf = self.file_links.read(4)
            if len(buf) < 4:
                break
            cnt = struct.unpack(">i", buf)[0]
            for _ in range(cnt):
                new_link = self._read_link()
                self.seqs[coarse_seq_id]._add_link_locked(new_link)
            coarse_seq_id += 1

        vprintf("\t\tDone reading %s (%.3fs).\n", FILE_COARSE_LINKS, time.time() - start)

    def _save_links(self) -> None:
        vprintf("Writing %s...\n", FILE_COARSE_LINKS)
        vprintf("Writing %s...\n", FILE_COARSE_LINKS_INDEX)
        start = time.time()

        byte_off = 0
        for seq in self.seqs:
            buf = bytearray()

            cnt = 0
            link = seq.links
            while link is not None:
                cnt += 1
                link = link.next
            buf.extend(struct.pack(">i", cnt))

            link = seq.links
            while link is not None:
                buf.extend(struct.pack(">I", link.org_seq_id))
                buf.extend(struct.pack(">H", link.coarse_start))
                buf.extend(struct.pack(">H", link.coarse_end))
                link = link.next

            self.file_links.write(buf)
            self.file_links_index.write(struct.pack(">q", byte_off))
            byte_off += len(buf)

        vprintf("Done writing %s (%.3fs).\n", FILE_COARSE_LINKS, time.time() - start)
        vprintf("Done writing %s (%.3fs).\n", FILE_COARSE_LINKS_INDEX, time.time() - start)

    def _save_links_plain(self) -> None:
        vprintf("Writing %s...\n", FILE_COARSE_PLAIN_LINKS)
        start = time.time()

        text_buf = io.StringIO()
        writer = csv.writer(text_buf, delimiter=",", lineterminator="\n")
        for seq in self.seqs:
            record: List[str] = []
            link = seq.links
            while link is not None:
                record.append(str(link.org_seq_id))
                record.append(str(link.coarse_start))
                record.append(str(link.coarse_end))
                link = link.next
            writer.writerow(record)
        self._plain_links.write(text_buf.getvalue().encode("utf-8"))

        vprintf("Done writing %s (%.3fs).\n", FILE_COARSE_PLAIN_LINKS, time.time() - start)


def _open_for_write(append: bool, path: str) -> IO:
    if append:
        f = open(path, "r+b")
        f.seek(0, os.SEEK_END)
        return f
    return open(path, "w+b")


def _read_exact(stream, n: int) -> bytes:
    data = stream.read(n)
    if len(data) != n:
        raise IOError("Unexpected EOF (wanted %d bytes, got %d)" % (n, len(data)))
    return data
