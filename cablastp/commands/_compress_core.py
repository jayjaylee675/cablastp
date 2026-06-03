"""Core compression worker (port of cmd/cablastp-compress/compression.go)."""

from __future__ import annotations

import threading
from concurrent.futures import ThreadPoolExecutor
from queue import Queue
from typing import TYPE_CHECKING

from cablastp.compressed import CompressedSeq, new_compressed_seq
from cablastp.links import (
    LinkToCompressed,
    new_link_to_coarse,
    new_link_to_coarse_no_diff,
)
from cablastp.seq import OriginalSeq, seq_identity
from cablastp.commands._align import align_len, align_ungapped
from cablastp.commands._memory import Memory, new_memory
from cablastp.commands._nw import nw_align

if TYPE_CHECKING:
    from cablastp.coarse import CoarseDB
    from cablastp.db import DB


_DONE = object()


class CompressPool:
    """Pool of compression workers (one Memory arena per worker thread)."""

    def __init__(self, db: "DB", num_workers: int):
        self.db = db
        self._jobs: "Queue" = Queue(maxsize=200)
        self._closed = False
        self._workers = []
        self._memories = []
        for _ in range(max(1, num_workers)):
            mem = new_memory()
            self._memories.append(mem)
            t = threading.Thread(target=self._worker, args=(mem,), daemon=True)
            self._workers.append(t)
            t.start()

    def _worker(self, mem: Memory) -> None:
        while True:
            job = self._jobs.get()
            if job is _DONE:
                break
            org_seq_id, org_seq = job
            cseq = compress(self.db, org_seq_id, org_seq, mem)
            self.db.com_db.write(cseq)

    def compress(self, id: int, seq: OriginalSeq) -> int:
        self._jobs.put((id, seq))
        return id + 1

    def done(self) -> None:
        if self._closed:
            return
        self._closed = True
        for _ in self._workers:
            self._jobs.put(_DONE)
        for t in self._workers:
            t.join()


def compress(
    db: "DB",
    org_seq_id: int,
    org_seq: OriginalSeq,
    mem: Memory,
) -> CompressedSeq:
    """Compress a single original sequence.

    Walks the residues looking for K-mer hits in the seeds table, attempts
    extension on each, accepts the first match that meets the configured
    sequence-identity threshold, and adds a coarse-DB sequence for any
    leftover residues that don't match.
    """
    cseq = new_compressed_seq(org_seq_id, org_seq.name)

    coarsedb: "CoarseDB" = db.coarse_db
    map_seed_size = db.MapSeedSize
    ext_seed_size = db.ExtSeedSize
    olen = org_seq.Len()
    residues = org_seq.residues

    last_match = 0
    current = -1
    limit = olen - map_seed_size - ext_seed_size

    while True:
        current += 1
        if current >= limit:
            break

        kmer = residues[current:current + map_seed_size]
        seeds = coarsedb.seeds.lookup(kmer, mem.seeds)

        if db.LowComplexity > 0:
            skip = skip_low_complexity(
                residues[current:], db.MinMatchLen, db.LowComplexity,
            )
            if skip > 0:
                current += skip
                continue

        for seed_loc in seeds:
            cor_seq_id = int(seed_loc[0])
            cor_res_ind = int(seed_loc[1])
            cor_seq = coarsedb.coarse_seq_get(cor_seq_id)

            ext_cor_start = cor_res_ind + map_seed_size
            ext_org_start = current + map_seed_size
            if ext_cor_start + ext_seed_size >= cor_seq.Len():
                continue

            cseq_ext = cor_seq.residues[ext_cor_start:ext_cor_start + ext_seed_size]
            oseq_ext = residues[ext_org_start:ext_org_start + ext_seed_size]
            if cseq_ext != oseq_ext:
                continue

            cor_match, org_match = extend_match(
                cor_seq.residues[cor_res_ind:],
                residues[current:],
                db.GappedWindowSize, db.UngappedWindowSize,
                db.MatchKmerSize, db.ExtSeqIdThreshold,
                mem,
            )
            if len(org_match) < db.MinMatchLen:
                continue

            # Gap-free fast path: extend_match advances cor/org pointers in
            # lockstep through align_ungapped, so when no gapped step was taken
            # both halves are equal length and gap-free. NW on equal-length
            # inputs maximises score, not identity — and any gapped alignment
            # for equal-length inputs has identity <= the 1:1 pairing's (gaps
            # only dilute the denominator). So the 1:1 alignment is the
            # right alignment for the identity gate; skip the DP.
            if len(cor_match) == len(org_match):
                alignment = (cor_match, org_match)
            else:
                alignment = nw_align(cor_match, org_match, mem)
            id_pct = seq_identity(alignment[0], alignment[1])
            if id_pct < db.MatchSeqIdThreshold:
                continue

            changed = False
            if len(org_match) + db.MatchExtend >= olen - current:
                org_match = residues[current:]
                changed = True

            if current - last_match <= db.MatchExtend:
                end = current + len(org_match)
                org_match = residues[last_match:end]
                current = last_match
                changed = True

            if changed:
                alignment = nw_align(cor_match, org_match, mem)

            cor_start = cor_res_ind
            cor_end = cor_start + len(cor_match)
            org_start = current
            org_end = org_start + len(org_match)

            if org_start - last_match > 0:
                org_sub = org_seq.new_sub_sequence(last_match, current)
                _add_without_match(cseq, coarsedb, org_seq_id, org_sub)

            cseq.add(new_link_to_coarse(
                cor_seq_id, cor_start, cor_end, alignment,
            ))
            cor_seq.add_link(LinkToCompressed(
                org_seq_id, cor_start, cor_end,
            ))

            last_match = org_end
            current = org_end - 1
            break

    if olen - last_match > 0:
        org_sub = org_seq.new_sub_sequence(last_match, olen)
        _add_without_match(cseq, coarsedb, org_seq_id, org_sub)

    return cseq


def extend_match(
    cor_res: bytes,
    org_res: bytes,
    gapped_window_size: int,
    ungapped_window_size: int,
    kmer_size: int,
    id_threshold: int,
    mem: Memory,
):
    """Greedily extend a match using ungapped + gapped extension passes."""
    cor_match_len = 0
    org_match_len = 0
    while True:
        if cor_match_len == len(cor_res) or org_match_len == len(org_res):
            break

        match_len = align_ungapped(
            cor_res[cor_match_len:],
            org_res[org_match_len:],
            ungapped_window_size, kmer_size, id_threshold,
        )
        cor_match_len += match_len
        org_match_len += match_len

        gapped_cor = cor_res[
            cor_match_len:min(len(cor_res), cor_match_len + gapped_window_size)
        ]
        gapped_org = org_res[
            org_match_len:min(len(org_res), org_match_len + gapped_window_size)
        ]
        if not gapped_cor or not gapped_org:
            break
        alignment = nw_align(gapped_cor, gapped_org, mem)

        id_pct = seq_identity(alignment[0], alignment[1])
        if id_pct < id_threshold:
            break

        cor_match_len += align_len(alignment[0])
        org_match_len += align_len(alignment[1])

    return cor_res[:cor_match_len], org_res[:org_match_len]


def skip_low_complexity(seq: bytes, window_size: int, region_size: int) -> int:
    """Return the offset past a low-complexity region starting near `seq[0]`."""
    upto = min(len(seq), window_size + region_size)
    last = 0
    repeats = 1
    found = False
    i = 0
    while i < upto:
        if seq[i] == last:
            repeats += 1
            if repeats >= region_size:
                found = True
                break
            i += 1
            continue
        last = seq[i]
        repeats = 1
        i += 1

    if not found:
        return 0

    while i < len(seq):
        if seq[i] != last:
            break
        i += 1
    return i


def _add_without_match(
    cseq: CompressedSeq,
    coarsedb: "CoarseDB",
    org_seq_id: int,
    org_sub: OriginalSeq,
) -> None:
    sub_copy = bytes(org_sub.residues)
    cor_seq_id, cor_seq = coarsedb.add(sub_copy)
    cor_seq.add_link(LinkToCompressed(org_seq_id, 0, len(sub_copy)))
    cseq.add(new_link_to_coarse_no_diff(cor_seq_id, 0, len(sub_copy)))
