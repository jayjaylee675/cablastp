"""cablastp-compress CLI (port of cmd/cablastp-compress/main.go)."""

from __future__ import annotations

import argparse
import os
import shutil
import signal
import sys
import threading
import time
from typing import List, Optional, Set

from cablastp import misc
from cablastp.dbconf import default_db_conf
from cablastp.db import new_write_db, DB
from cablastp.fasta import read_original_seqs

from cablastp.commands._compress_core import CompressPool


# Residues NOT in BLOSUM62 are mapped to 'X'.
IGNORED_RESIDUES = (ord("J"), ord("O"), ord("U"))

INTERVAL = 1000


def _build_parser(default_conf) -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="cablastp-compress",
        usage=(
            "%(prog)s [flags] database-directory "
            "fasta-file [fasta-file ...]"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        add_help=True,
    )
    add = parser.add_argument

    add("--min-match-len", type=int, default=default_conf.MinMatchLen,
        dest="MinMatchLen", help="The minimum size of a match.")
    add("--match-kmer-size", type=int, default=default_conf.MatchKmerSize,
        dest="MatchKmerSize",
        help="The size of kmer fragments to match in ungapped extension.")
    add("--gapped-window-size", type=int, default=default_conf.GappedWindowSize,
        dest="GappedWindowSize", help="The size of the gapped match window.")
    add("--ungapped-window-size", type=int, default=default_conf.UngappedWindowSize,
        dest="UngappedWindowSize", help="The size of the ungapped match window.")
    add("--ext-seq-id-threshold", type=int, default=default_conf.ExtSeqIdThreshold,
        dest="ExtSeqIdThreshold",
        help="The sequence identity threshold of [un]gapped extension. "
             "(An integer in the inclusive range from 0 to 100.)")
    add("--match-seq-id-threshold", type=int, default=default_conf.MatchSeqIdThreshold,
        dest="MatchSeqIdThreshold",
        help="The sequence identity threshold of an entire match.")
    add("--match-extend", type=int, default=default_conf.MatchExtend,
        dest="MatchExtend",
        help="The maximum number of residues to blindly extend a match "
             "without regard to sequence identity.")
    add("--map-seed-size", type=int, default=default_conf.MapSeedSize,
        dest="MapSeedSize",
        help="The size of a seed in the K-mer map.")
    add("--ext-seed-size", type=int, default=default_conf.ExtSeedSize,
        dest="ExtSeedSize",
        help="The additional residues to require for each seed match.")
    add("--low-complexity", type=int, default=default_conf.LowComplexity,
        dest="LowComplexity",
        help="The window size used to detect regions of low complexity.")
    add("--seed-low-complexity", type=int, default=default_conf.SeedLowComplexity,
        dest="SeedLowComplexity",
        help="The seed window size used to detect regions of low complexity.")
    add("--plain", action="store_true", default=default_conf.SavePlain,
        dest="SavePlain",
        help="Also save plain-text versions of binary-encoded files.")
    add("--read-only", action="store_true", default=default_conf.ReadOnly,
        dest="ReadOnly", help="Create a read-only database (smaller, no append).")
    add("--makeblastdb", default=default_conf.BlastMakeBlastDB,
        dest="BlastMakeBlastDB",
        help="The location of the 'makeblastdb' executable.")

    add("-p", type=int, default=os.cpu_count() or 1, dest="num_workers",
        help="The maximum number of CPUs that can execute simultaneously.")
    add("--append", action="store_true", default=False,
        help="Append compressed sequences to an existing database.")
    add("--overwrite", action="store_true", default=False,
        help="Destroy any existing database before compressing.")
    add("--quiet", action="store_true", default=False,
        help="Only emit errors to stderr.")
    add("--max-seeds", type=float, default=8.0,
        help="When set, the in-memory seeds table is wiped when its memory "
             "usage exceeds this many gigabytes. 0 disables.")
    add("--cpuprofile", default="", help="Ignored (no Go pprof equivalent).")
    add("--memprofile", default="", help="Ignored (no Go pprof equivalent).")
    add("--memstats", default="", help="Path to write memstats to.")
    add("--mem-interval", action="store_true", default=False,
        help="Periodically write memstats during compression.")

    add("paths", nargs="+",
        help="database-directory followed by one or more fasta files.")
    return parser


def _explicit_flag_names(argv: List[str]) -> Set[str]:
    """Return the set of long-flag names actually present on the command line."""
    out: Set[str] = set()
    for arg in argv:
        if arg.startswith("--"):
            name = arg[2:].split("=", 1)[0]
            out.add(name)
    return out


def _apply_args_to_conf(args, conf) -> None:
    for field in (
        "MinMatchLen", "MatchKmerSize", "GappedWindowSize",
        "UngappedWindowSize", "ExtSeqIdThreshold", "MatchSeqIdThreshold",
        "MatchExtend", "MapSeedSize", "ExtSeedSize", "LowComplexity",
        "SeedLowComplexity", "SavePlain", "ReadOnly", "BlastMakeBlastDB",
    ):
        setattr(conf, field, getattr(args, field))


def _write_memstats(path: str) -> None:
    try:
        import resource
        usage = resource.getrusage(resource.RUSAGE_SELF)
        body = (
            "MaxRSS: %d KB\n"
            "UserTime: %.3fs\n"
            "SystemTime: %.3fs\n"
        ) % (usage.ru_maxrss, usage.ru_utime, usage.ru_stime)
    except ImportError:
        try:
            import psutil  # type: ignore
            mem = psutil.Process().memory_info()
            body = "RSS: %d KB\nVMS: %d KB\n" % (mem.rss // 1024, mem.vms // 1024)
        except ImportError:
            body = "(memstats unavailable on this platform)\n"
    with open(path, "w") as f:
        f.write(body)


def _verbose_output(args, org_seq_id: int, timer_state: dict) -> None:
    if org_seq_id % INTERVAL != 0:
        return
    if not args.quiet:
        elapsed = time.time() - timer_state["start"]
        seqs_per_sec = INTERVAL / elapsed if elapsed > 0 else 0.0
        sys.stdout.write(
            "\r%d sequences compressed (%.4f seqs/sec)"
            % (org_seq_id, seqs_per_sec)
        )
        sys.stdout.flush()
        timer_state["start"] = time.time()
    if args.mem_interval:
        if args.memstats:
            _write_memstats("%s.%d" % (args.memstats, org_seq_id))


def _cleanup(args, db: DB, pool: CompressPool) -> None:
    if args.memstats:
        _write_memstats("%s.last" % args.memstats)
    pool.done()
    try:
        db.save()
    except Exception as err:
        sys.stderr.write("Could not save database: %s\n" % err)
        sys.exit(1)
    db.write_close()


def main(argv: Optional[List[str]] = None) -> int:
    if argv is None:
        argv = sys.argv[1:]

    conf = default_db_conf()
    parser = _build_parser(conf)
    args = parser.parse_args(argv)

    if len(args.paths) < 2:
        parser.print_help(sys.stderr)
        return 1

    if args.append and args.overwrite:
        sys.stderr.write(
            "Both the 'append' and 'overwrite' flags are set. It does not "
            "make sense to set both of these flags.\n"
        )
        return 1

    if not args.quiet:
        misc.set_verbose(True)

    db_dir = args.paths[0]
    fasta_files = args.paths[1:]

    if args.overwrite and os.path.exists(db_dir):
        shutil.rmtree(db_dir)

    _apply_args_to_conf(args, conf)
    explicit = _explicit_flag_names(argv)

    db = new_write_db(args.append, conf, db_dir, explicit)
    misc.vprintln("")

    pool = CompressPool(db, args.num_workers)
    org_seq_id = db.com_db.num_sequences()

    main_quit = threading.Event()

    def signal_handler(signum, frame):
        main_quit.set()
        _cleanup(args, db, pool)
        os._exit(0)

    try:
        signal.signal(signal.SIGINT, signal_handler)
        if hasattr(signal, "SIGTERM"):
            signal.signal(signal.SIGTERM, signal_handler)
    except (ValueError, OSError):
        # Signals can't be set in non-main threads.
        pass

    timer_state = {"start": time.time()}

    for path in fasta_files:
        seq_queue = read_original_seqs(path, IGNORED_RESIDUES)
        if org_seq_id == 0:
            timer_state["start"] = time.time()
        while True:
            if main_quit.is_set():
                return 0
            item = seq_queue.get()
            if item.err is not None:
                sys.stderr.write("%s\n" % item.err)
                return 1
            if item.seq is None:
                break
            conf.BlastDBSize += item.seq.Len()
            org_seq_id = pool.compress(org_seq_id, item.seq)
            _verbose_output(args, org_seq_id, timer_state)
            if args.max_seeds > 0 and org_seq_id % 10000 == 0:
                db.coarse_db.seeds.maybe_wipe(args.max_seeds)

    misc.vprintln("\n")
    misc.vprintf("Wrote %s.\n", "compressed")
    misc.vprintf("Wrote %s.\n", "compressed.index")

    _cleanup(args, db, pool)
    return 0


if __name__ == "__main__":
    sys.exit(main())
