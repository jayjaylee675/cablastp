"""Profile a cablastp or hs-cablastp compress run.

Captures both CPU time (cProfile) and memory allocations (tracemalloc), then
prints the top hot spots and dumps the .prof file for later inspection
with snakeviz / pstats.

Usage:
    python scripts/profile_compress.py {cablastp|hs_cablastp} <db_dir> <fasta>...

Notes
-----
* For cablastp the run is forced to `-p 1 --quiet`. The worker pool uses a
  ThreadPoolExecutor, so cProfile (main-thread only) will under-count time
  spent inside the worker. Look at *cumulative* time in the main thread to
  see orchestration overhead; per-residue alignment cost lives in the worker.
* hs_cablastp is single-threaded, so its profile is complete.
"""
from __future__ import annotations

import cProfile
import os
import pstats
import shutil
import sys
import time
import tracemalloc


def main() -> int:
    if len(sys.argv) < 4:
        print(__doc__)
        return 1

    target = sys.argv[1]
    db_dir = sys.argv[2]
    fastas = sys.argv[3:]

    if os.path.exists(db_dir):
        shutil.rmtree(db_dir)

    if target == "cablastp":
        from cablastp.commands.compress import main as compress_main
        argv = [db_dir, *fastas, "-p", "1", "--quiet"]
    elif target == "hs_cablastp":
        from hs_cablastp.commands.compress import main as compress_main
        argv = [db_dir, *fastas]
    else:
        print(f"unknown target: {target!r}; use 'cablastp' or 'hs_cablastp'",
              file=sys.stderr)
        return 1

    profile_path = f"{db_dir}.prof"

    tracemalloc.start(5)
    profiler = cProfile.Profile()
    wall_start = time.perf_counter()
    profiler.enable()
    rc = 0
    try:
        rc = compress_main(argv) or 0
    finally:
        profiler.disable()
        wall = time.perf_counter() - wall_start
        snap = tracemalloc.take_snapshot()
        current, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        profiler.dump_stats(profile_path)

    print()
    print("=" * 72)
    print(f"Profiled: {target}  argv: {' '.join(argv)}")
    print(f"Exit code: {rc}")
    print(f"Wall time: {wall:.2f} s")
    print(f"Peak memory (tracemalloc): {peak / 1024 / 1024:.2f} MB"
          f"   (current at end: {current / 1024 / 1024:.2f} MB)")
    print(f"cProfile dump: {profile_path}")
    print("=" * 72)

    stats = pstats.Stats(profiler)
    stats.strip_dirs()

    print("\n--- TOP 25 BY CUMULATIVE TIME ---")
    stats.sort_stats("cumulative").print_stats(25)

    print("\n--- TOP 25 BY SELF TIME (tottime) ---")
    stats.sort_stats("tottime").print_stats(25)

    print("\n--- TOP 15 ALLOCATIONS (by line) ---")
    for stat in snap.statistics("lineno")[:15]:
        print(stat)

    print("\n--- TOP 10 ALLOCATIONS (by file) ---")
    for stat in snap.statistics("filename")[:10]:
        print(stat)

    return rc


if __name__ == "__main__":
    sys.exit(main())
