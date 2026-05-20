#!/usr/bin/env python3
"""Run cablastp-compress + cablastp-search (or the hs_cablastp variant) on a
FASTA file and append timing, peak-memory, disk-footprint, and match stats to
a log file."""

import argparse
import ctypes
import datetime
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path


# ----- platform-specific peak RSS -----------------------------------------

def peak_rss_bytes(popen):
    """Peak resident-set memory used by the finished child, in bytes (or None)."""
    if sys.platform == "win32":
        try:
            from ctypes import wintypes

            class PROCESS_MEMORY_COUNTERS(ctypes.Structure):
                _fields_ = [
                    ("cb", wintypes.DWORD),
                    ("PageFaultCount", wintypes.DWORD),
                    ("PeakWorkingSetSize", ctypes.c_size_t),
                    ("WorkingSetSize", ctypes.c_size_t),
                    ("QuotaPeakPagedPoolUsage", ctypes.c_size_t),
                    ("QuotaPagedPoolUsage", ctypes.c_size_t),
                    ("QuotaPeakNonPagedPoolUsage", ctypes.c_size_t),
                    ("QuotaNonPagedPoolUsage", ctypes.c_size_t),
                    ("PagefileUsage", ctypes.c_size_t),
                    ("PeakPagefileUsage", ctypes.c_size_t),
                ]

            counters = PROCESS_MEMORY_COUNTERS()
            counters.cb = ctypes.sizeof(PROCESS_MEMORY_COUNTERS)
            psapi = ctypes.WinDLL("psapi")
            handle = wintypes.HANDLE(int(popen._handle))
            if psapi.GetProcessMemoryInfo(handle, ctypes.byref(counters), counters.cb):
                return int(counters.PeakWorkingSetSize)
        except Exception:
            return None
        return None
    try:
        import resource
        usage = resource.getrusage(resource.RUSAGE_CHILDREN)
        return usage.ru_maxrss * (1 if sys.platform == "darwin" else 1024)
    except Exception:
        return None


def fmt_bytes(n):
    if n is None:
        return "n/a"
    size = float(n)
    for unit in ("B", "KiB", "MiB", "GiB"):
        if size < 1024:
            return f"{size:.1f} {unit}"
        size /= 1024
    return f"{size:.1f} TiB"


def dir_size(path):
    total = 0
    for root, _dirs, files in os.walk(path):
        for f in files:
            try:
                total += os.path.getsize(os.path.join(root, f))
            except OSError:
                pass
    return total


# ----- input introspection ------------------------------------------------

def count_input_sequences(fasta_path):
    n = 0
    with open(fasta_path) as f:
        for line in f:
            if line.startswith(">"):
                n += 1
    return n


# ----- subprocess runner --------------------------------------------------

def run_and_measure(cmd):
    start = time.perf_counter()
    proc = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    )
    out, err = proc.communicate()
    elapsed = time.perf_counter() - start
    rss = peak_rss_bytes(proc)
    return proc.returncode, elapsed, rss, out, err


# ----- variant-specific helpers ------------------------------------------

class _Variant:
    name = ""
    compress_module = ""
    search_module = ""

    def compress_cmd(self, db_dir, fasta, makeblastdb, blastp):
        raise NotImplementedError

    def search_cmd(self, db_dir, fasta, makeblastdb, blastp):
        raise NotImplementedError

    def parse_compress(self, out, err):
        return {}

    def parse_search(self, out, err):
        return {}


class _CablastpVariant(_Variant):
    name = "cablastp"
    compress_module = "cablastp.commands.compress"
    search_module = "cablastp.commands.search"

    def compress_cmd(self, db_dir, fasta, makeblastdb, blastp):
        cmd = [sys.executable, "-m", self.compress_module, "--overwrite", "--quiet"]
        if makeblastdb:
            cmd += ["--makeblastdb", makeblastdb]
        cmd += [str(db_dir), str(fasta)]
        return cmd

    def search_cmd(self, db_dir, fasta, makeblastdb, blastp):
        cmd = [sys.executable, "-m", self.search_module, "--quiet"]
        if makeblastdb:
            cmd += ["--makeblastdb", makeblastdb]
        if blastp:
            cmd += ["--blastp", blastp]
        cmd += [str(db_dir), str(fasta)]
        return cmd

    def parse_compress(self, out, err):
        text = out + err
        m = re.search(r"added (\d+) sequences", text)
        return {"coarse seqs": int(m.group(1)) if m else None}

    def parse_search(self, out, err):
        subjects = len(re.findall(r"^> ", out, flags=re.MULTILINE))
        hsps = len(re.findall(r"^ Score = ", out, flags=re.MULTILINE))
        no_hits = "***** No hits found *****" in out
        return {
            "subject hits": 0 if no_hits else subjects,
            "HSPs": 0 if no_hits else hsps,
        }


class _HsCablastpVariant(_Variant):
    name = "hs_cablastp"
    compress_module = "hs_cablastp.commands.compress"
    search_module = "hs_cablastp.commands.search"

    def compress_cmd(self, db_dir, fasta, makeblastdb, blastp):
        cmd = [sys.executable, "-m", self.compress_module, "--overwrite"]
        if makeblastdb:
            cmd += ["--makeblastdb", makeblastdb]
        cmd += [str(db_dir), str(fasta)]
        return cmd

    def search_cmd(self, db_dir, fasta, makeblastdb, blastp):
        cmd = [sys.executable, "-m", self.search_module]
        if makeblastdb:
            cmd += ["--makeblastdb", makeblastdb]
        if blastp:
            cmd += ["--blastp", blastp]
        cmd += [str(db_dir), str(fasta)]
        return cmd

    def parse_compress(self, out, err):
        text = out + err
        m = re.search(r"Forest: (\d+) total nodes, (\d+) roots, (\d+) children", text)
        if not m:
            return {}
        return {
            "total nodes":   int(m.group(1)),
            "roots":         int(m.group(2)),
            "children":      int(m.group(3)),
        }

    def parse_search(self, out, err):
        result = {}
        for key, pat in [
            ("coarse HSPs",   r"coarse HSPs:\s+(\d+)"),
            ("matched roots", r"matched root nodes:\s+(\d+)"),
            ("candidates",    r"candidates after prune:\s*(\d+)"),
            ("fine HSPs",     r"fine HSPs:\s+(\d+)"),
        ]:
            m = re.search(pat, out)
            if m:
                result[key] = int(m.group(1))
        return result


_VARIANTS = {
    "cablastp":    _CablastpVariant(),
    "hs_cablastp": _HsCablastpVariant(),
}


# ----- main ---------------------------------------------------------------

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("fasta", help="Input FASTA (also used as the search query)")
    p.add_argument("--variant", choices=sorted(_VARIANTS.keys()),
                   default="cablastp", help="Which pipeline to run")
    p.add_argument("--log", default="cablastp_runs.log", help="Log file (appended to)")
    p.add_argument("--db-dir", default=None, help="Database dir (default: fresh temp dir)")
    p.add_argument("--makeblastdb", default=None, help="Path to makeblastdb executable")
    p.add_argument("--blastp", default=None, help="Path to blastp executable")
    p.add_argument("--keep-db", action="store_true", help="Keep the database after running")
    args = p.parse_args()

    variant = _VARIANTS[args.variant]
    fasta = Path(args.fasta).resolve()
    if not fasta.is_file():
        sys.exit(f"FASTA file not found: {fasta}")

    if args.db_dir:
        db_dir = Path(args.db_dir).resolve()
        db_dir.mkdir(parents=True, exist_ok=True)
        cleanup_db = False
    else:
        db_dir = Path(tempfile.mkdtemp(prefix=f"{variant.name}_db_"))
        cleanup_db = not args.keep_db

    compress_cmd = variant.compress_cmd(db_dir, fasta, args.makeblastdb, args.blastp)
    search_cmd = variant.search_cmd(db_dir, fasta, args.makeblastdb, args.blastp)

    input_seqs = count_input_sequences(fasta)
    file_size = fasta.stat().st_size
    timestamp = datetime.datetime.now().isoformat(timespec="seconds")

    print(f"== [{variant.name}] compress: {fasta.name} ==", flush=True)
    c_rc, c_time, c_rss, c_out, c_err = run_and_measure(compress_cmd)
    compress_stats = variant.parse_compress(c_out, c_err)

    db_size_after_compress = dir_size(db_dir) if c_rc == 0 else None

    print(f"== [{variant.name}] search: {fasta.name} ==", flush=True)
    s_rc, s_time, s_rss, s_out, s_err = run_and_measure(search_cmd)
    search_stats = variant.parse_search(s_out, s_err)

    log_path = Path(args.log).resolve()
    with open(log_path, "a", encoding="utf-8") as f:
        f.write("=" * 70 + "\n")
        f.write(f"timestamp:        {timestamp}\n")
        f.write(f"variant:          {variant.name}\n")
        f.write(f"input fasta:      {fasta}\n")
        f.write(f"file size:        {file_size} bytes\n")
        f.write(f"input sequences:  {input_seqs}\n")
        f.write(f"database dir:     {db_dir}\n")
        f.write("--- compress ---\n")
        f.write(f"  exit code:      {c_rc}\n")
        f.write(f"  wall time:      {c_time:.3f} s\n")
        f.write(f"  peak memory:    {fmt_bytes(c_rss)}\n")
        f.write(f"  db disk size:   {fmt_bytes(db_size_after_compress)}\n")
        for k, v in compress_stats.items():
            f.write(f"  {k+':':<16}{v}\n")
        if c_rc != 0:
            tail = (c_err or "").strip().splitlines()
            f.write(f"  stderr tail:    {tail[-1] if tail else ''}\n")
        f.write("--- search ---\n")
        f.write(f"  exit code:      {s_rc}\n")
        f.write(f"  wall time:      {s_time:.3f} s\n")
        f.write(f"  peak memory:    {fmt_bytes(s_rss)}\n")
        for k, v in search_stats.items():
            f.write(f"  {k+':':<16}{v}\n")
        if s_rc != 0:
            tail = (s_err or "").strip().splitlines()
            f.write(f"  stderr tail:    {tail[-1] if tail else ''}\n")
        f.write("\n")

    print(f"Log appended: {log_path}")
    if cleanup_db:
        shutil.rmtree(db_dir, ignore_errors=True)
    sys.exit(c_rc or s_rc)


if __name__ == "__main__":
    main()
