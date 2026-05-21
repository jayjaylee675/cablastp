"""Parse timing files and write bench/performance.csv."""
import argparse
import re
import subprocess
import sys
from pathlib import Path


def parse_timing(path: Path):
    try:
        text = path.read_text()
    except FileNotFoundError:
        return None, None

    wall = re.search(r"Elapsed \(wall clock\) time.*?:\s*(?:(\d+):)?(\d+):(\d+\.?\d*)", text)
    rss = re.search(r"Maximum resident set size \(kbytes\):\s*(\d+)", text)

    wall_s = None
    if wall:
        h = int(wall.group(1) or 0)
        m = int(wall.group(2))
        s = float(wall.group(3))
        wall_s = h * 3600 + m * 60 + s

    rss_mb = None
    if rss:
        rss_mb = int(rss.group(1)) / 1024.0

    return wall_s, rss_mb


def dir_size_bytes(path: Path) -> int:
    try:
        result = subprocess.run(
            ["du", "-sb", str(path)], capture_output=True, text=True
        )
        if result.returncode == 0:
            return int(result.stdout.split()[0])
    except Exception:
        pass
    total = 0
    if path.exists():
        for f in path.rglob("*"):
            if f.is_file():
                try:
                    total += f.stat().st_size
                except OSError:
                    pass
    return total


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cab-timing", type=Path, default=Path("bench/cab_compress.timing"))
    ap.add_argument("--hs-timing", type=Path, default=Path("bench/hs_compress.timing"))
    ap.add_argument("--cab-db", type=Path, default=Path("bench/cab_db"))
    ap.add_argument("--hs-db", type=Path, default=Path("bench/hs_db"))
    ap.add_argument("--hs-timed-out", action="store_true")
    ap.add_argument("--out", type=Path, default=Path("bench/performance.csv"))
    args = ap.parse_args()

    cab_wall, cab_rss = parse_timing(args.cab_timing)
    hs_wall, hs_rss = parse_timing(args.hs_timing)

    cab_db_bytes = dir_size_bytes(args.cab_db)
    hs_db_bytes = dir_size_bytes(args.hs_db)

    lines = ["pipeline,dataset,n_sequences,compress_wall_seconds,compress_peak_rss_mb,db_size_bytes"]

    def fmt(v, timeout=False):
        if timeout:
            return "TIMEOUT"
        return f"{v:.2f}" if v is not None else "UNKNOWN"

    lines.append(
        f"cablastp,sprot_5k,5000,"
        f"{fmt(cab_wall)},"
        f"{fmt(cab_rss)},"
        f"{cab_db_bytes}"
    )
    lines.append(
        f"hs-cablastp,sprot_5k,5000,"
        f"{fmt(hs_wall, args.hs_timed_out)},"
        f"{fmt(hs_rss, args.hs_timed_out)},"
        f"{hs_db_bytes if not args.hs_timed_out else 'TIMEOUT'}"
    )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text("\n".join(lines) + "\n")
    print(f"Written {args.out}")
    for line in lines:
        print(" ", line)


if __name__ == "__main__":
    main()
