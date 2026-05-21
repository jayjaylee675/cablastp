"""Parse /usr/bin/time -v output to extract wall_seconds and peak_rss_mb."""
import re
import sys

def parse_timing(path):
    with open(path) as f:
        text = f.read()

    # Extract wall clock time (m:ss.ss or h:mm:ss formats)
    m = re.search(r'Elapsed \(wall clock\) time.*?: ([0-9:]+\.?[0-9]*)', text)
    if not m:
        return None, None
    elapsed_str = m.group(1)
    parts = elapsed_str.split(':')
    if len(parts) == 2:
        wall_seconds = float(parts[0]) * 60 + float(parts[1])
    elif len(parts) == 3:
        wall_seconds = float(parts[0]) * 3600 + float(parts[1]) * 60 + float(parts[2])
    else:
        return None, None

    # Extract peak RSS
    m = re.search(r'Maximum resident set size \(kbytes\): (\d+)', text)
    if not m:
        return wall_seconds, None
    peak_rss_mb = int(m.group(1)) / 1024.0

    return wall_seconds, peak_rss_mb

if __name__ == '__main__':
    path = sys.argv[1]
    ws, rss = parse_timing(path)
    print(f"wall_seconds={ws:.2f}")
    print(f"peak_rss_mb={rss:.2f}")
