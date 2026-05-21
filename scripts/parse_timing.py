"""Parse /usr/bin/time -v output and print wall_seconds and peak_rss_mb."""
import re
import sys

with open(sys.argv[1]) as f:
    text = f.read()

wall = re.search(r"Elapsed \(wall clock\) time.*?:\s*(?:(\d+):)?(\d+):(\d+\.?\d*)", text)
rss = re.search(r"Maximum resident set size \(kbytes\):\s*(\d+)", text)

if wall:
    h = int(wall.group(1) or 0)
    m = int(wall.group(2))
    s = float(wall.group(3))
    total = h * 3600 + m * 60 + s
    print(f"wall_seconds={total:.2f}")
else:
    print("wall_seconds=UNKNOWN")

if rss:
    kb = int(rss.group(1))
    print(f"peak_rss_mb={kb/1024:.2f}")
else:
    print("peak_rss_mb=UNKNOWN")
