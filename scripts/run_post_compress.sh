#!/bin/bash
# Run after hs-cablastp compress completes.
# Usage: bash scripts/run_post_compress.sh [--hs-timed-out]
set -e

HS_TIMED_OUT=0
if [[ "$1" == "--hs-timed-out" ]]; then
    HS_TIMED_OUT=1
fi

echo "=== Post-compress pipeline starting at $(date) ==="

# Parse timings
CAB_TIMING=$(python3 scripts/parse_timing.py bench/cab_compress.timing 2>/dev/null)
CAB_WALL=$(echo "$CAB_TIMING" | grep wall | cut -d= -f2)
echo "cablastp wall time: ${CAB_WALL}s"

if [[ $HS_TIMED_OUT -eq 0 && -s bench/hs_compress.timing ]]; then
    HS_TIMING=$(python3 scripts/parse_timing.py bench/hs_compress.timing 2>/dev/null)
    HS_WALL=$(echo "$HS_TIMING" | grep wall | cut -d= -f2)
    echo "hs-cablastp wall time: ${HS_WALL}s"
    TIMED_OUT_FLAG=""
else
    HS_WALL="TIMEOUT"
    TIMED_OUT_FLAG="--hs-timed-out"
    echo "hs-cablastp TIMED OUT"
fi

# Write performance.csv
python3 scripts/write_performance_csv.py \
    --cab-timing bench/cab_compress.timing \
    --hs-timing bench/hs_compress.timing \
    --cab-db bench/cab_db \
    --hs-db bench/hs_db \
    $TIMED_OUT_FLAG \
    --out bench/performance.csv
echo "Written bench/performance.csv"

# Run ROC benchmark
if [[ $HS_TIMED_OUT -eq 0 ]]; then
    echo "Running full ROC benchmark..."
    python3 scripts/roc_benchmark.py \
        --fasta data/sprot_5k.fasta \
        --cab-db bench/cab_db \
        --hs-db bench/hs_db \
        --n-queries 50 \
        --out bench \
        --summary-csv bench/summary.csv \
        --roc-csv bench/roc_points.csv \
        --cab-compress-seconds "$CAB_WALL" \
        --hs-compress-seconds "$HS_WALL" 2>&1 | tee bench/roc_run.log
else
    echo "Running cablastp-only ROC (hs-cablastp timed out)..."
    python3 scripts/roc_cab_only.py \
        --fasta data/sprot_5k.fasta \
        --cab-db bench/cab_db \
        --n-queries 50 \
        --out bench \
        --summary-csv bench/summary.csv \
        --roc-csv bench/roc_points.csv \
        --cab-compress-seconds "$CAB_WALL" 2>&1 | tee bench/roc_run.log
fi

echo "=== ROC benchmark complete ==="

# Git commit and push
git add -f bench/performance.csv bench/summary.csv bench/roc_points.csv bench/roc.png \
    bench/cab_compress.log bench/cab_compress.timing \
    bench/hs_compress.log bench/hs_compress.timing \
    bench/roc_run.log 2>/dev/null || true
git add -f bench/missed_*.tsv bench/extra_*.tsv 2>/dev/null || true
git commit -m 'bench: cablastp vs hs-cablastp on SwissProt 5k subset (2-hr hs timeout)' || echo "Nothing to commit"
git push origin bench/sprot-2026-05-21-5k-v2 || {
    sleep 2
    git push origin bench/sprot-2026-05-21-5k-v2 || {
        sleep 4
        git push origin bench/sprot-2026-05-21-5k-v2
    }
}

echo "=== DONE at $(date) ==="
echo "Branch: bench/sprot-2026-05-21-5k-v2"
