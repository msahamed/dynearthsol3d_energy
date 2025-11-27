#!/bin/bash
echo "=========================================="
echo "LARGE-SCALE OpenMP Benchmark"
echo "5000 steps, 400m resolution"
echo "=========================================="
echo ""

CONFIG="examples/shear_large_benchmark.cfg"
RESULTS_DIR="benchmark_results/large_scale_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$RESULTS_DIR"

echo "Starting benchmarks... (this will take a few minutes)"
echo ""

for THREADS in 1 2 4 8; do
    echo "----------------------------------------"
    echo "Running with $THREADS thread(s)..."
    echo "----------------------------------------"
    START=$(date +%s)
    
    OMP_NUM_THREADS=$THREADS ./bin/dynearthsol2d $CONFIG > "$RESULTS_DIR/threads_${THREADS}.log" 2>&1
    
    END=$(date +%s)
    ELAPSED=$((END - START))
    
    echo "✓ Completed in ${ELAPSED} seconds"
    
    # Move timing CSV
    if [ -f "shear_large_timing.csv" ]; then
        mv shear_large_timing.csv "$RESULTS_DIR/threads_${THREADS}.csv"
    fi
    
    # Clean up output directories
    rm -rf output/shear_large_*
    
    echo ""
done

echo "=========================================="
echo "RESULTS SUMMARY"
echo "=========================================="
echo ""

# Extract and display results
echo "Thread Count | Wall Time (s) | Speedup | Efficiency"
echo "-------------|---------------|---------|------------"

BASELINE=""
for THREADS in 1 2 4 8; do
    LOG="$RESULTS_DIR/threads_${THREADS}.log"
    if [ -f "$LOG" ]; then
        # Extract wall time from log
        WALL_TIME=$(grep "total$" "$LOG" | tail -1 | awk '{print $1}' | sed 's/s//')
        
        if [ -z "$BASELINE" ]; then
            BASELINE=$WALL_TIME
            SPEEDUP="1.00x"
            EFFICIENCY="100%"
        else
            SPEEDUP=$(echo "scale=2; $BASELINE / $WALL_TIME" | bc)
            EFF=$(echo "scale=1; ($SPEEDUP / $THREADS) * 100" | bc)
            SPEEDUP="${SPEEDUP}x"
            EFFICIENCY="${EFF}%"
        fi
        
        printf "%12s | %13s | %7s | %10s\n" "$THREADS" "$WALL_TIME" "$SPEEDUP" "$EFFICIENCY"
    fi
done

echo ""
echo "Detailed results saved to: $RESULTS_DIR"
echo "=========================================="
