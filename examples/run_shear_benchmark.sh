#!/bin/bash
echo "=========================================="
echo "Shear Zone OpenMP Scaling Benchmark"
echo "Configuration: 500 steps, resolution=800"
echo "=========================================="
echo ""

CONFIG="examples/shear_zone_benchmark.cfg"
RESULTS_DIR="benchmark_results/shear_zone_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$RESULTS_DIR"

for THREADS in 1 2 4 8; do
    echo "Running with $THREADS thread(s)..."
    START=$(date +%s)
    
    OMP_NUM_THREADS=$THREADS ./bin/dynearthsol2d $CONFIG > "$RESULTS_DIR/threads_${THREADS}.log" 2>&1
    
    END=$(date +%s)
    ELAPSED=$((END - START))
    
    echo "  Completed in ${ELAPSED} seconds"
    
    # Move timing CSV
    if [ -f "shear_zone_benchmark_timing.csv" ]; then
        mv shear_zone_benchmark_timing.csv "$RESULTS_DIR/threads_${THREADS}.csv"
    fi
    
    # Clean up output directories
    rm -rf output/shear_zone_benchmark_*
    
    echo ""
done

echo "=========================================="
echo "Results Summary"
echo "=========================================="
echo ""

# Extract timing data
echo "Thread Count | Total Time (s) | Speedup"
echo "-------------|----------------|--------"

BASELINE=""
for THREADS in 1 2 4 8; do
    CSV="$RESULTS_DIR/threads_${THREADS}.csv"
    if [ -f "$CSV" ]; then
        # Get Time Step total time
        TIME=$(grep "Time Step" "$CSV" | cut -d',' -f3)
        
        if [ -z "$BASELINE" ]; then
            BASELINE=$TIME
            SPEEDUP="1.00x"
        else
            SPEEDUP=$(echo "scale=2; $BASELINE / $TIME" | bc)
            SPEEDUP="${SPEEDUP}x"
        fi
        
        printf "%12s | %14.3f | %s\n" "$THREADS" "$TIME" "$SPEEDUP"
    fi
done

echo ""
echo "Full results saved to: $RESULTS_DIR"
