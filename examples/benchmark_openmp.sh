#!/bin/bash
# OpenMP Performance Benchmark Script
# Compares performance with different thread counts

echo "========================================"
echo "OpenMP Performance Benchmark"
echo "========================================"
echo ""

# Configuration
CONFIG_FILE="examples/benchmark.cfg"
EXECUTABLE="./bin/dynearthsol2d"

# Check if executable exists
if [ ! -f "$EXECUTABLE" ]; then
    echo "Error: $EXECUTABLE not found. Please build first with ./build.sh"
    exit 1
fi

# Check if config exists
if [ ! -f "$CONFIG_FILE" ]; then
    echo "Error: $CONFIG_FILE not found."
    exit 1
fi

# Create results directory
RESULTS_DIR="benchmark_results/openmp_scaling_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$RESULTS_DIR"

echo "Results will be saved to: $RESULTS_DIR"
echo ""

# Thread counts to test
THREAD_COUNTS=(1 2 4 8)

# Run benchmarks
for THREADS in "${THREAD_COUNTS[@]}"; do
    echo "----------------------------------------"
    echo "Running with OMP_NUM_THREADS=$THREADS"
    echo "----------------------------------------"
    
    OUTPUT_FILE="$RESULTS_DIR/threads_${THREADS}.log"
    CSV_FILE="$RESULTS_DIR/threads_${THREADS}.csv"
    
    # Run simulation
    OMP_NUM_THREADS=$THREADS $EXECUTABLE $CONFIG_FILE > "$OUTPUT_FILE" 2>&1
    
    # Move CSV file if it exists
    if [ -f "benchmark_2d_timing.csv" ]; then
        mv benchmark_2d_timing.csv "$CSV_FILE"
    fi
    
    # Extract timing summary
    echo "Timing results:"
    grep -A 10 "Performance Summary" "$OUTPUT_FILE" || echo "No summary found"
    echo ""
done

# Generate comparison report
echo "========================================"
echo "Performance Comparison Summary"
echo "========================================"
echo ""

REPORT_FILE="$RESULTS_DIR/comparison_report.txt"
echo "OpenMP Scaling Analysis" > "$REPORT_FILE"
echo "Generated: $(date)" >> "$REPORT_FILE"
echo "========================================" >> "$REPORT_FILE"
echo "" >> "$REPORT_FILE"

for THREADS in "${THREAD_COUNTS[@]}"; do
    CSV_FILE="$RESULTS_DIR/threads_${THREADS}.csv"
    if [ -f "$CSV_FILE" ]; then
        echo "Threads: $THREADS" >> "$REPORT_FILE"
        cat "$CSV_FILE" >> "$REPORT_FILE"
        echo "" >> "$REPORT_FILE"
        
        # Display on console
        echo "Threads: $THREADS"
        cat "$CSV_FILE"
        echo ""
    fi
done

echo "Full report saved to: $REPORT_FILE"
echo ""
echo "To analyze speedup, compare the 'Total_s' values for each section."
echo "Speedup = Time(1 thread) / Time(N threads)"
