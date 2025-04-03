#!/bin/bash
# Build script for DynEarthSol

# Default configuration
BUILD_TYPE="Release"
DIMENSIONS="2D"
USE_OPENMP=ON
USE_ADAPT=OFF
BUILD_DIR="build"

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --debug)
            BUILD_TYPE="Debug"
            shift
            ;;
        --3d)
            DIMENSIONS="3D"
            shift
            ;;
        --openmp)
            USE_OPENMP=ON
            shift
            ;;
        --no-openmp)
            USE_OPENMP=OFF
            shift
            ;;
        --adapt)
            USE_ADAPT=ON
            shift
            ;;
        --clean)
            echo "Cleaning build directory..."
            rm -rf $BUILD_DIR
            shift
            ;;
        *)
            echo "Unknown option: $1"
            echo "Usage: $0 [--debug] [--3d] [--openmp|--no-openmp] [--adapt] [--clean]"
            exit 1
            ;;
    esac
done

# Create build directory if it doesn't exist
mkdir -p $BUILD_DIR
cd $BUILD_DIR

# Configure with CMake
echo "Configuring with CMake..."
cmake .. \
    -DCMAKE_BUILD_TYPE=$BUILD_TYPE \
    -DWITH_3D=$([ "$DIMENSIONS" == "3D" ] && echo "ON" || echo "OFF") \
    -DWITH_OPENMP=$USE_OPENMP \
    -DWITH_ADAPT=$USE_ADAPT

# Build
echo "Building..."
# Determine number of CPU cores
if [ "$(uname)" = "Darwin" ]; then
    # macOS
    NUM_CORES=$(sysctl -n hw.ncpu 2>/dev/null || echo 4)
else
    # Linux
    NUM_CORES=$(nproc 2>/dev/null || echo 4)
fi

cmake --build . -- -j$NUM_CORES

# Return to original directory
cd ..

# Print success message
if [ $? -eq 0 ]; then
    echo "Build successful!"
    echo "  Build type: $BUILD_TYPE"
    echo "  Dimensions: $DIMENSIONS"
    echo "  OpenMP: $USE_OPENMP"
    echo "  Adaptivity: $USE_ADAPT"
    
    # Print run instructions
    if [ "$DIMENSIONS" == "2D" ]; then
        echo ""
        echo "To run, use: ./build/dynearthsol2d <config_file>"
        echo "Example: ./build/dynearthsol2d test.cfg"
    else
        echo ""
        echo "To run, use: ./build/dynearthsol3d <config_file>"
    fi
else
    echo "Build failed!"
fi