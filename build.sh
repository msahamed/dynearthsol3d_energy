#!/bin/bash
# DynEarthSol Enhanced Build Script
# This script provides a convenient CLI for configuring and building DynEarthSol
# while integrating with the CMake-generated Makefile

# Default configuration
BUILD_TYPE="Release"
DIMENSIONS="2D"
USE_OPENMP=ON
USE_ADAPT=OFF
USE_THERMAL_STRESS=OFF  # Example of how to add a new feature flag
CLEAN=0
DEEPCLEAN=0
RECONFIGURE=0

# Print help message
function show_help {
    echo "DynEarthSol Build Script"
    echo "Usage: $0 [options]"
    echo ""
    echo "Options:"
    echo "  --help            Show this help message"
    echo "  --debug           Build with debug symbols"
    echo "  --release         Build optimized release version (default)"
    echo "  --2d              Build 2D version (default)"
    echo "  --3d              Build 3D version"
    echo "  --openmp          Enable OpenMP parallelization (default)"
    echo "  --no-openmp       Disable OpenMP parallelization"
    echo "  --adapt           Enable mesh adaptation support (requires VTK)"
    echo "  --no-adapt        Disable mesh adaptation (default)"
    echo "  --thermal-stress  Enable thermal stress calculations"
    echo "  --no-thermal-stress Disable thermal stress calculations (default)"
    echo "  --clean           Clean build files"
    echo "  --deepclean       Remove all generated files and start fresh"
    echo "  --reconfigure     Force CMake reconfiguration"
    echo ""
    echo "Example:"
    echo "  $0 --3d --openmp --debug    # Build 3D debug version with OpenMP"
    echo "  $0 --clean                  # Clean build files"
    echo "  $0                          # Build with default settings"
    exit 0
}

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --help)
            show_help
            ;;
        --debug)
            BUILD_TYPE="Debug"
            shift
            ;;
        --release)
            BUILD_TYPE="Release"
            shift
            ;;
        --2d)
            DIMENSIONS="2D"
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
        --no-adapt)
            USE_ADAPT=OFF
            shift
            ;;
        --thermal-stress)
            USE_THERMAL_STRESS=ON
            shift
            ;;
        --no-thermal-stress)
            USE_THERMAL_STRESS=OFF
            shift
            ;;
        --clean)
            CLEAN=1
            shift
            ;;
        --deepclean)
            DEEPCLEAN=1
            shift
            ;;
        --reconfigure)
            RECONFIGURE=1
            shift
            ;;
        *)
            echo "Unknown option: $1"
            echo "Use --help to see available options"
            exit 1
            ;;
    esac
done

# Handle cleaning operations
if [ $DEEPCLEAN -eq 1 ]; then
    echo "Performing deep clean..."
    # Remove all CMake-generated files and directories
    rm -rf CMakeFiles bin lib
    rm -f CMakeCache.txt cmake_install.cmake Makefile
    rm -f */CMakeCache.txt */cmake_install.cmake */Makefile
    rm -f */*/CMakeCache.txt */*/cmake_install.cmake */*/Makefile
    echo "Deep clean completed."
    
    # If deep clean was requested without building, exit here
    if [ $CLEAN -eq 1 ] && [ $RECONFIGURE -eq 0 ]; then
        exit 0
    fi
elif [ $CLEAN -eq 1 ]; then
    echo "Cleaning build files..."
    if [ -f Makefile ]; then
        make clean
    else
        echo "No Makefile found. Skipping clean operation."
    fi
    
    # If just cleaning was requested, exit here
    if [ $RECONFIGURE -eq 0 ]; then
        exit 0
    fi
fi

# Configure with CMake
# Only reconfigure if:
# 1. Explicitly requested with --reconfigure
# 2. No Makefile exists yet
# 3. After a deep clean
if [ $RECONFIGURE -eq 1 ] || [ ! -f Makefile ] || [ $DEEPCLEAN -eq 1 ]; then
    echo "Configuring with CMake..."
    
    # Construct CMake command with all options
    CMAKE_CMD="cmake . -DCMAKE_BUILD_TYPE=$BUILD_TYPE"
    CMAKE_CMD+=" -DWITH_3D=$([ "$DIMENSIONS" == "3D" ] && echo "ON" || echo "OFF")"
    CMAKE_CMD+=" -DWITH_OPENMP=$USE_OPENMP"
    CMAKE_CMD+=" -DWITH_ADAPT=$USE_ADAPT"
    
    # Add any new feature flags here
    if [ "$USE_THERMAL_STRESS" == "ON" ]; then
        CMAKE_CMD+=" -DWITH_THERMAL_STRESS=ON"
    fi
    
    # Execute CMake command
    echo "Executing: $CMAKE_CMD"
    eval $CMAKE_CMD
    
    if [ $? -ne 0 ]; then
        echo "CMake configuration failed!"
        exit 1
    fi
fi

# Determine number of CPU cores for parallel build
if [ "$(uname)" = "Darwin" ]; then
    # macOS
    NUM_CORES=$(sysctl -n hw.ncpu 2>/dev/null || echo 4)
else
    # Linux
    NUM_CORES=$(nproc 2>/dev/null || echo 4)
fi

# Build with make
echo "Building with $NUM_CORES cores..."
make -j$NUM_CORES

# Print result based on exit code
if [ $? -eq 0 ]; then
    echo "======================="
    echo "Build successful!"
    echo "======================="
    echo "  Build type: $BUILD_TYPE"
    echo "  Dimensions: $DIMENSIONS"
    echo "  OpenMP: $USE_OPENMP"
    echo "  Adaptivity: $USE_ADAPT"
    if [ "$USE_THERMAL_STRESS" == "ON" ]; then
        echo "  Thermal Stress: Enabled"
    fi
    
    # Print run instructions
    echo ""
    echo "Ready to run simulations:"
    if [ "$DIMENSIONS" == "2D" ]; then
        echo "  ./bin/dynearthsol2d examples/shear_zone_2d.cfg  # For 2D strike-slip model"
        echo "  ./bin/dynearthsol2d examples/simple_test.cfg    # For simple 2D test"
    else
        echo "  ./bin/dynearthsol3d examples/simple_3d_model.cfg  # For 3D model"
    fi
    echo ""
    echo "Use --help for build options"
else
    echo "Build failed! Check error messages above."
    exit 1
fi