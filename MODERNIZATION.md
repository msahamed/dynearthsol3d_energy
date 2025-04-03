# DynEarthSol3D Modernization Guide

This guide documents the modernization changes made to the DynEarthSol3D codebase to improve usability, maintainability, and performance.

## Overview of Changes

1. **Modern Build System**
   - Replaced Makefile with CMake for better cross-platform support
   - Added configurable build options with CMake
   - Created a simple build script (build.sh) for common build configurations

2. **Code Modernization**
   - Updated to C++17 standard
   - Better organization of code into subdirectories
   - Added proper unit testing framework

3. **Improved Visualization Tools**
   - Enhanced Python utilities for analysis and visualization
   - Added performance analysis tools 
   - Direct VTK export from the simulation
   - Organized output directory structure for visualization outputs

## Building the Code

### Prerequisites

- C++17 compiler (GCC 7+, Clang 5+, MSVC 2017+)
- CMake 3.14+
- Boost libraries (program_options)
- Python 3.6+ with NumPy and Matplotlib for visualization

### Using the Build Script

The simplest way to build the code is using the provided build script:

```bash
# Build 2D version (default)
./build.sh

# Build 3D version
./build.sh --3d

# Debug build
./build.sh --debug

# Without OpenMP
./build.sh --no-openmp

# With adaptivity
./build.sh --adapt

# Clean and rebuild
./build.sh --clean
```

### Manual CMake Build

You can also use CMake directly:

```bash
# Create build directory
mkdir build && cd build

# Configure
cmake .. -DCMAKE_BUILD_TYPE=Release -DWITH_3D=OFF -DWITH_OPENMP=ON

# Build
cmake --build .
```

## Running Simulations

The executable name depends on the dimension:

```bash
# For 2D simulations
./build/dynearthsol2d test.cfg

# For 3D simulations
./build/dynearthsol3d test.cfg
```

## Visualization and Analysis

The code now outputs VTK files directly during simulation, which can be found in the output directory. Several Python utilities are also provided for additional visualization and analysis:

```bash
# Visualize domain and temperature field
python plot_domain.py test_run

# Analyze performance metrics
python utils/analyze_performance.py test_run

# View VTK files directly in ParaView
# Open output/test_run_vtk/test_run_000010.vtk in ParaView
```

All visualization outputs are now organized in the `output/` directory, with VTK files in model-specific subdirectories.

## Directory Structure

The modernized codebase is organized as follows:

```
.
├── 3x3-C/               # 3x3 matrix operations library
├── ann/                 # Approximate Nearest Neighbor library
├── build/               # Build directory (created by build script)
├── examples/            # Example input configurations
├── libadaptivity/       # Mesh adaptation library (optional)
├── output/              # Organized output directory
│   └── <modelname>_vtk/ # VTK files for each simulation
├── tetgen/              # 3D mesh generation
├── tests/               # Unit tests
├── triangle/            # 2D mesh generation
├── utils/               # Utility scripts for visualization and analysis
├── build.sh             # Build script
├── CMakeLists.txt       # Main CMake configuration
└── README.md            # Documentation
```

## Future Improvements

Future modernization work could include:

1. Modularize the code further into component libraries
2. Add continuous integration with GitHub Actions
3. Improve parallel performance with MPI
4. Develop a GUI for simulation setup and visualization
5. Add HDF5 output format for better compatibility with analysis tools