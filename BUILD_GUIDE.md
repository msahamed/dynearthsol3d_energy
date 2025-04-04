# DynEarthSol Build and Usage Guide

This guide demonstrates how to build DynEarthSol from scratch and run simulations using the example configuration files.

## 1. Building DynEarthSol From Scratch

DynEarthSol uses CMake as its build system, which makes it easy to compile on different platforms. Follow these steps to build the code from scratch:

### Prerequisites

Make sure you have the following dependencies installed:

- C++ compiler (GCC 7+ or Clang 10+)
- CMake (version 3.14+)
- Boost libraries (specifically boost_program_options)
- Python 3 with NumPy (for visualization)
- Optional: VTK libraries (for mesh adaptation)

### Step-by-Step Build Process

1. **Clone the repository** (if you haven't already):
   ```bash
   git clone https://github.com/your-username/dynearthsol3d_energy.git
   cd dynearthsol3d_energy
   ```

2. **Use the build script** (recommended for most users):
   ```bash
   # For 2D simulation (default)
   ./build.sh

   # For 3D simulation
   ./build.sh --3d

   # For debug build
   ./build.sh --debug

   # For build with OpenMP parallelization (on by default)
   ./build.sh --openmp

   # For build without OpenMP
   ./build.sh --no-openmp

   # For build with mesh adaptation (requires VTK)
   ./build.sh --adapt

   # To clean before building
   ./build.sh --clean
   ```

3. **Manual build** (alternative approach):
   ```bash
   # Create and enter build directory
   mkdir -p build
   cd build

   # Configure with CMake
   cmake .. -DWITH_3D=OFF -DWITH_OPENMP=ON -DWITH_ADAPT=OFF -DWITH_DEBUG=OFF

   # Build with multiple cores
   make -j$(nproc)  # On Linux
   # or
   make -j$(sysctl -n hw.ncpu)  # On macOS
   ```

4. **Verify the build**:
   After successful compilation, you should have the executable in the `bin` directory:
   ```bash
   # For 2D build
   ls -l bin/dynearthsol2d

   # For 3D build
   ls -l bin/dynearthsol3d
   ```

## 2. Running a Simulation Using Example Configurations

DynEarthSol comes with several example configuration files in the `examples` directory. These are ready-to-use setups for different geological scenarios.

### Available Example Configurations

Here are the key example configurations you can use:

1. **`simple_test.cfg`**: Basic 2D simulation for testing
2. **`shear_zone_2d.cfg`**: Strike-slip fault simulation
3. **`extension_model.cfg`**: Continental extension/rifting simulation
4. **`subduction_model.cfg`**: Subduction zone simulation
5. **`thermal_convection.cfg`**: Thermal convection test
6. **`simple_3d_model.cfg`**: Simple 3D simulation example

### Running a Simulation with an Example Configuration

To run a simulation using one of the example configurations:

```bash
# Go to the project root directory
cd /path/to/dynearthsol3d_energy

# For 2D simulation
./bin/dynearthsol2d examples/shear_zone_2d.cfg

# For 3D simulation
./bin/dynearthsol3d examples/simple_3d_model.cfg
```

The simulation will start and create a timestamped output directory based on the model name specified in the configuration file. You'll see output like this:

```
Checking consistency of input parameters...
Initializing mesh and field data...
  Created directory: /path/to/dynearthsol3d_energy/output/shear_zone_2d_YYYYMMDD_HHMMSS
  Created directory: /path/to/dynearthsol3d_energy/output/shear_zone_2d_YYYYMMDD_HHMMSS/runs
  Created directory: /path/to/dynearthsol3d_energy/output/shear_zone_2d_YYYYMMDD_HHMMSS/vtk
  Created directory: /path/to/dynearthsol3d_energy/output/shear_zone_2d_YYYYMMDD_HHMMSS/viz
  Created output directories for model: shear_zone_2d (timestamp: YYYYMMDD_HHMMSS)
  Wrote VTK file: /path/to/dynearthsol3d_energy/output/shear_zone_2d_YYYYMMDD_HHMMSS/vtk/shear_zone_2d_000000.vtk
  Output # 0, step = 0, time = 0 yr, dt = 0.264073 yr.
Starting simulation...
  ...
  Output # 10, step = 100, time = 26.4073 yr, dt = 0.264073 yr.
Ending simulation.
```

### Using a Custom Configuration

You can also create a new configuration file by modifying an existing example:

```bash
# Copy an example config to modify
cp examples/shear_zone_2d.cfg my_custom_model.cfg

# Edit the configuration file
nano my_custom_model.cfg  # or use your preferred text editor

# Run with your custom configuration
./bin/dynearthsol2d my_custom_model.cfg
```

Make sure to change the `sim.modelname` parameter in your custom configuration to give your simulation a unique name.

## 3. Full End-to-End Workflow Example

Here's a complete workflow example, from building to visualization:

```bash
# 1. Clean build from scratch
./build.sh --clean

# 2. Run a simulation with an example config
./bin/dynearthsol2d examples/shear_zone_2d.cfg

# 3. Find the output directory (will have a timestamp)
latest_output=$(find output -name "shear_zone_2d_*" -type d | sort | tail -n 1)
echo "Latest output: $latest_output"

# 4. Convert output to VTK format if needed (typically done automatically)
python 2vtk.py $latest_output/runs/shear_zone_2d.save.* $latest_output/vtk/

# 5. Visualize using a simple Python script
python utils/plot_domain.py shear_zone_2d
```

## 4. Customizing Build Parameters

DynEarthSol build parameters can be customized to fit your needs:

| Parameter | Options | Description |
|-----------|---------|-------------|
| WITH_3D | ON/OFF | Build for 3D simulations (default: OFF) |
| WITH_OPENMP | ON/OFF | Enable OpenMP parallelization (default: ON) |
| WITH_ADAPT | ON/OFF | Enable mesh adaptation (requires VTK, default: OFF) |
| WITH_DEBUG | ON/OFF | Build with debug symbols (default: OFF) |

You can set these parameters either using the `build.sh` script or directly with CMake:

```bash
# Using build.sh
./build.sh --3d --no-openmp --debug

# Using CMake directly
cmake -B build -S . -DWITH_3D=ON -DWITH_OPENMP=OFF -DWITH_DEBUG=ON
cmake --build build -- -j$(nproc)
```

## 5. Output File Structure

Each simulation creates a timestamped output directory with the following structure:

```
output/
└── shear_zone_2d_YYYYMMDD_HHMMSS/
    ├── runs/              # Binary simulation data
    │   ├── shear_zone_2d.info
    │   ├── shear_zone_2d.chkpt.*
    │   └── shear_zone_2d.save.*
    ├── vtk/               # VTK files for visualization
    │   ├── shear_zone_2d_000000.vtk
    │   ├── shear_zone_2d_000001.vtk
    │   └── ...
    └── viz/               # Visualization outputs
        ├── shear_zone_2d_frame_*.png
        └── shear_zone_2d_time_series.png
```

## 6. Example Configuration: Understanding Key Parameters

Let's examine the key parameters in the example configuration files:

### Basic Simulation Parameters

```
# Simulation parameters
sim.modelname = shear_zone_2d    # Name of model (used for output files)
sim.max_steps = 100              # Maximum number of time steps
sim.output_step_interval = 10    # Save output every X steps
sim.max_time_in_yr = 1e5         # Maximum simulation time in years
```

### Mesh Parameters

```
# Mesh parameters
mesh.meshing_option = 1          # Mesh generation method (1 = box mesh)
mesh.xlength = 50e3              # Width of model domain in meters
mesh.ylength = 1e3               # Thickness (required for 2D) in meters
mesh.zlength = 20e3              # Depth of model domain in meters
mesh.resolution = 1e3            # Base resolution in meters
```

### Material Properties

```
# Material properties
mat.rheology_type = elasto-plastic   # Rheology model
mat.num_materials = 1                # Number of materials in model
mat.rho0 = [2800]                    # Density (kg/m³)
mat.bulk_modulus = [50e9]            # Bulk modulus (Pa)
mat.shear_modulus = [30e9]           # Shear modulus (Pa)
```

### Boundary Conditions

```
# Boundary conditions (example: strike-slip)
bc.vbc_x0 = 2                        # Left boundary condition type
bc.vbc_x1 = 2                        # Right boundary condition type
bc.vbc_val_x0 = 5e-9                 # Left boundary velocity (m/s)
bc.vbc_val_x1 = -5e-9                # Right boundary velocity (m/s)
bc.vbc_z0 = 3                        # Bottom boundary condition type
bc.vbc_z1 = 0                        # Top boundary condition type
```

For a complete parameter reference, see the `CONFIG_GUIDE.md` document.

## 7. Visualization Options

DynEarthSol outputs can be visualized in several ways:

### Basic Visualization with Included Scripts

```bash
# Plot domain overview
python utils/plot_domain.py shear_zone_2d

# Animate simulation progression
python utils/animate_simulation.py shear_zone_2d
```

### Advanced Visualization with ParaView

For more detailed visualization:

1. Install ParaView from [https://www.paraview.org/download/](https://www.paraview.org/download/)
2. Open the VTK files in ParaView:
   ```bash
   # Find VTK files
   find output -name "*.vtk"
   
   # Launch ParaView (varies by platform)
   paraview &
   ```
3. In ParaView: File → Open → Navigate to VTK files and select them

## 8. Troubleshooting Common Build Issues

If you encounter build problems, try these solutions:

- **Missing Boost libraries**: Install boost-program-options package for your system
  ```bash
  # Ubuntu/Debian
  sudo apt install libboost-program-options-dev
  
  # macOS
  brew install boost
  ```

- **Compilation errors**: Make sure your compiler supports C++17
  ```bash
  # Check GCC version
  g++ --version  # Should be 7.0 or higher
  
  # Check Clang version
  clang++ --version  # Should be 10.0 or higher
  ```

- **Build directory issues**: Try cleaning and rebuilding
  ```bash
  ./build.sh --clean
  ```

- **CMake configuration errors**: Make sure CMake is recent enough
  ```bash
  cmake --version  # Should be 3.14 or higher
  ```

## 9. Running Multiple Simulations for Parameter Studies

For parameter studies, create multiple configuration files with variations and run them sequentially:

```bash
# Create a script to run multiple simulations
cat > run_parameter_study.sh << 'EOF'
#!/bin/bash

# Array of config files to run
configs=(
  "examples/shear_zone_2d.cfg"
  "examples/extension_model.cfg"
  "examples/thermal_convection.cfg"
)

# Run each simulation
for config in "${configs[@]}"; do
  echo "Running simulation with config: $config"
  ./bin/dynearthsol2d "$config"
  echo "Completed simulation with config: $config"
  echo "--------------------------------------------"
done

echo "All simulations complete."
EOF

# Make the script executable
chmod +x run_parameter_study.sh

# Run the parameter study
./run_parameter_study.sh
```

## 10. Next Steps

Once you're comfortable with the basics, consider:

- Creating custom configuration files for your specific research questions
- Extending the code with new rheological models or boundary conditions
- Developing custom post-processing and visualization scripts
- Running larger, more complex 3D models
- Performing parameter sensitivity analyses

Refer to the `CONFIG_GUIDE.md` document for detailed information about all available configuration parameters.

## 11. Modifying Source Code and Rebuilding

If you modify the source code, rebuild the project:

```bash
# After making changes to source files
./build.sh

# Or to rebuild specific components
cd build
make dynearthsol2d
```

Remember to test your changes with simple configurations before running complex simulations.