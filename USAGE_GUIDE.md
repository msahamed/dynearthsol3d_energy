# DynEarthSol3D Usage Guide

This guide provides step-by-step instructions for running DynEarthSol3D simulations and visualizing the results.

## 1. Building the Code

```bash
# Clean any previous builds
cd /Users/sabberahamed/Documents/dynearthsol3d_energy
make clean

# Build the code (2D version by default)
make
```

For debugging builds or different configurations:
```bash
# Debug build
make opt=0

# 3D build
# First edit Makefile and set ndims=3
make

# Without OpenMP
make openmp=0

# Deep clean (remove all built libraries)
make deepclean
```

## 2. Creating a Configuration File

Create a file named `sim_config.cfg` with the following content:

```
[mesh]
poly_filename = examples/rifting-2d.poly
meshing_option = 90
xlength = 500e3
ylength = 1
zlength = 150e3
resolution = 5e3

[sim]
max_steps = 100
output_step_interval = 10
output_averaged_fields = 0
modelname = my_simulation

[control]
characteristic_speed = 1e-10
is_quasi_static = 1
dt_fraction = 0.5

[bc]
vbc_x0 = 1
vbc_val_x0 = -1e-10
vbc_x1 = 1
vbc_val_x1 = 1e-10
vbc_z0 = 1
vbc_val_z0 = 0
vbc_z1 = 1
vbc_val_z1 = 0
surface_temperature = 273
mantle_temperature = 1600

[ic]
temperature_option = 0
oceanic_plate_age_in_yr = 60e6

[mat]
rheology_type = elasto-plastic
num_materials = 2
rho0 = [3210, 3300]
shear_modulus = [80.5e9, 80.5e9]
bulk_modulus = [128.2e9, 128.2e9]
```

## 3. Running the Simulation

```bash
# Run with the configuration file
./dynearthsol2d sim_config.cfg

# To see available options
./dynearthsol2d -h
```

Output files will be created with the prefix specified in `modelname` (e.g., `my_simulation.save.000000`, `my_simulation.info`, etc.)

## 4. Setting Up Python Environment for Visualization

```bash
# Create a Python virtual environment
python3 -m venv venv

# Activate the environment
source venv/bin/activate

# Install required packages
pip install numpy matplotlib

# For 3D visualization (optional)
pip install vtk
```

## 5. Converting Results to VTK Format

The provided `2vtk.py` script converts the binary output to VTK format, which can be read by visualization tools like ParaView.

```bash
# Convert all frames
python 2vtk.py my_simulation

# Convert specific frames
python 2vtk.py my_simulation -s 5 -e 10
```

This will create VTK files for each frame (e.g., `my_simulation.000005.vtk`).

## 6. Visualizing Results with ParaView

1. Install ParaView:
   ```bash
   brew install paraview
   ```

2. Open ParaView and:
   - Click "File" → "Open"
   - Navigate to your VTK files
   - Select one or more VTK files and click "OK"
   - Click "Apply" in the "Properties" panel
   - Choose visualization options in the "Display" panel

3. For time series:
   - Select all VTK files for the series
   - In "Properties" panel, check "Group files by: Smallest Sequence Number"
   - Click "Apply"
   - Use the animation controls to play through the sequence

## 7. Alternative: Simple Python Visualization

For quick visualization, you can use matplotlib:

```python
import matplotlib.pyplot as plt
import numpy as np
import os

def plot_info(prefix):
    """Plot basic info from the simulation"""
    info_file = prefix + '.info'
    if not os.path.exists(info_file):
        print(f"Error: {info_file} not found")
        return
        
    steps = []
    times = []
    
    with open(info_file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 4 and parts[1].isdigit():
                steps.append(int(parts[2]))
                times.append(float(parts[3]))
    
    plt.figure(figsize=(10, 6))
    plt.plot(steps, times, 'o-')
    plt.xlabel('Step')
    plt.ylabel('Simulation Time (s)')
    plt.title(f'Simulation Progress: {prefix}')
    plt.grid(True)
    plt.savefig(f"{prefix}_time_progress.png", dpi=300)
    plt.show()

# Example usage
plot_info('my_simulation')
```

## 8. Troubleshooting

- **Compilation errors**: Try `make opt=0` for a debug build with fewer optimizations
- **Binary format errors**: The binary format may vary with software versions. Check header formats with:
  ```bash
  hexdump -C -n 256 my_simulation.save.000000
  ```
- **Visualization errors**: Try generating dummy VTK files with domain bounds for debugging:
  ```bash
  python 2vtk.py my_simulation --dummy
  ```

## 9. Advanced Configuration

For more complex simulations, you can modify these parameters:

**Mesh refinement**:
```
[mesh]
meshing_option = 2
refined_zonex = [0.4, 0.6]
refined_zonez = [0.8, 1.0]
resolution = 2e3
```

**Rheology models**:
```
[mat]
rheology_type = elasto-visco-plastic
visc_exponent = [3.05, 3.05]
visc_coefficient = [1.25e-1, 1.25e-1]
visc_activation_energy = [3.76e5, 3.76e5]
```

**Boundary conditions**:
```
[bc]
vbc_z1 = 10
has_wrinkler_foundation = 1
wrinkler_delta_rho = 30
```

## 10. Future Updates

For upgrading the codebase:
1. Update to C++17 with CMake build system
2. Modernize visualization with Python VTK bindings 
3. Add proper test suites
4. Implement parallel solvers