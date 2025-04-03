# DynEarthSol End-to-End Usage Guide

This guide demonstrates the complete process of running a simulation with DynEarthSol3D, from configuration to visualization.

## 1. Setting Up Configuration

First, create a configuration file (`.cfg`) that defines your simulation parameters. Here's an example of a 2D shear zone model:

```
# 2D Shear Zone Model Configuration
# This simulates strike-slip faulting in a simple box

# Simulation parameters
sim.modelname = shear_zone_2d
sim.max_steps = 100
sim.output_step_interval = 10
sim.output_averaged_fields = 0
sim.max_time_in_yr = 1e5  # 100,000 years

# Mesh parameters
mesh.meshing_option = 1
mesh.xlength = 50e3      # 50 km width
mesh.ylength = 1e3       # 1 km thickness (required for 2D)
mesh.zlength = 20e3      # 20 km depth
mesh.resolution = 1e3    # 1 km resolution
mesh.min_angle = 30

# Material properties - elasto-plastic rheology
mat.rheology_type = elasto-plastic
mat.num_materials = 1
mat.rho0 = [2800]                   # Density (kg/m³)
mat.thermal_coefficient = [3e-5]    # Thermal expansion coefficient (1/K)
mat.bulk_modulus = [50e9]           # Bulk modulus (Pa)
mat.shear_modulus = [30e9]          # Shear modulus (Pa)
mat.pls0 = [0]                      # Plastic strain start
mat.pls1 = [0.1]                    # Plastic strain saturation
mat.cohesion0 = [2e7]               # Initial cohesion (Pa)
mat.cohesion1 = [1e6]               # Weakened cohesion (Pa)
mat.friction_angle0 = [30]          # Initial friction angle (degrees)
mat.friction_angle1 = [5]           # Weakened friction angle (degrees)

# Initial conditions - a central weak zone to localize deformation
ic.weakzone_option = 1
ic.weakzone_plstrain = 0.05
ic.weakzone_inclination = 90        # Vertical weak zone
ic.weakzone_halfwidth = 3           # Width in mesh resolution units
ic.weakzone_xcenter = 0.5           # Center of model
ic.weakzone_depth_min = 0           # From surface
ic.weakzone_depth_max = 1.0         # To bottom

# Control parameters
control.gravity = 9.81              # m/s²
control.is_quasi_static = true
control.dt_fraction = 0.5           # Reduced time step for stability

# Boundary conditions - strike-slip (y-direction velocity in a 2D x-z model)
bc.vbc_x0 = 2                       # Left boundary: fix shear, free normal
bc.vbc_x1 = 2                       # Right boundary: fix shear, free normal
bc.vbc_val_x0 = 5e-9                # 5 cm/yr right-lateral (into the plane)
bc.vbc_val_x1 = -5e-9               # 5 cm/yr right-lateral (out of the plane)
bc.vbc_z0 = 3                       # Fix bottom boundary (both normal and shear)
bc.vbc_z1 = 0                       # Free top boundary
bc.vbc_val_z0 = 0                   # No motion at bottom

# Thermal conditions
ic.temperature_option = 0           # Half-space cooling
bc.surface_temperature = 273        # Surface at 0°C
bc.mantle_temperature = 600         # Bottom at 327°C (shallow model)
```

Save this file in the `examples` directory as `shear_zone_2d.cfg`.

## 2. Building the Code

Navigate to the build directory and compile the code:

```bash
cd /path/to/dynearthsol3d_energy/build
make
```

This will build the necessary executables, including `dynearthsol2d` for 2D simulations and `dynearthsol3d` for 3D simulations.

## 3. Running the Simulation

Run the simulation using the configuration file:

```bash
cd /path/to/dynearthsol3d_energy/build
./bin/dynearthsol2d ../examples/shear_zone_2d.cfg
```

The simulation will run and produce output like this:

```
Checking consisitency of input parameters...
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

## 4. Output File Structure

The simulation creates a timestamped output directory with the following structure:

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

## 5. Basic Visualization

You can use the provided plotting script to visualize the simulation results:

```bash
cd /path/to/dynearthsol3d_energy/build
python bin/plot_domain.py shear_zone_2d
```

This generates two types of plots:
1. A frame plot showing the simulation domain at the latest frame
2. A time series plot showing the progression of simulation time vs. step

These plots are saved in the `viz/` directory with filenames like:
- `shear_zone_2d_frame_10.png`
- `shear_zone_2d_time_series.png`

## 6. Customized Visualization

For more specific visualizations, you can create custom scripts that read the VTK files. Here's an example script to visualize specific fields from the VTK output:

```python
#!/usr/bin/env python
"""
Script to plot fields from DynEarthSol VTK files
"""
import sys
import os
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy
import matplotlib.pyplot as plt
import matplotlib.tri as tri

def get_output_root():
    """Get the absolute path to the output directory"""
    project_root = '/path/to/dynearthsol3d_energy'
    return os.path.join(project_root, 'output')

def find_latest_run(model_name):
    """Find the latest run directory for a given model name"""
    output_root = get_output_root()
    matching_dirs = []
    for item in os.listdir(output_root):
        if item.startswith(model_name + "_") and os.path.isdir(os.path.join(output_root, item)):
            matching_dirs.append(item)
    
    if not matching_dirs:
        return None
    
    matching_dirs.sort(reverse=True)
    return os.path.join(output_root, matching_dirs[0])

def read_vtk_file(filename):
    """Read VTK file and extract data"""
    reader = vtk.vtkUnstructuredGridReader()
    reader.SetFileName(filename)
    reader.Update()
    
    data = reader.GetOutput()
    
    # Extract coordinates
    points = data.GetPoints()
    vertices = vtk_to_numpy(points.GetData())
    
    # Extract connectivity
    cells = data.GetCells()
    connectivity = []
    for i in range(data.GetNumberOfCells()):
        cell = data.GetCell(i)
        if cell.GetNumberOfPoints() == 3:  # Triangle cell
            connectivity.append([cell.GetPointId(0), cell.GetPointId(1), cell.GetPointId(2)])
    
    # Extract scalar data
    scalar_data = {}
    point_data = data.GetPointData()
    for i in range(point_data.GetNumberOfArrays()):
        name = point_data.GetArrayName(i)
        array = point_data.GetArray(i)
        scalar_data[name] = vtk_to_numpy(array)
    
    return {
        'coordinates': vertices,
        'connectivity': np.array(connectivity),
        'scalar_data': scalar_data
    }

def plot_scalar_field(data, field_name, output_path=None):
    """Plot a scalar field from the VTK data"""
    coords = data['coordinates']
    
    # For 2D, use only x and z components
    x = coords[:, 0]
    z = coords[:, 2] if coords.shape[1] > 2 else coords[:, 1]
    
    # Create triangulation
    triang = tri.Triangulation(x, z, data['connectivity'])
    
    # Get the scalar field data
    if field_name in data['scalar_data']:
        scalar_field = data['scalar_data'][field_name]
    else:
        print(f"Field '{field_name}' not found. Available fields: {list(data['scalar_data'].keys())}")
        return
    
    # Create figure
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Plot triangulation with scalar field
    tpc = ax.tripcolor(triang, scalar_field, cmap='viridis')
    
    # Add colorbar
    cbar = plt.colorbar(tpc, ax=ax)
    cbar.set_label(field_name)
    
    # Set labels and title
    ax.set_xlabel('X (m)')
    ax.set_ylabel('Z (m)')
    ax.set_title(f"DynEarthSol Simulation - {field_name.replace('_', ' ').title()}")
    
    # Equal aspect ratio
    ax.set_aspect('equal')
    
    # Save figure if output path provided
    if output_path:
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        print(f"Plot saved to {output_path}")
    else:
        plt.show()
    
    plt.close()

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python plot_fields.py <model_name> [frame_number] [field_name]")
        sys.exit(1)
    
    model_name = sys.argv[1]
    frame_number = int(sys.argv[2]) if len(sys.argv) > 2 else 10  # Default to frame 10
    field_name = sys.argv[3] if len(sys.argv) > 3 else 'temperature'  # Default to temperature
    
    # Find the latest run directory
    run_dir = find_latest_run(model_name)
    if not run_dir:
        print(f"No run directory found for model: {model_name}")
        sys.exit(1)
    
    # Find VTK file for the requested frame
    vtk_dir = os.path.join(run_dir, 'vtk')
    vtk_filename = os.path.join(vtk_dir, f"{model_name}_{frame_number:06d}.vtk")
    
    if not os.path.exists(vtk_filename):
        print(f"VTK file not found: {vtk_filename}")
        sys.exit(1)
    
    # Read VTK file
    data = read_vtk_file(vtk_filename)
    
    # Create output directory for custom plots
    viz_dir = os.path.join(run_dir, 'viz')
    os.makedirs(viz_dir, exist_ok=True)
    
    # Plot the scalar field
    output_path = os.path.join(viz_dir, f"{model_name}_{field_name}_frame_{frame_number:02d}.png")
    plot_scalar_field(data, field_name, output_path)
```

Run this script to visualize specific fields from the simulation:

```bash
python plot_fields.py shear_zone_2d 10 temperature
python plot_fields.py shear_zone_2d 10 strain_rate
```

## 7. Working with VTK Files Directly

You can inspect VTK files to understand what data is available:

```bash
# Create a simple script to inspect VTK files
cat > inspect_vtk.py << 'EOF'
#!/usr/bin/env python
import sys
import os

def read_vtk_header(filename):
    """Read and display information about a VTK file"""
    info = {}
    with open(filename, 'r') as f:
        lines = f.readlines()
        
        # Extract header information
        info['version'] = lines[0].strip()
        info['title'] = lines[1].strip()
        info['format'] = lines[2].strip()
        info['dataset'] = lines[3].strip()
        
        # Find POINTS line
        for i, line in enumerate(lines):
            if line.startswith('POINTS'):
                parts = line.strip().split()
                info['num_points'] = int(parts[1])
                break
                
        # Find CELLS line
        for i, line in enumerate(lines):
            if line.startswith('CELLS'):
                parts = line.strip().split()
                info['num_cells'] = int(parts[1])
                break
                
        # Find scalar data
        scalar_fields = []
        for i, line in enumerate(lines):
            if line.startswith('SCALARS'):
                parts = line.strip().split()
                scalar_fields.append(parts[1])
        
        info['scalar_fields'] = scalar_fields
        
    return info

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python inspect_vtk.py <vtk_file>")
        sys.exit(1)
        
    filename = sys.argv[1]
    if not os.path.exists(filename):
        print(f"Error: File {filename} not found.")
        sys.exit(1)
        
    info = read_vtk_header(filename)
    
    print(f"VTK File Information: {os.path.basename(filename)}")
    print(f"Version: {info['version']}")
    print(f"Title: {info['title']}")
    print(f"Format: {info['format']}")
    print(f"Dataset: {info['dataset']}")
    print(f"Number of Points: {info['num_points']}")
    print(f"Number of Cells: {info['num_cells']}")
    print("Scalar Fields:", ", ".join(info['scalar_fields']))
EOF

# Run the script on a VTK file
python inspect_vtk.py /path/to/dynearthsol3d_energy/output/shear_zone_2d_YYYYMMDD_HHMMSS/vtk/shear_zone_2d_000010.vtk
```

This will show you information about the VTK file and the available data fields:

```
VTK File Information: shear_zone_2d_000010.vtk
Version: # vtk DataFile Version 3.0
Title: DynEarthSol output, frame 10, time 26.4073 years
Format: ASCII
Dataset: DATASET UNSTRUCTURED_GRID
Number of Points: 583
Number of Cells: 1077
Scalar Fields: temperature, velocity, strain_rate, stress
```

## 8. Advanced Visualization

For more advanced visualization, you can use tools like ParaView to load and analyze the VTK files. ParaView provides a graphical interface for interactive visualization and analysis of scientific data.

1. Install ParaView from [https://www.paraview.org/download/](https://www.paraview.org/download/)
2. Launch ParaView and open a VTK file:
   - File → Open → Navigate to the vtk directory and select a file
   - Click "Apply" to load the data
3. Visualize different fields:
   - Select a field from the dropdown in the Properties panel
   - Apply color maps, filters, and other visualization techniques

## 9. Running Multiple Simulations

One of the key features of DynEarthSol3D is the ability to run multiple simulations without overwriting previous results. Each simulation creates a unique timestamped directory:

```
output/
├── shear_zone_2d_20250403_093651/  # First run
├── shear_zone_2d_20250403_102543/  # Second run
└── shear_zone_2d_20250403_115824/  # Third run
```

This allows you to easily compare results from different parameter sets or model configurations.

## 10. Next Steps

Once you're comfortable with the basics, consider exploring:

- Longer simulation runs with more time steps
- Different model configurations (extension, compression, subduction)
- 3D simulations using `dynearthsol3d`
- Custom post-processing scripts for advanced analysis
- Parameter studies by running multiple simulations with varying parameters

Refer to the `CONFIG_GUIDE.md` document for detailed information about all available configuration parameters.