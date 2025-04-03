# DynEarthSol Configuration Guide

This document provides a comprehensive guide to creating configuration files (.cfg) for DynEarthSol3D simulations.

## Configuration File Format

DynEarthSol3D uses the Boost Program Options library to parse configuration files. The format is simple:

```
parameter_name = value
```

For array values, use square brackets:

```
array_parameter = [value1, value2, value3]
```

Comments start with `#`:

```
# This is a comment
parameter_name = value  # This is also a comment
```

## Basic Configuration Example

Here's a minimal configuration file to run a simple 2D simulation:

```
# Basic simulation parameters
sim.modelname = simple_test
sim.max_steps = 1000
sim.output_step_interval = 10

# Mesh configuration
mesh.meshing_option = 1  # Rectangular box with uniform resolution
mesh.xlength = 100e3     # 100 km width
mesh.zlength = 50e3      # 50 km depth
mesh.resolution = 5e3    # 5 km resolution
mesh.min_angle = 32      # Minimum angle for triangles (2D only)

# Material properties
mat.rheology_type = elastic
mat.num_materials = 1
mat.rho0 = [3210]                   # Density (kg/m³)
mat.thermal_coefficient = [3e-5]    # Thermal expansion coefficient (1/K)
mat.bulk_modulus = [128.2e9]        # Bulk modulus (Pa)
mat.shear_modulus = [80.5e9]        # Shear modulus (Pa)
```

For 3D simulations, you need to add the y-dimension:

```
mesh.ylength = 100e3  # 100 km length in y-direction
```

## Configuration Sections

The configuration is divided into logical sections:

1. **sim** - Simulation parameters
2. **mesh** - Mesh generation and remeshing
3. **markers** - Marker (material point) settings
4. **control** - Physical and numerical control parameters
5. **bc** - Boundary conditions
6. **ic** - Initial conditions
7. **mat** - Material properties

## Detailed Parameter Reference

### Simulation Parameters (sim)

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| sim.modelname | string | *required* | Prefix for output files |
| sim.max_steps | int | *required* | Maximum number of time steps |
| sim.max_time_in_yr | double | *optional* | Maximum simulation time in years |
| sim.output_step_interval | int | *required* | Output frequency in steps |
| sim.output_time_interval_in_yr | double | *optional* | Output frequency in years |
| sim.checkpoint_frame_interval | int | 10 | Frequency of checkpoint writing |
| sim.is_restarting | bool | false | Whether to restart from a checkpoint |
| sim.restarting_from_modelname | string | *optional* | Model name for restarting |
| sim.restarting_from_frame | int | *optional* | Output frame for restarting |
| sim.has_initial_checkpoint | bool | false | Output checkpoint at 0th step |
| sim.has_marker_output | bool | false | Output marker coordinates and material |
| sim.has_output_during_remeshing | bool | false | Output before and after remeshing |
| sim.output_averaged_fields | int | 1 | Output time-averaged field variables (0=no, 1=yes, N>2=average over N steps) |

### Mesh Parameters (mesh)

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| mesh.meshing_option | int | 1 | Mesh generation method: 1=uniform box, 2=refined box, 90=poly file |
| mesh.xlength | double | *required* | Domain length in x-direction (m) |
| mesh.ylength | double | *required for 3D* | Domain length in y-direction (m) |
| mesh.zlength | double | *required* | Domain length in z-direction (m) |
| mesh.resolution | double | *required* | Target spatial resolution (m) |
| mesh.smallest_size | double | 0.01 | Size of smallest element relative to resolution |
| mesh.largest_size | double | 30 | Size of largest element relative to resolution |
| mesh.min_angle | double | 32 | Min. angle of triangles in degrees (2D only) |
| mesh.min_tet_angle | double | 22 | Min. dihedral angle of tetrahedra in degrees (3D only) |
| mesh.max_ratio | double | 2 | Max. radius/length ratio of tetrahedra (3D only) |
| mesh.quality_check_step_interval | int | 100 | How often to check mesh quality |
| mesh.min_quality | double | 0.4 | Min. mesh quality before remeshing (0-1) |
| mesh.max_boundary_distortion | double | 0.25 | Max. boundary distortion before remeshing |
| mesh.remeshing_option | int | 0 | How to handle boundaries during remeshing |
| mesh.refined_zonex | string | "[0.4, 0.6]" | Refinement zone in x (for meshing_option=2) |
| mesh.refined_zoney | string | "[0.4, 0.6]" | Refinement zone in y (for meshing_option=2, 3D only) |
| mesh.refined_zonez | string | "[0.8, 1]" | Refinement zone in z (for meshing_option=2) |
| mesh.poly_filename | string | "mesh.poly" | Filename for polygon input (for meshing_option=90) |
| mesh.is_discarding_internal_segments | bool | true | Discard internal segments after initial mesh |

### Marker Parameters (markers)

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| markers.init_marker_option | int | 1 | Marker generation: 1=random, 2=regular |
| markers.markers_per_element | int | 4 | Number of markers per element (for option 1) |
| markers.init_marker_spacing | double | 0.3 | Marker spacing relative to resolution (for option 2) |
| markers.min_num_markers_in_element | int | 3 | Threshold for marker replenishment |
| markers.replenishment_option | int | 2 | Material type for new markers: 0=always 0, 1=probabilistic, 2=nearest |
| markers.random_seed | uint | 1 | Random seed (0=use current time) |

### Control Parameters (control)

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| control.gravity | double | 10 | Gravity magnitude (m/s²) |
| control.characteristic_speed | double | 0 | Characteristic tectonic speed (m/s) |
| control.is_quasi_static | bool | true | Use quasi-static mode (inertial scaling) |
| control.dt_fraction | double | 1.0 | Fraction of max stable time step (0-1) |
| control.fixed_dt | double | 0 | Fixed time step in seconds (0=dynamic) |
| control.inertial_scaling | double | 1e5 | Scaling factor for inertial (quasi-static) |
| control.damping_factor | double | 0.8 | Force damping factor (0-1) |
| control.ref_pressure_option | int | 0 | Reference pressure method |
| control.surface_process_option | int | 0 | Surface process: 0=none, 1=diffusion |
| control.surface_diffusivity | double | 1e-6 | Surface diffusion coefficient (m²/s) |
| control.has_thermal_diffusion | bool | true | Enable thermal diffusion |
| control.has_hydration_processes | bool | false | Enable hydration processes |
| control.hydration_migration_speed | double | 3e-9 | Hydrous fluid migration speed (m/s) |

### Boundary Conditions (bc)

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| bc.surface_temperature | double | 273 | Surface temperature (K) |
| bc.mantle_temperature | double | 1600 | Mantle temperature (K) |
| bc.has_wrinkler_foundation | bool | true | Enable Wrinkler foundation at bottom |
| bc.wrinkler_delta_rho | double | 0 | Excess density of Wrinkler foundation (kg/m³) |
| bc.has_elastic_foundation | bool | false | Enable elastic foundation at bottom |
| bc.elastic_foundation_constant | double | 1e11 | Elastic foundation constant |
| bc.has_water_loading | bool | true | Apply water loading on submerged boundaries |
| bc.vbc_x0 | int | 1 | Boundary condition type (left/west) |
| bc.vbc_x1 | int | 1 | Boundary condition type (right/east) |
| bc.vbc_val_x0 | double | -1e-9 | BC value (left/west) - m/s for velocity, Pa for stress |
| bc.vbc_val_x1 | double | 1e-9 | BC value (right/east) - m/s for velocity, Pa for stress |
| bc.vbc_y0 | int | 0 | Boundary condition type (south, 3D only) |
| bc.vbc_y1 | int | 0 | Boundary condition type (north, 3D only) |
| bc.vbc_val_y0 | double | 0 | BC value (south) |
| bc.vbc_val_y1 | double | 0 | BC value (north) |
| bc.vbc_z0 | int | 0 | Boundary condition type (bottom) |
| bc.vbc_z1 | int | 0 | Boundary condition type (top) |
| bc.vbc_val_z0 | double | 0 | BC value (bottom) |
| bc.vbc_val_z1 | double | 0 | BC value (top) |

Boundary condition types:
- 0: All velocity components free
- 1: Normal component fixed, shear components free
- 2: Normal component free, shear components fixed at 0
- 3: Normal component fixed, shear components fixed at 0
- 4-7: Complex 3D boundary conditions (see source documentation)

### Initial Conditions (ic)

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| ic.mattype_option | int | 0 | Material type initialization: 0=region, 1=layered |
| ic.num_mattype_layers | int | *required if option=1* | Number of material layers |
| ic.layer_mattypes | int[] | *required if option=1* | Material types for each layer |
| ic.mattype_layer_depths | double[] | *required if option=1* | Depths of material interfaces |
| ic.weakzone_option | int | 1 | Weak zone type: 0=none, 1=planar, 2=ellipsoidal |
| ic.weakzone_plstrain | double | 0.1 | Initial plastic strain in weak zone |
| ic.weakzone_azimuth | double | 0 | Azimuth angle relative to +y (degrees) |
| ic.weakzone_inclination | double | 90 | Inclination angle (degrees) |
| ic.weakzone_halfwidth | double | 1.5 | Half-width in units of resolution |
| ic.weakzone_depth_min | double | 0 | Min depth (0-1 relative to zlength) |
| ic.weakzone_depth_max | double | 1 | Max depth (0-1 relative to zlength) |
| ic.weakzone_y_min | double | 0 | Min y position (0-1 relative to ylength) |
| ic.weakzone_y_max | double | 1 | Max y position (0-1 relative to ylength) |
| ic.weakzone_xcenter | double | 0.5 | Center x position (0-1 relative to xlength) |
| ic.weakzone_ycenter | double | 0.5 | Center y position (0-1 relative to ylength) |
| ic.weakzone_zcenter | double | 0.5 | Center z position (0-1 relative to zlength) |
| ic.weakzone_xsemi_axis | double | 1000 | X semi-axis length for ellipsoidal weak zone (m) |
| ic.weakzone_ysemi_axis | double | 1000 | Y semi-axis length for ellipsoidal weak zone (m) |
| ic.weakzone_zsemi_axis | double | 1000 | Z semi-axis length for ellipsoidal weak zone (m) |
| ic.temperature_option | int | 0 | Temperature initialization: 0=half-space cooling, 90=from file |
| ic.oceanic_plate_age_in_yr | double | 60e6 | Age of oceanic plate (years) |
| ic.Temp_filename | string | "Thermal.dat" | Temperature file (option=90) |
| ic.Nodes_filename | string | "Coord.dat" | Coordinates file (option=90) |
| ic.Connectivity_filename | string | "Connectivity.dat" | Connectivity file (option=90) |
| ic.isostasy_adjustment_time_in_yr | double | 0 | Time for isostasy adjustment (years) |

### Material Properties (mat)

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| mat.rheology_type | string | *required* | "elastic", "viscous", "maxwell", "elasto-plastic", or "elasto-visco-plastic" |
| mat.is_plane_strain | bool | false | Use plane-strain formulation (2D only) |
| mat.phase_change_option | int | 0 | Phase changes: 0=none, 1=subduction |
| mat.num_materials | int | 1 | Number of material types |
| mat.visc_min | double | 1e18 | Minimum viscosity (Pa·s) |
| mat.visc_max | double | 1e24 | Maximum viscosity (Pa·s) |
| mat.tension_max | double | 1e9 | Maximum tensile stress (Pa) |
| mat.therm_diff_max | double | 5e-6 | Maximum thermal diffusivity (m²/s) |
| mat.rho0 | double[] | [3210] | Reference densities at 0 Pa, 273 K (kg/m³) |
| mat.thermal_coefficient | double[] | [3e-5] | Thermal expansion coefficients (1/K) |
| mat.bulk_modulus | double[] | [128.2e9] | Bulk moduli (Pa) |
| mat.shear_modulus | double[] | [80.5e9] | Shear moduli (Pa) |
| mat.visc_exponent | double[] | [3.05] | Non-linear viscosity exponents |
| mat.visc_coefficient | double[] | [1.25e-1] | Non-linear viscosity coefficients |
| mat.visc_activation_energy | double[] | [3.76e5] | Activation energies (J/mol) |
| mat.heat_capacity | double[] | [1e3] | Heat capacities (J/kg/K) |
| mat.therm_cond | double[] | [3] | Thermal conductivities (W/m/K) |
| mat.pls0 | double[] | [0] | Plastic strain where weakening starts |
| mat.pls1 | double[] | [0.1] | Plastic strain where weakening saturates |
| mat.cohesion0 | double[] | [4e7] | Initial cohesion (Pa) |
| mat.cohesion1 | double[] | [4e6] | Weakened cohesion (Pa) |
| mat.friction_angle0 | double[] | [30] | Initial friction angle (degrees) |
| mat.friction_angle1 | double[] | [5] | Weakened friction angle (degrees) |
| mat.dilation_angle0 | double[] | [0] | Initial dilation angle (degrees) |
| mat.dilation_angle1 | double[] | [0] | Weakened dilation angle (degrees) |

## Example Configurations

### Simple 2D Elastic Model

```
# Simple 2D elastic model
sim.modelname = elastic_test
sim.max_steps = 1000
sim.output_step_interval = 10
sim.output_averaged_fields = 0

# Mesh
mesh.meshing_option = 1
mesh.xlength = 100e3
mesh.zlength = 50e3
mesh.resolution = 5e3
mesh.min_angle = 32

# Material
mat.rheology_type = elastic
mat.num_materials = 1
mat.rho0 = [3210]
mat.thermal_coefficient = [3e-5]
mat.bulk_modulus = [128.2e9]
mat.shear_modulus = [80.5e9]

# Boundary conditions - extension
bc.vbc_x0 = 1
bc.vbc_x1 = 1
bc.vbc_val_x0 = -1e-9  # 1 mm/yr compression
bc.vbc_val_x1 = 1e-9   # 1 mm/yr extension
```

### 2D Model with Weak Zone and Custom Boundary

```
# 2D model with custom geometry and weak zone
sim.modelname = rifting_model
sim.max_steps = 5000
sim.output_step_interval = 50
sim.output_averaged_fields = 0

# Mesh from polygon file
mesh.meshing_option = 90
mesh.poly_filename = rifting-2d.poly
mesh.xlength = 500e3
mesh.zlength = 150e3
mesh.resolution = 5e3
mesh.min_angle = 30

# Material
mat.rheology_type = elasto-visco-plastic
mat.num_materials = 2
mat.rho0 = [2800, 3300]
mat.thermal_coefficient = [3e-5, 2.5e-5]
mat.bulk_modulus = [50e9, 130e9] 
mat.shear_modulus = [30e9, 70e9]
mat.visc_coefficient = [1e-15, 1e-13]
mat.visc_exponent = [3.0, 3.5]
mat.visc_activation_energy = [3.5e5, 5.4e5]
mat.heat_capacity = [1000, 1200]
mat.therm_cond = [3.0, 3.5]
mat.pls0 = [0, 0]
mat.pls1 = [0.1, 0.1]
mat.cohesion0 = [2e7, 5e7]
mat.cohesion1 = [1e6, 2e7]
mat.friction_angle0 = [30, 30]
mat.friction_angle1 = [5, 15]

# Initial conditions
ic.temperature_option = 0
ic.oceanic_plate_age_in_yr = 50e6
ic.weakzone_option = 1
ic.weakzone_plstrain = 0.2
ic.weakzone_azimuth = 0
ic.weakzone_inclination = 90
ic.weakzone_halfwidth = 5
ic.weakzone_xcenter = 0.5

# Boundary conditions - rifting
bc.vbc_x0 = 1
bc.vbc_x1 = 1
bc.vbc_val_x0 = -5e-10  # 0.5 cm/yr divergence
bc.vbc_val_x1 = 5e-10   # 0.5 cm/yr divergence
bc.surface_temperature = 273
bc.mantle_temperature = 1600
```

## Running a Simulation

To run a simulation with your configuration file, use:

```
./dynearthsol2d your_config.cfg    # For 2D simulations
./dynearthsol3d your_config.cfg    # For 3D simulations
```

## Tips and Troubleshooting

1. Always specify required parameters. The most critical ones are:
   - sim.modelname
   - sim.max_steps
   - sim.output_step_interval
   - mesh dimensions and resolution
   - mat.rheology_type

2. When running 3D simulations, remember to include mesh.ylength.

3. For array parameters with multiple materials, ensure the array length matches mat.num_materials.

4. Avoid setting sim.output_step_interval < sim.output_averaged_fields as this will cause an error.

5. When using polygon files (mesh.meshing_option = 90), ensure the file exists and follows the correct format.

6. For very fine meshes, reduce the time step size with control.dt_fraction (e.g., 0.5) to ensure stability.

7. When restarting a simulation, ensure the checkpoint files exist and set all the necessary restarting parameters.