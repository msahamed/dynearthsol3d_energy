# DynEarthSol3D Example Configurations

This directory contains example configuration files and documentation for running DynEarthSol3D simulations.

## Configuration Guide

The comprehensive configuration guide `CONFIG_GUIDE.md` provides detailed information about all possible configuration parameters and their usage.

## Example Configurations

The following example configuration files showcase different simulation scenarios:

1. **simple_test.cfg** - Minimal configuration for a basic 2D elastic simulation
2. **extension_model.cfg** - Continental extension/rifting model
3. **subduction_model.cfg** - Oceanic plate subduction beneath continental lithosphere
4. **simple_3d_model.cfg** - Basic 3D box model with extension and compression
5. **shear_zone_2d.cfg** - Strike-slip fault/shear zone development
6. **thermal_convection.cfg** - Thermal convection in the mantle

## Polygon Files

The directory also includes polygon geometry files for complex model domains:

1. **rifting-2d.poly** - 2D geometry for rifting simulations

## Running Simulations

To run a simulation using one of these configurations:

```bash
# For 2D simulations
./dynearthsol2d examples/simple_test.cfg

# For 3D simulations
./dynearthsol3d examples/simple_3d_model.cfg
```

Simulation outputs will be stored in timestamp-based directories under `output/` in the format:
`output/modelname_YYYYMMDD_HHMMSS/`

Each output directory contains:
- `runs/` - Binary output files and checkpoint data
- `vtk/` - VTK visualization files
- `viz/` - Generated visualizations

## Visualization

To visualize a completed simulation, use the plotting scripts in the `utils/` directory:

```bash
# Generate domain visualization for the latest run
python utils/plot_domain.py modelname

# Specify a particular run by timestamp
python utils/plot_domain.py modelname 20250403_000230
```