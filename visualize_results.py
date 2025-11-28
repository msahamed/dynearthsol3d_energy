#!/usr/bin/env python3
"""
Simple visualization of DynEarthSol VTK output using matplotlib
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from matplotlib.colors import Normalize
import sys


def read_vtk(filename):
    """Read a VTK file and extract mesh and data"""
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    # Parse VTK file
    coords = []
    connectivity = []
    temperature = []
    plastic_strain = []
    von_mises = []
    
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        
        # Read points
        if line.startswith('POINTS'):
            npoints = int(line.split()[1])
            i += 1
            for j in range(npoints):
                parts = lines[i].split()
                coords.append([float(parts[0]), float(parts[1])])
                i += 1
            continue
        
        # Read cells
        if line.startswith('CELLS'):
            ncells = int(line.split()[1])
            i += 1
            for j in range(ncells):
                parts = lines[i].split()
                connectivity.append([int(parts[1]), int(parts[2]), int(parts[3])])
                i += 1
            continue
        
        # Read temperature
        if line.startswith('SCALARS temperature'):
            i += 2  # Skip LOOKUP_TABLE line
            for j in range(npoints):
                temperature.append(float(lines[i]))
                i += 1
            continue
        
        # Read plastic strain
        if line.startswith('SCALARS plastic_strain'):
            i += 2  # Skip LOOKUP_TABLE line
            for j in range(ncells):
                plastic_strain.append(float(lines[i]))
                i += 1
            continue
        
        # Read von Mises stress
        if line.startswith('SCALARS von_mises_stress'):
            i += 2  # Skip LOOKUP_TABLE line
            for j in range(ncells):
                von_mises.append(float(lines[i]))
                i += 1
            continue
        
        i += 1
    
    coords = np.array(coords)
    connectivity = np.array(connectivity)
    temperature = np.array(temperature) if temperature else None
    plastic_strain = np.array(plastic_strain) if plastic_strain else None
    von_mises = np.array(von_mises) if von_mises else None
    
    return coords, connectivity, temperature, plastic_strain, von_mises


def plot_results(vtk_files, output_file='visualization.png'):
    """Create a visualization of the simulation results"""
    nframes = len(vtk_files)
    
    fig, axes = plt.subplots(nframes, 3, figsize=(18, 5*nframes))
    if nframes == 1:
        axes = axes.reshape(1, -1)
    
    for idx, vtk_file in enumerate(vtk_files):
        print(f"Plotting {vtk_file}...")
        coords, connectivity, temperature, plastic_strain, von_mises = read_vtk(vtk_file)
        
        # Create triangulation
        triang = tri.Triangulation(coords[:, 0]/1000, coords[:, 1]/1000, connectivity)
        
        # Plot temperature
        ax = axes[idx, 0]
        if temperature is not None:
            tcf = ax.tricontourf(triang, temperature, levels=20, cmap='hot')
            plt.colorbar(tcf, ax=ax, label='Temperature (K)')
        ax.set_xlabel('X (km)')
        ax.set_ylabel('Z (km)')
        ax.set_title(f'Frame {idx}: Temperature')
        ax.set_aspect('equal')
        ax.grid(True, alpha=0.3)
        
        # Plot plastic strain
        ax = axes[idx, 1]
        if plastic_strain is not None:
            # Use tripcolor for cell data
            pcf = ax.tripcolor(triang, plastic_strain, cmap='YlOrRd', 
                              vmin=0, vmax=plastic_strain.max() if plastic_strain.max() > 0 else 1)
            plt.colorbar(pcf, ax=ax, label='Plastic Strain')
        ax.set_xlabel('X (km)')
        ax.set_ylabel('Z (km)')
        ax.set_title(f'Frame {idx}: Plastic Strain')
        ax.set_aspect('equal')
        ax.grid(True, alpha=0.3)
        
        # Plot von Mises stress
        ax = axes[idx, 2]
        if von_mises is not None:
            vmf = ax.tripcolor(triang, von_mises/1e6, cmap='viridis',
                              vmin=0, vmax=np.percentile(von_mises/1e6, 95))
            plt.colorbar(vmf, ax=ax, label='Von Mises Stress (MPa)')
        ax.set_xlabel('X (km)')
        ax.set_ylabel('Z (km)')
        ax.set_title(f'Frame {idx}: Von Mises Stress')
        ax.set_aspect('equal')
        ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"\nSaved visualization to {output_file}")
    plt.close()


def main():
    if len(sys.argv) < 2:
        print("Usage: python visualize_results.py <vtk_file1> [vtk_file2] ...")
        print("Example: python visualize_results.py output/*/runs/*.vtk")
        return 1
    
    vtk_files = sys.argv[1:]
    plot_results(vtk_files)
    return 0


if __name__ == '__main__':
    sys.exit(main())
