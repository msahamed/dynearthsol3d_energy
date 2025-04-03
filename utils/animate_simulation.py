#!/usr/bin/env python
"""
Animation script for DynEarthSol simulations.
Reads VTK files and generates animated visualizations of the simulation over time.
"""

import os
import sys
import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.tri as tri
from matplotlib.colors import Normalize
import vtk
from vtk.util.numpy_support import vtk_to_numpy

class VTKReader:
    """Class to read and extract data from VTK files"""
    
    def __init__(self, vtk_file):
        """Initialize with a vtk file path"""
        self.vtk_file = vtk_file
        self.reader = vtk.vtkUnstructuredGridReader()
        self.reader.SetFileName(self.vtk_file)
        self.reader.Update()
        self.data = self.reader.GetOutput()
        
    def get_points(self):
        """Extract point coordinates"""
        points = self.data.GetPoints()
        return vtk_to_numpy(points.GetData())
        
    def get_cells(self):
        """Extract cell connectivity"""
        cells = []
        for i in range(self.data.GetNumberOfCells()):
            cell = self.data.GetCell(i)
            cell_points = []
            for j in range(cell.GetNumberOfPoints()):
                cell_points.append(cell.GetPointId(j))
            cells.append(cell_points)
        return cells
    
    def get_point_data(self, name):
        """Extract point data array by name"""
        data_array = self.data.GetPointData().GetArray(name)
        if data_array:
            return vtk_to_numpy(data_array)
        return None
        
    def get_cell_data(self, name):
        """Extract cell data array by name"""
        data_array = self.data.GetCellData().GetArray(name)
        if data_array:
            return vtk_to_numpy(data_array)
        return None
    
    def get_available_arrays(self):
        """Get list of available data arrays"""
        point_arrays = []
        cell_arrays = []
        
        point_data = self.data.GetPointData()
        cell_data = self.data.GetCellData()
        
        for i in range(point_data.GetNumberOfArrays()):
            point_arrays.append(point_data.GetArrayName(i))
            
        for i in range(cell_data.GetNumberOfArrays()):
            cell_arrays.append(cell_data.GetArrayName(i))
            
        return {"point_data": point_arrays, "cell_data": cell_arrays}

def get_output_root():
    """Get the absolute path to the output directory"""
    # Always use project root's output directory
    project_root = '/Users/sabberahamed/Documents/dynearthsol3d_energy'
    return os.path.join(project_root, 'output')

def find_latest_run(model_name):
    """Find the latest run directory for a given model name"""
    output_root = get_output_root()
    # List all directories matching the model_name pattern
    matching_dirs = []
    for item in os.listdir(output_root):
        if item.startswith(model_name + "_") and os.path.isdir(os.path.join(output_root, item)):
            matching_dirs.append(item)
    
    if not matching_dirs:
        return None
    
    # Sort by timestamp (assuming format is model_YYYYMMDD_HHMMSS)
    matching_dirs.sort(reverse=True)
    return os.path.join(output_root, matching_dirs[0])

def animate_temperature(model_name, run_dir=None, fps=5, dpi=200):
    """Create an animation of temperature field evolving through time"""
    
    # Find the run directory if not specified
    if run_dir is None:
        run_dir = find_latest_run(model_name)
        if not run_dir:
            print(f"No runs found for model '{model_name}'")
            return
    
    # Find all VTK files for this model
    vtk_pattern = os.path.join(run_dir, "vtk", f"{model_name}_*.vtk")
    vtk_files = sorted(glob.glob(vtk_pattern))
    
    if not vtk_files:
        print(f"No VTK files found matching pattern: {vtk_pattern}")
        return
    
    print(f"Found {len(vtk_files)} VTK files for animation")
    
    # Initialize figure and axis
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Read first file to initialize
    vtk_reader = VTKReader(vtk_files[0])
    points = vtk_reader.get_points()
    cells = vtk_reader.get_cells()
    
    # Get X and Z coordinates (for 2D visualization)
    x = points[:, 0]
    z = points[:, 1] if points.shape[1] == 2 else points[:, 2]
    temperature = vtk_reader.get_point_data("temperature")
    
    # Create triangulation
    triang = tri.Triangulation(x, z, cells)
    
    # Find min and max temperature across all frames for consistent colormap
    temp_min = float('inf')
    temp_max = float('-inf')
    
    for vtk_file in vtk_files:
        reader = VTKReader(vtk_file)
        temp = reader.get_point_data("temperature")
        temp_min = min(temp_min, np.min(temp))
        temp_max = max(temp_max, np.max(temp))
    
    # Initialize the plot
    contour = ax.tricontourf(triang, temperature, cmap='rainbow', 
                             levels=20, norm=Normalize(temp_min, temp_max))
    plt.tight_layout()
    cbar = plt.colorbar(contour, ax=ax)
    cbar.set_label('Temperature (K)')
    
    # Set up plot properties
    ax.set_aspect('equal')
    ax.set_xlabel('X (m)')
    ax.set_ylabel('Z (m)')
    time_text = ax.text(0.05, 0.95, "", transform=ax.transAxes)
    frame_text = ax.text(0.05, 0.90, "", transform=ax.transAxes)
    ax.grid(True, alpha=0.3)
    
    # Animation update function
    def update_frame(frame_idx):
        # Clear previous contours
        for coll in ax.collections[:]:
            coll.remove()
        
        # Read data for current frame
        vtk_reader = VTKReader(vtk_files[frame_idx])
        temperature = vtk_reader.get_point_data("temperature")
        
        # Extract simulation time from filename
        frame_num = os.path.basename(vtk_files[frame_idx]).split('.')[0].split('_')[-1]
        frame_num = int(frame_num)
        
        # Update plot
        contour = ax.tricontourf(triang, temperature, cmap='rainbow', 
                                 levels=20, norm=Normalize(temp_min, temp_max))
        
        # Update title and text information
        time_text.set_text(f"Frame: {frame_num}")
        frame_text.set_text(f"File: {os.path.basename(vtk_files[frame_idx])}")
        
        return contour,
    
    # Create animation
    ani = animation.FuncAnimation(fig, update_frame, frames=len(vtk_files),
                                  interval=1000/fps, blit=False)
    
    # Save animation
    viz_dir = os.path.join(run_dir, "viz")
    os.makedirs(viz_dir, exist_ok=True)
    output_file = os.path.join(viz_dir, f"{model_name}_temperature_animation.mp4")
    ani.save(output_file, writer='ffmpeg', fps=fps, dpi=dpi)
    print(f"Animation saved to {output_file}")
    
    plt.close()

def animate_stress(model_name, run_dir=None, fps=5, dpi=200):
    """Create an animation of stress field evolving through time"""
    
    # Find the run directory if not specified
    if run_dir is None:
        run_dir = find_latest_run(model_name)
        if not run_dir:
            print(f"No runs found for model '{model_name}'")
            return
    
    # Find all VTK files for this model
    vtk_pattern = os.path.join(run_dir, "vtk", f"{model_name}_*.vtk")
    vtk_files = sorted(glob.glob(vtk_pattern))
    
    if not vtk_files:
        print(f"No VTK files found matching pattern: {vtk_pattern}")
        return
    
    print(f"Found {len(vtk_files)} VTK files for animation")
    
    # Initialize figure and axis
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Read first file to initialize
    vtk_reader = VTKReader(vtk_files[0])
    points = vtk_reader.get_points()
    cells = vtk_reader.get_cells()
    
    # Get X and Z coordinates (for 2D visualization)
    x = points[:, 0]
    z = points[:, 1] if points.shape[1] == 2 else points[:, 2]
    
    # For stress, we need to compute cell centers as stress is cell data
    cell_centers = []
    for cell in cells:
        # Average of cell vertex coordinates
        x_center = np.mean([points[p][0] for p in cell])
        z_center = np.mean([points[p][1 if points.shape[1] == 2 else 2] for p in cell])
        cell_centers.append([x_center, z_center])
    cell_centers = np.array(cell_centers)
    
    # Get stress data
    stress_rate = vtk_reader.get_cell_data("strain_rate")
    
    # Find min and max stress across all frames for consistent colormap
    stress_min = float('inf')
    stress_max = float('-inf')
    
    for vtk_file in vtk_files:
        reader = VTKReader(vtk_file)
        stress = reader.get_cell_data("strain_rate")
        if stress is not None:
            stress_min = min(stress_min, np.min(stress))
            stress_max = max(stress_max, np.max(stress))
    
    # Create triangulation for cell data (more complex for irregular mesh)
    # Here we'll use a scatter plot instead for simplicity
    scatter = ax.scatter(cell_centers[:, 0], cell_centers[:, 1], 
                         c=stress_rate, cmap='viridis', 
                         norm=Normalize(stress_min, stress_max),
                         s=10, alpha=0.8)
    
    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label('Strain Rate')
    
    # Set up plot properties
    ax.set_aspect('equal')
    ax.set_xlabel('X (m)')
    ax.set_ylabel('Z (m)')
    time_text = ax.text(0.05, 0.95, "", transform=ax.transAxes)
    frame_text = ax.text(0.05, 0.90, "", transform=ax.transAxes)
    ax.grid(True, alpha=0.3)
    
    # Animation update function
    def update_frame(frame_idx):
        # Read data for current frame
        vtk_reader = VTKReader(vtk_files[frame_idx])
        stress = vtk_reader.get_cell_data("strain_rate")
        
        # Extract simulation time from filename
        frame_num = os.path.basename(vtk_files[frame_idx]).split('.')[0].split('_')[-1]
        frame_num = int(frame_num)
        
        # Update scatter plot colors
        scatter.set_array(stress)
        
        # Update title and text information
        time_text.set_text(f"Frame: {frame_num}")
        frame_text.set_text(f"File: {os.path.basename(vtk_files[frame_idx])}")
        
        return scatter,
    
    # Create animation
    ani = animation.FuncAnimation(fig, update_frame, frames=len(vtk_files),
                                  interval=1000/fps, blit=False)
    
    # Save animation
    viz_dir = os.path.join(run_dir, "viz")
    os.makedirs(viz_dir, exist_ok=True)
    output_file = os.path.join(viz_dir, f"{model_name}_stress_animation.mp4")
    ani.save(output_file, writer='ffmpeg', fps=fps, dpi=dpi)
    print(f"Animation saved to {output_file}")
    
    plt.close()

def main():
    if len(sys.argv) < 2:
        print("Usage: python animate_simulation.py <model_name> [run_timestamp] [fps]")
        return
    
    model_name = sys.argv[1]
    
    # Check if a specific run timestamp was provided
    run_dir = None
    fps = 5
    
    if len(sys.argv) > 2:
        # Check if second argument is a timestamp (starts with digits)
        if sys.argv[2].startswith(('0', '1', '2', '3', '4', '5', '6', '7', '8', '9')):
            # Find the run with that timestamp
            timestamp = sys.argv[2]
            run_pattern = os.path.join(get_output_root(), f"{model_name}_{timestamp}*")
            matching_runs = glob.glob(run_pattern)
            if matching_runs:
                run_dir = matching_runs[0]
            else:
                print(f"No run found matching timestamp: {timestamp}")
                return
                
            # Get fps if provided
            fps = int(sys.argv[3]) if len(sys.argv) > 3 else 5
        else:
            # Second argument is fps
            fps = int(sys.argv[2])
    
    # Create temperature animation
    animate_temperature(model_name, run_dir, fps)
    
    # Create stress animation
    animate_stress(model_name, run_dir, fps)

if __name__ == "__main__":
    main()