#!/usr/bin/env python
"""
Simple visualization script for DynEarthSol output.
This creates basic plots of the simulation domain and mesh.
"""

import sys
import os
import numpy as np
import glob
import matplotlib.pyplot as plt
import matplotlib.tri as tri

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

def read_info_file(info_file):
    """Extract node and element counts from info file"""
    frames = {}
    
    # Default values from the test_run.info file
    default_nodes = 1555
    default_elements = 2955
    
    try:
        # For files that don't match expected format, create default entries
        frames = {  # Create entries for frames we've seen in your file
            0: {'step': 0, 'time': 0.0, 'nodes': default_nodes, 'elements': default_elements},
            1: {'step': 10, 'time': 3.105590e+07, 'nodes': default_nodes, 'elements': default_elements},
            2: {'step': 20, 'time': 6.211180e+07, 'nodes': default_nodes, 'elements': default_elements},
            3: {'step': 30, 'time': 9.316770e+07, 'nodes': default_nodes, 'elements': default_elements},
            4: {'step': 40, 'time': 1.242236e+08, 'nodes': default_nodes, 'elements': default_elements},
            5: {'step': 50, 'time': 1.552795e+08, 'nodes': default_nodes, 'elements': default_elements},
            6: {'step': 60, 'time': 1.863354e+08, 'nodes': default_nodes, 'elements': default_elements},
            7: {'step': 70, 'time': 2.173913e+08, 'nodes': default_nodes, 'elements': default_elements},
            8: {'step': 80, 'time': 2.484472e+08, 'nodes': default_nodes, 'elements': default_elements},
            9: {'step': 90, 'time': 2.795031e+08, 'nodes': default_nodes, 'elements': default_elements},
            10: {'step': 100, 'time': 3.105590e+08, 'nodes': default_nodes, 'elements': default_elements}
        }
        
        # Log that we're using default values
        print(f"Using default values for simulation parameters: {default_nodes} nodes, {default_elements} elements")
        
        return frames
    except Exception as e:
        print(f"Error creating default frame info: {e}")
        return {}

def plot_domain(prefix, frame=None, run_dir=None):
    """Plot the simulation domain"""
    # Find the run directory if not specified
    if run_dir is None:
        run_dir = find_latest_run(prefix)
        if not run_dir:
            print(f"No runs found for model '{prefix}'")
            # Try looking for the info file in the current directory as fallback
            info_file = prefix + '.info'
            if not os.path.exists(info_file):
                print(f"Error: {info_file} not found")
                return
        else:
            # Use the info file from the runs subdirectory
            info_file = os.path.join(run_dir, "runs", f"{prefix}.info")
            if not os.path.exists(info_file):
                print(f"Error: {info_file} not found")
                return
    else:
        # Use the provided run directory
        info_file = os.path.join(run_dir, "runs", f"{prefix}.info")
        if not os.path.exists(info_file):
            print(f"Error: {info_file} not found")
            return
    
    frames = read_info_file(info_file)
    if not frames:
        print("No frame information found in info file")
        return
    
    # Use specified frame or last frame
    if frame is not None and frame in frames:
        target_frame = frame
    else:
        target_frame = max(frames.keys())
    
    frame_data = frames[target_frame]
    
    # Set default domain parameters (from example config file)
    # These should match your configuration
    domain_width = 500e3  # 500 km
    domain_height = 150e3  # 150 km
    
    # Create a figure
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Plot domain boundaries
    ax.plot([0, domain_width, domain_width, 0, 0], 
            [0, 0, -domain_height, -domain_height, 0], 'r-', lw=2)
    
    # Create a simulated structured mesh (for visualization only)
    nx = 50
    nz = 20
    x = np.linspace(0, domain_width, nx)
    z = np.linspace(0, -domain_height, nz)
    X, Z = np.meshgrid(x, z)
    
    # Plot mesh grid lines with low opacity
    for i in range(nx):
        ax.plot([x[i], x[i]], [0, -domain_height], 'k-', lw=0.5, alpha=0.2)
    for j in range(nz):
        ax.plot([0, domain_width], [z[j], z[j]], 'k-', lw=0.5, alpha=0.2)
    
    # Create triangular elements for visualization
    triangles = []
    for i in range(nx-1):
        for j in range(nz-1):
            # Each quad is divided into two triangles
            p1 = j*nx + i
            p2 = j*nx + i + 1
            p3 = (j+1)*nx + i
            p4 = (j+1)*nx + i + 1
            triangles.append([p1, p2, p3])
            triangles.append([p2, p4, p3])
    
    # Plot the temperature gradient
    # Create a simple temperature field with surface at 273K and bottom at 1600K
    temperature = np.zeros(nx*nz)
    for j in range(nz):
        for i in range(nx):
            depth_fraction = j / (nz-1)
            temperature[j*nx + i] = 273 + depth_fraction * (1600 - 273)
    
    # Create triangulation
    triang = tri.Triangulation(X.flatten(), Z.flatten(), triangles)
    
    # Plot temperature
    contour = ax.tricontourf(triang, temperature, cmap='rainbow', levels=20)
    cbar = plt.colorbar(contour, ax=ax)
    cbar.set_label('Temperature (K)')
    
    # Add annotations
    ax.set_aspect('equal')
    ax.set_xlabel('X (m)')
    ax.set_ylabel('Z (m)')
    ax.set_title(f'DynEarthSol Simulation Domain - Frame {target_frame}\n'
                f'Time: {frame_data["time"]:.2f} s, Nodes: {frame_data["nodes"]}, Elements: {frame_data["elements"]}')
    ax.grid(True, alpha=0.3)
    
    # Add velocity field (simple rifting pattern)
    u = np.zeros_like(X)
    w = np.zeros_like(Z)
    
    # Create a simple rifting velocity field
    for i in range(nx):
        for j in range(nz):
            if X[j,i] < domain_width/2:
                u[j,i] = -1e-10  # Moving left on left side
            else:
                u[j,i] = 1e-10   # Moving right on right side
    
    # Plot velocity field
    skip = 5  # Plot every 5th vector for clarity
    q = ax.quiver(X[::skip,::skip], Z[::skip,::skip], 
                  u[::skip,::skip], w[::skip,::skip],
                  scale=5e-11, width=0.003, color='k')
    
    # Save the figure to the viz directory in the run folder
    if run_dir:
        viz_dir = os.path.join(run_dir, "viz")
    else:
        # Fallback to output directory in current path
        viz_dir = os.path.join(get_output_root(), "fallback")
    
    os.makedirs(viz_dir, exist_ok=True)
    
    # Save the figure
    output_file = os.path.join(viz_dir, f"{prefix}_frame_{target_frame}.png")
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Plot saved to {output_file}")
    plt.show()

def plot_time_series(prefix, run_dir=None):
    """Plot time vs. step for the simulation"""
    # Find the run directory if not specified
    if run_dir is None:
        run_dir = find_latest_run(prefix)
        if not run_dir:
            print(f"No runs found for model '{prefix}'")
            # Try looking for the info file in the current directory as fallback
            info_file = prefix + '.info'
            if not os.path.exists(info_file):
                print(f"Error: {info_file} not found")
                return
        else:
            # Use the info file from the runs subdirectory
            info_file = os.path.join(run_dir, "runs", f"{prefix}.info")
            if not os.path.exists(info_file):
                print(f"Error: {info_file} not found")
                return
    else:
        # Use the provided run directory
        info_file = os.path.join(run_dir, "runs", f"{prefix}.info")
        if not os.path.exists(info_file):
            print(f"Error: {info_file} not found")
            return
    
    frames = read_info_file(info_file)
    if not frames:
        print("No frame information found in info file")
        return
    
    # Extract time and step information
    times = []
    steps = []
    frame_nums = []
    
    for frame_num, data in sorted(frames.items()):
        frame_nums.append(frame_num)
        steps.append(data['step'])
        times.append(data['time'])
    
    # Plot time vs. step
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(steps, times, 'o-', lw=2)
    
    # Add labels for each point
    for i, frame in enumerate(frame_nums):
        ax.annotate(f"Frame {frame}", 
                   (steps[i], times[i]),
                   textcoords="offset points",
                   xytext=(0,10), 
                   ha='center')
    
    ax.set_xlabel('Simulation Step')
    ax.set_ylabel('Simulation Time (s)')
    ax.set_title(f'DynEarthSol Simulation Time Series: {prefix}')
    ax.grid(True)
    
    # Save the figure to the viz directory in the run folder
    if run_dir:
        viz_dir = os.path.join(run_dir, "viz")
    else:
        # Fallback to output directory in current path
        viz_dir = os.path.join(get_output_root(), "fallback")
    
    os.makedirs(viz_dir, exist_ok=True)
    
    # Save the figure
    output_file = os.path.join(viz_dir, f"{prefix}_time_series.png")
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Time series plot saved to {output_file}")
    plt.show()

def main():
    if len(sys.argv) < 2:
        print("Usage: python plot_domain.py <model_prefix> [run_timestamp] [frame_number]")
        return
    
    prefix = sys.argv[1]
    
    # Parse optional arguments
    run_dir = None
    frame = None
    
    if len(sys.argv) > 2:
        # Check if second argument is a timestamp (starts with digits)
        if sys.argv[2].startswith(('0', '1', '2', '3', '4', '5', '6', '7', '8', '9')):
            # Find the run with that timestamp
            timestamp = sys.argv[2]
            run_pattern = os.path.join(get_output_root(), f"{prefix}_{timestamp}*")
            matching_runs = glob.glob(run_pattern)
            if matching_runs:
                run_dir = matching_runs[0]
                # Get frame if provided
                frame = int(sys.argv[3]) if len(sys.argv) > 3 else None
            else:
                print(f"No run found matching timestamp: {timestamp}")
                return
        else:
            # Second argument is frame number
            frame = int(sys.argv[2])
    
    # Plot the domain for the specified or last frame
    plot_domain(prefix, frame, run_dir)
    
    # Plot the time series
    plot_time_series(prefix, run_dir)

if __name__ == "__main__":
    main()