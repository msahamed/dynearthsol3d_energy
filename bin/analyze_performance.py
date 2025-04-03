#!/usr/bin/env python
"""
Performance analysis tool for DynEarthSol simulations.
Analyzes simulation time, mesh quality, and performance metrics.
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

def read_info_file(filename):
    """Read simulation info file and extract performance data"""
    if not os.path.exists(filename):
        print(f"Error: File '{filename}' not found")
        return None
        
    frames = []
    steps = []
    times = []
    dt_values = []
    mesh_quality = []
    nodes = []
    elements = []
    markers = []
    
    try:
        with open(filename, 'r') as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 9 and parts[0].isdigit():
                    frames.append(int(parts[1]))
                    steps.append(int(parts[2]))
                    times.append(float(parts[3]))
                    dt_values.append(float(parts[4]))
                    mesh_quality.append(float(parts[5]))
                    nodes.append(int(parts[6]))
                    elements.append(int(parts[7]))
                    markers.append(int(parts[8]) if len(parts) > 8 else 0)
    except Exception as e:
        print(f"Error reading file: {e}")
        return None
        
    return {
        'frames': np.array(frames),
        'steps': np.array(steps),
        'times': np.array(times),
        'dt': np.array(dt_values),
        'quality': np.array(mesh_quality),
        'nodes': np.array(nodes),
        'elements': np.array(elements),
        'markers': np.array(markers)
    }

def analyze_performance(prefix):
    """Analyze performance metrics from simulation data"""
    info_file = prefix + '.info'
    data = read_info_file(info_file)
    
    if data is None:
        return
        
    # Create figure with multiple subplots
    fig, axs = plt.subplots(2, 2, figsize=(12, 10))
    
    # Time step size vs. step number
    axs[0, 0].plot(data['steps'], data['dt'], 'o-', lw=2)
    axs[0, 0].set_xlabel('Step')
    axs[0, 0].set_ylabel('Time Step (s)')
    axs[0, 0].set_title('Time Step Evolution')
    axs[0, 0].grid(True)
    
    # Convert to years formatter
    year_formatter = FuncFormatter(lambda x, pos: f'{x/(365*24*3600):.1f}')
    
    # Simulation time vs. step number
    axs[0, 1].plot(data['steps'], data['times'], 'o-', lw=2)
    axs[0, 1].set_xlabel('Step')
    axs[0, 1].set_ylabel('Simulation Time (s)')
    axs[0, 1].yaxis.set_major_formatter(year_formatter)
    axs[0, 1].set_title('Simulation Progress (Years)')
    axs[0, 1].grid(True)
    
    # Mesh quality vs. step number
    axs[1, 0].plot(data['steps'], data['quality'], 'o-', lw=2)
    axs[1, 0].set_xlabel('Step')
    axs[1, 0].set_ylabel('Mesh Quality')
    axs[1, 0].set_title('Mesh Quality Evolution')
    axs[1, 0].grid(True)
    
    # Nodes and elements vs. step number
    ax1 = axs[1, 1]
    ax1.plot(data['steps'], data['nodes'], 'o-', label='Nodes', lw=2)
    ax1.set_xlabel('Step')
    ax1.set_ylabel('Number of Nodes')
    ax1.tick_params(axis='y', labelcolor='blue')
    
    # Create a second y-axis
    ax2 = ax1.twinx()
    ax2.plot(data['steps'], data['elements'], 'ro-', label='Elements', lw=2)
    ax2.set_ylabel('Number of Elements')
    ax2.tick_params(axis='y', labelcolor='red')
    
    # Add legend
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper left')
    
    ax1.set_title('Mesh Size Evolution')
    ax1.grid(True)
    
    # Adjust layout and add title
    plt.tight_layout()
    fig.suptitle(f'Performance Analysis: {prefix}', fontsize=16, y=1.02)
    
    # Create output directory if it doesn't exist
    output_dir = "output"
    os.makedirs(output_dir, exist_ok=True)
    
    # Save the figure
    output_file = os.path.join(output_dir, f"{prefix}_performance.png")
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Performance analysis saved to {output_file}")
    plt.show()

def main():
    if len(sys.argv) < 2:
        print("Usage: python analyze_performance.py <model_prefix>")
        return
        
    prefix = sys.argv[1]
    analyze_performance(prefix)

if __name__ == "__main__":
    main()