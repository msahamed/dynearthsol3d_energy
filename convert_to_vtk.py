#!/usr/bin/env python3
"""
Convert DynEarthSol binary output to VTK format
"""

import sys
import os
import struct
import numpy as np
import argparse


def read_info_file(prefix):
    """Read the info file to get frame information"""
    info_file = prefix + '.info'
    if not os.path.exists(info_file):
        print(f"Error: info file '{info_file}' not found!")
        return {}
    
    frame_info = {}
    with open(info_file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 8:
                try:
                    frame_idx = int(parts[0])
                    frame_info[frame_idx] = {
                        'step': int(parts[1]),
                        'time': float(parts[2]),
                        'dt': float(parts[3]),
                        'quality': float(parts[4]),
                        'nodes': int(parts[5]),
                        'elements': int(parts[6]),
                        'markers': int(parts[7])
                    }
                except (ValueError, IndexError):
                    continue
    
    return frame_info


def read_binary_header(filename):
    """Read the header from a binary save file"""
    sections = {}
    ndims = 2
    
    with open(filename, 'rb') as f:
        # Read first line to get ndims
        first_line = f.readline().decode('utf-8')
        if 'ndims=' in first_line:
            ndims = int(first_line.split('ndims=')[1].split()[0])
        
        # Read section offsets
        while True:
            line = f.readline().decode('utf-8', errors='ignore')
            if not line or line == '\n' or '\x00' in line:
                break
            
            parts = line.strip().split('\t')
            if len(parts) == 2:
                name = parts[0]
                try:
                    offset = int(parts[1])
                    sections[name] = offset
                except ValueError:
                    pass
    
    return sections, ndims


def read_array_2d(f, count, cols):
    """Read a 2D array from binary file"""
    data = np.zeros((count, cols))
    for i in range(count):
        for j in range(cols):
            data[i, j] = struct.unpack('d', f.read(8))[0]
    return data


def read_array_1d(f, count):
    """Read a 1D array from binary file"""
    data = np.zeros(count)
    for i in range(count):
        data[i] = struct.unpack('d', f.read(8))[0]
    return data


def read_int_array_2d(f, count, cols):
    """Read a 2D integer array from binary file"""
    data = np.zeros((count, cols), dtype=int)
    for i in range(count):
        for j in range(cols):
            data[i, j] = struct.unpack('i', f.read(4))[0]
    return data


def convert_frame(prefix, frame_idx, frame_info):
    """Convert a single frame to VTK"""
    savefile = f"{prefix}.save.{frame_idx:06d}"
    
    if not os.path.exists(savefile):
        print(f"Error: {savefile} not found")
        return False
    
    # Get frame information
    if frame_idx not in frame_info:
        print(f"Error: Frame {frame_idx} not in info file")
        return False
    
    info = frame_info[frame_idx]
    nnode = info['nodes']
    nelem = info['elements']
    
    print(f"Converting frame {frame_idx}: {nnode} nodes, {nelem} elements")
    
    # Read binary header
    sections, ndims = read_binary_header(savefile)
    print(f"  Dimensions: {ndims}D")
    print(f"  Found {len(sections)} data sections")
    
    # Read data from binary file
    coords = None
    conn = None
    velocity = None
    temperature = None
    stress = None
    strain = None
    plstrain = None
    
    with open(savefile, 'rb') as f:
        # Read coordinates
        if 'coordinate' in sections:
            f.seek(sections['coordinate'])
            coords = read_array_2d(f, nnode, ndims)
            print(f"  Read coordinates: {coords.shape}")
        
        # Read connectivity
        if 'connectivity' in sections:
            f.seek(sections['connectivity'])
            conn = read_int_array_2d(f, nelem, ndims + 1)
            print(f"  Read connectivity: {conn.shape}")
        
        # Read velocity
        if 'velocity' in sections:
            f.seek(sections['velocity'])
            velocity = read_array_2d(f, nnode, ndims)
            print(f"  Read velocity: {velocity.shape}")
        
        # Read temperature
        if 'temperature' in sections:
            f.seek(sections['temperature'])
            temperature = read_array_1d(f, nnode)
            print(f"  Read temperature: {temperature.shape}")
        
        # Read stress
        if 'stress' in sections:
            f.seek(sections['stress'])
            nstr = 3 if ndims == 2 else 6
            stress = read_array_2d(f, nelem, nstr)
            print(f"  Read stress: {stress.shape}")
        
        # Read strain
        if 'strain' in sections:
            f.seek(sections['strain'])
            nstr = 3 if ndims == 2 else 6
            strain = read_array_2d(f, nelem, nstr)
            print(f"  Read strain: {strain.shape}")
        
        # Read plastic strain
        if 'plastic strain' in sections:
            f.seek(sections['plastic strain'])
            plstrain = read_array_1d(f, nelem)
            print(f"  Read plastic strain: {plstrain.shape}")
    
    # Write VTK file
    vtkfile = f"{prefix}.{frame_idx:06d}.vtk"
    
    with open(vtkfile, 'w') as f:
        # Header
        f.write('# vtk DataFile Version 3.0\n')
        f.write(f'DynEarthSol output, frame {frame_idx}, step {info["step"]}, time {info["time"]:.3e} s\n')
        f.write('ASCII\n')
        f.write('DATASET UNSTRUCTURED_GRID\n')
        
        # Points
        f.write(f'POINTS {nnode} double\n')
        for i in range(nnode):
            if ndims == 2:
                f.write(f'{coords[i,0]:.6e} {coords[i,1]:.6e} 0.0\n')
            else:
                f.write(f'{coords[i,0]:.6e} {coords[i,1]:.6e} {coords[i,2]:.6e}\n')
        
        # Cells
        points_per_cell = ndims + 1
        f.write(f'CELLS {nelem} {nelem * (points_per_cell + 1)}\n')
        for i in range(nelem):
            if ndims == 2:
                f.write(f'3 {conn[i,0]} {conn[i,1]} {conn[i,2]}\n')
            else:
                f.write(f'4 {conn[i,0]} {conn[i,1]} {conn[i,2]} {conn[i,3]}\n')
        
        # Cell types
        cell_type = 5 if ndims == 2 else 10  # Triangle or Tetrahedron
        f.write(f'CELL_TYPES {nelem}\n')
        for i in range(nelem):
            f.write(f'{cell_type}\n')
        
        # Point data
        if temperature is not None or velocity is not None:
            f.write(f'POINT_DATA {nnode}\n')
            
            if temperature is not None:
                f.write('SCALARS temperature double 1\n')
                f.write('LOOKUP_TABLE default\n')
                for i in range(nnode):
                    f.write(f'{temperature[i]:.6e}\n')
            
            if velocity is not None:
                f.write('VECTORS velocity double\n')
                for i in range(nnode):
                    if ndims == 2:
                        f.write(f'{velocity[i,0]:.6e} {velocity[i,1]:.6e} 0.0\n')
                    else:
                        f.write(f'{velocity[i,0]:.6e} {velocity[i,1]:.6e} {velocity[i,2]:.6e}\n')
        
        # Cell data
        if stress is not None or strain is not None or plstrain is not None:
            f.write(f'CELL_DATA {nelem}\n')
            
            if plstrain is not None:
                f.write('SCALARS plastic_strain double 1\n')
                f.write('LOOKUP_TABLE default\n')
                for i in range(nelem):
                    f.write(f'{plstrain[i]:.6e}\n')
            
            if stress is not None:
                # Write von Mises stress
                f.write('SCALARS von_mises_stress double 1\n')
                f.write('LOOKUP_TABLE default\n')
                for i in range(nelem):
                    if ndims == 2:
                        # 2D: sxx, szz, sxz
                        sxx, szz, sxz = stress[i]
                        vm = np.sqrt(sxx**2 + szz**2 - sxx*szz + 3*sxz**2)
                    else:
                        # 3D: sxx, syy, szz, sxy, sxz, syz
                        sxx, syy, szz, sxy, sxz, syz = stress[i]
                        vm = np.sqrt(0.5*((sxx-syy)**2 + (syy-szz)**2 + (szz-sxx)**2 + 6*(sxy**2 + sxz**2 + syz**2)))
                    f.write(f'{vm:.6e}\n')
    
    print(f"Created {vtkfile}")
    return True


def main():
    parser = argparse.ArgumentParser(description='Convert DynEarthSol output to VTK')
    parser.add_argument('prefix', help='Prefix of input files (e.g., output/run/model)')
    parser.add_argument('-s', '--start', type=int, default=0, help='Starting frame index')
    parser.add_argument('-e', '--end', type=int, default=-1, help='Ending frame index (-1 for all)')
    args = parser.parse_args()
    
    # Read info file
    frame_info = read_info_file(args.prefix)
    if not frame_info:
        print("Error: Could not read info file or no frames found")
        return 1
    
    # Determine frames to convert
    if args.end == -1:
        frames = sorted(frame_info.keys())
    else:
        frames = [f for f in range(args.start, args.end + 1) if f in frame_info]
    
    if not frames:
        print(f"No frames found in range {args.start} to {args.end}")
        return 1
    
    print(f"Converting {len(frames)} frames: {frames[0]} to {frames[-1]}")
    
    # Convert each frame
    for frame_idx in frames:
        convert_frame(args.prefix, frame_idx, frame_info)
    
    print(f"\nConversion complete! Created {len(frames)} VTK files")
    return 0


if __name__ == '__main__':
    sys.exit(main())
