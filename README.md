## DynEarthSol3D_energy : A numerical tool for tectonic modeling with thermomechanics

### Basic Information:
DynEarthSol3D is a finite element code that solves the momentum balance and 
the heat transfer in Lagrangian form using unstructured meshes. It can be
used to study the long-term deformation of Earth's lithosphere and problems
alike.

The purpose of this fork of DES3D (https://bitbucket.org/tan2/dynearthsol3d) 
is to show that libadaptivity from Fluidity
(https://github.com/FluidityProject/fluidity), 
a self-contained library for anisotropic adaptive mesh refinement, 
works well as a mesh optimizer for DES3D's remeshing.

### Thermomechanics module:
In this code i have added the full energy balance equation with mass 
conservation. For anything about these implementations please feel free
to contact :

Sabber Ahamed
sabbers@gmail.com 
Senior Data Scientist
Bridgestone America
Nashville, TN, USA

### Installation

<li> You will need a C++ compiler that supports C++17 standard or higher. (GNU g++
  7.0 or newer version will suffice.)
<li> You will need CMake 3.14 or newer.
<li> You will need a recent version of Boost::Program_options library (1.42 or
  newer version).
<li> You will need Python 3.6+ with Numpy and Matplotlib packages.

For libadaptiviy, you further need
<li> VTK (v.5.10 tested) built from source or development packages.
<li> MPI (openmpi-1.6.1 tested).

Build procedure:
- Using the build script:
<li> Run the included build script: `./build.sh` 
<li> This will configure and build the code with default options.
<li> For different configurations, run `./build.sh --help` to see available options.

- Manual CMake build:
<li> Create a build directory: `mkdir build && cd build`
<li> Configure with CMake: `cmake .. -DCMAKE_BUILD_TYPE=Release -DWITH_3D=OFF`
<li> Build the code: `make -j4`

- Build options:
<li> `-DWITH_3D=ON/OFF` to build 2D or 3D version
<li> `-DWITH_OPENMP=ON/OFF` to enable/disable OpenMP
<li> `-DWITH_ADAPTIVITY=ON/OFF` to enable/disable mesh adaptivity

### How to run

<li> Execute "dynearthsol2d inputfile".
<li> Several example input files are provided under 'examples/' directory. The
  format of the input file is described in 'examples/defaults.cfg'.
<li> Benchmark cases with analytical solution can be found under 'benchmarks/'
  directory.
<li> Execute the executable with '-h' flag to see the available input parameters
  and their descriptions.

### How to plot

<li> The simulation now outputs VTK files directly in a well-organized output directory structure.
<li> Output is organized into three main directories:
  <ul>
  <li> `output/runs/<modelname>/` - Contains binary simulation data (.save, .info, .chkpt files)
  <li> `output/vtk/<modelname>/` - Contains VTK visualization files
  <li> `output/viz/<modelname>/` - Contains visualization outputs (images, animations)
  </ul>
<li> Run "python plot_domain.py modelname" to create visualization plots of the domain.
<li> Run "python utils/animate_simulation.py modelname" to create animations of the temperature and stress fields over time.
<li> Plot the VTK files with Paraview or LLNL's Visit program.
<li> See the README.md in the output directory for more details on the output organization.


### Bug reports
      
Bug reports, comments, and suggestions are always welcome. The best 
channel is to create an issue on the Issue Tracker here at the original code:
   <a>http://bitbucket.org/tan2/dynearthsol3d</a>
   or this version with thermomechanics :
   <a>https://github.com/msahamed/dynearthsol3d_energy</a>


### License

This program is free software: you can redistribute it and/or modify
it under the terms of the MIT / X Windows System license (see the
file LICENSE for the full text).

The files under the subdirectories 3x3-C/, ann/, tetgen/, and
triangles/ are distributed by their own license(s).
