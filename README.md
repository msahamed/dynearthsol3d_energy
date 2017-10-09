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
msahamed@memphis.edu
Center For Earthquake Research and Information (CERI)
The University of Memphis
3890 Central Ave
Memphis, TN 38152, USA

### Installation

<li> You will need a recent C++ compiler that supports C++11 standard. (GNU g++
  4.4 or newer version will suffice.)
<li> You will need a recent version of Boost::Program_options library (1.42 or
  newer version). Instructions for building the library:
  -- Download the source code from www.boost.org
  -- In the untarred source directory, run "./bootstrap.sh"
  -- In the same directory, run "./b2 --with-program_options -q" to build
     the library.
<li> You will need Python 2.6+ or 3.2+ and the Numpy package.

For libadaptiviy, you further need
<li> VTK (v.5.10 tested) built from source or development packages.
<li> MPI (openmpi-1.6.1 tested).

Build procedure:
- libadaptivity
<li> Run "LDFLAGS=-L${VTK_LIBDIR} ./configure" in libadaptivity
<li> Run "make".
<li> Run "make" in tests/. Depending on the flavor of MPI, one might need to add '-lmpi_f77' to LIBS in Makefile.

- DES3D
<li> Edit 'Makefile', 
  1) modify BOOST_ROOT_DIR if you manually built or installed 
  boost library. If you followed the instructions above to build 
  Boost::Program_options library, set BOOST_ROOT_DIR to the untarred boost
  directory.
  2) turn on 'useadapt' to use libadaptivity.
  3) Set 'ndims = 3'. libadaptivity works only for 3D.
  4) Check if VTK_INC path is correct.
  5) Copy LIBS in libadaptivity/tests/Makefile to LIBADAPTIVITY_LIBS.
  6) Make sure that the path to libadaptivity.a is correctly set in the beginning of LIBADAPTIVITY_LIBS. 
     $(LIBADAPTIVITY_LIB)/libadaptivity.a will do.

<li> Run "make opt=0" to build a debugging executable.
<li> Run "make openmp=0" to build the executable without OpenMP. This is
  necessary to debug the code under valgrind.

### How to run

<li> Execute "dynearthsol2d inputfile".
<li> Several example input files are provided under 'examples/' directory. The
  format of the input file is described in 'examples/defaults.cfg'.
<li> Benchmark cases with analytical solution can be found under 'benchmarks/'
  directory.
<li> Execute the executable with '-h' flag to see the available input parameters
  and their descriptions.

### How to plot

<li> Run "2vtk.py modelname" to convert the binary output to VTK files.
<li> Execute 2vtk.py with '-h' flag to see more usage information.
<li> Some of the simulation outputs might be disabled. Edit 2vtk.py and
  output.cxx to disable/enable them.
<li> Plot the VTK files with Paraview or LLNL's Visit program.


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
