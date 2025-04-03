# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build Commands
- Build: `make` (production) or `make opt=0 openmp=0` (debugging)
- Clean: `make clean` or `make deepclean`
- Configuration options in Makefile: `ndims=2|3`, `opt=0-3`, `openmp=0|1`, `useadapt=0|1`
- Tests: Compile with `g++ --std=c++11 tests.cxx 3x3-C/lib3x3.a`

## Code Style Guidelines
- Naming: Functions/variables use `snake_case`, classes use `PascalCase_with_underscores`
- Constants: Use `UPPER_CASE` for constants (e.g., `NDIMS`, `NSTR`)
- Error codes: 1 (user input), 2 (IO), 10 (mesh), 11 (runtime), 12 (assertion)
- Avoid: C++ streams for bulk output, creating/destroying objects in inner loops, static/global variables
- Parallelization: Code uses OpenMP for parallelization when enabled

## Project Structure
- Core C++ simulation engine with C++11 support
- Dependencies: Boost::Program_options, Python with NumPy
- Optional: VTK and MPI (for libadaptivity)
- Main workflows: Define model in input file, run simulation, convert output with 2vtk.py, visualize with Paraview/Visit