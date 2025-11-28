cmake_minimum_required(VERSION 3.14)
project(DynEarthSol VERSION 1.0.0 LANGUAGES CXX C)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Set output directories to avoid duplication
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# Define project root directory for absolute paths
set(PROJECT_ROOT_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
add_definitions(-DPROJECT_ROOT_DIR="${PROJECT_ROOT_DIR}")

# Add compile options
option(WITH_OPENMP "Build with OpenMP support" ON)
option(WITH_OPENCL "Build with OpenCL support for GPU acceleration" OFF)
option(WITH_3D "Build 3D version" OFF)
option(WITH_ADAPT "Build with libadaptivity for mesh optimization" OFF)
option(WITH_DEBUG "Build with debug symbols" OFF)

# Set compiler flags
if(WITH_DEBUG)
    set(CMAKE_BUILD_TYPE Debug)
    add_compile_options(-O0 -Wall -Wno-unused-variable -Wno-unused-function -Wno-unknown-pragmas)
else()
    set(CMAKE_BUILD_TYPE Release)
    add_compile_options(-O2)
endif()

# Add dimension definition
if(WITH_3D)
    add_definitions(-DTHREED)
    set(DIM_SUFFIX "3d")
else()
    set(DIM_SUFFIX "2d")
endif()

# Set up OpenMP if requested
# Set up OpenMP if requested
if(WITH_OPENMP)
    # Hint for macOS Homebrew OpenMP
    if(APPLE)
        if(EXISTS "/opt/homebrew/opt/libomp")
            set(OpenMP_ROOT "/opt/homebrew/opt/libomp")
        elseif(EXISTS "/usr/local/opt/libomp")
            set(OpenMP_ROOT "/usr/local/opt/libomp")
        endif()
    endif()

    find_package(OpenMP)
    if(OpenMP_CXX_FOUND)
        message(STATUS "OpenMP Found!")
        
        # Override flags for Apple Clang if needed
        if(CMAKE_CXX_COMPILER_ID MATCHES "AppleClang")
            message(STATUS "Overriding OpenMP flags for AppleClang")
            set(OpenMP_CXX_FLAGS "-Xpreprocessor" "-fopenmp")
            set(OpenMP_C_FLAGS "-Xpreprocessor" "-fopenmp")
        endif()
        
        message(STATUS "  OpenMP_CXX_FLAGS: ${OpenMP_CXX_FLAGS}")
        message(STATUS "  OpenMP_CXX_LIBRARIES: ${OpenMP_CXX_LIBRARIES}")
        
        add_compile_options(${OpenMP_CXX_FLAGS})
        add_definitions(-DUSE_OMP)
        # We will link in target_link_libraries later
    else()
        message(WARNING "OpenMP requested but not found. Please install libomp (e.g. brew install libomp)")
    endif()
endif()

# Set up OpenCL if requested
if(WITH_OPENCL)
    find_package(OpenCL)
    if(OpenCL_FOUND)
        include_directories(${OpenCL_INCLUDE_DIRS})
        add_definitions(-DWITH_OPENCL)
        message(STATUS "Found OpenCL: ${OpenCL_LIBRARIES}")
    else()
        message(WARNING "OpenCL requested but not found")
    endif()
endif()

# Find Boost
find_package(Boost REQUIRED COMPONENTS program_options)
include_directories(${Boost_INCLUDE_DIRS})

# Find Eigen3 for vectorized math
find_package(Eigen3 3.3 QUIET)

# Fallback for Homebrew on Apple Silicon
if(NOT Eigen3_FOUND)
    if(EXISTS "/opt/homebrew/include/eigen3")
        set(EIGEN3_INCLUDE_DIR "/opt/homebrew/include/eigen3")
        set(Eigen3_FOUND TRUE)
    endif()
endif()

if(Eigen3_FOUND)
    message(STATUS "Found Eigen3: ${EIGEN3_INCLUDE_DIR}")
    add_definitions(-DUSE_EIGEN)
    include_directories(${EIGEN3_INCLUDE_DIR})
else()
    message(STATUS "Eigen3 not found. Please install it (e.g. 'brew install eigen') for better performance.")
endif()

# Add include directories
include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/ann/include
)

# Source files (common)
set(SOURCES
    barycentric-fn.cxx
    benchmark.cxx
    brc-interpolation.cxx
    bc.cxx
    binaryio.cxx
    dynearthsol.cxx
    fields.cxx
    geometry.cxx
    ic.cxx
    ic-read-temp.cxx
    input.cxx
    matprops.cxx
    mesh.cxx
    nn-interpolation.cxx
    output.cxx
    phasechanges.cxx
    rheology.cxx
    markerset.cxx
    vtk_output.cxx
)

# Add OpenCL sources if enabled
if(WITH_OPENCL)
    list(APPEND SOURCES rheology-opencl.cpp)
endif()

# Add remeshing file based on adaptivity option
if(WITH_ADAPT)
    list(APPEND SOURCES remeshing_adapt.cxx)
    
    # Add required libraries for adaptivity
    find_package(VTK REQUIRED)
    include_directories(
        ${CMAKE_CURRENT_SOURCE_DIR}/libadaptivity/include
        ${CMAKE_CURRENT_SOURCE_DIR}/libadaptivity/adapt3d/include
        ${CMAKE_CURRENT_SOURCE_DIR}/libadaptivity/metric_field/include
        ${CMAKE_CURRENT_SOURCE_DIR}/libadaptivity/load_balance/include
        ${VTK_INCLUDE_DIRS}
    )
    add_definitions(-DADAPT -DHAVE_VTK=1)
else()
    list(APPEND SOURCES remeshing.cxx)
endif()

# Create executable
add_executable(dynearthsol${DIM_SUFFIX} ${SOURCES})

# Libraries
# 3x3 matrix library
add_subdirectory(3x3-C)

# ANN library
add_subdirectory(ann)

# Triangle/Tetgen libraries
if(WITH_3D)
    # For 3D, add tetgen
    add_library(tetgen STATIC 
        tetgen/predicates.cxx
        tetgen/tetgen.cxx
    )
    set_source_files_properties(tetgen/predicates.cxx PROPERTIES COMPILE_FLAGS "-O0 -DTETLIBRARY")
    set_source_files_properties(tetgen/tetgen.cxx PROPERTIES COMPILE_FLAGS "-DNDEBUG -DTETLIBRARY")
    
    target_link_libraries(dynearthsol${DIM_SUFFIX} tetgen)
else()
    # For 2D, add triangle
    add_library(triangle STATIC triangle/triangle.c)
    set_source_files_properties(triangle/triangle.c PROPERTIES 
        COMPILE_FLAGS "-O1 -DTRILIBRARY -DREDUCED -DANSI_DECLARATORS"
        LANGUAGE CXX  # Force C++ compiler for triangle
    )
    
    target_link_libraries(dynearthsol${DIM_SUFFIX} triangle)
endif()

# Link libraries
target_link_libraries(dynearthsol${DIM_SUFFIX} 
    3x3
    ANN
    ${Boost_LIBRARIES}
)

if(OpenMP_CXX_FOUND)
    target_link_libraries(dynearthsol${DIM_SUFFIX} OpenMP::OpenMP_CXX)
endif()

# Add OpenCL if enabled
if(WITH_OPENCL)
    target_link_libraries(dynearthsol${DIM_SUFFIX} ${OpenCL_LIBRARIES})
    
    # Copy the OpenCL kernel file to the build directory
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/stress_update_kernel.cl 
                   ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/stress_update_kernel.cl COPYONLY)
endif()

# Add libadaptivity if enabled
if(WITH_ADAPT)
    add_subdirectory(libadaptivity)
    target_link_libraries(dynearthsol${DIM_SUFFIX} adaptivity ${VTK_LIBRARIES})
endif()

# Installation
install(TARGETS dynearthsol${DIM_SUFFIX} DESTINATION bin)
if(WITH_OPENCL)
    install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/stress_update_kernel.cl DESTINATION bin)
endif()

# Add Python utilities directory if it exists
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/utils/CMakeLists.txt")
    add_subdirectory(utils)
endif()

# Disable testing for now
# if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/tests/CMakeLists.txt")
#     enable_testing()
#     add_subdirectory(tests)
# endif()

# Output configuration info
message(STATUS "Configuration summary:")
message(STATUS "  Build type: ${CMAKE_BUILD_TYPE}")
if(WITH_3D)
    message(STATUS "  Dimensions: 3D")
else()
    message(STATUS "  Dimensions: 2D")
endif()
message(STATUS "  OpenMP: ${WITH_OPENMP}")
message(STATUS "  OpenCL: ${WITH_OPENCL}")
message(STATUS "  Adaptivity: ${WITH_ADAPT}")