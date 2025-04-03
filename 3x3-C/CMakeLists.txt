# CMakeLists.txt for 3x3-C library
add_library(3x3 STATIC
    dsytrd3.c
    dsyevc3.c
    dsyevq3.c
    dsyevh3.c
)

# Set compile flags
set_source_files_properties(
    dsytrd3.c
    dsyevc3.c
    dsyevq3.c
    dsyevh3.c
    PROPERTIES COMPILE_FLAGS "-Wall -O3 -ffast-math"
)

# Include directories
target_include_directories(3x3 PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})