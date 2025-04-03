# CMakeLists.txt for ANN library
add_library(ANN STATIC
    src/ANN.cpp
    src/brute.cpp
    src/kd_tree.cpp
    src/kd_util.cpp
    src/kd_split.cpp
    src/kd_dump.cpp
    src/kd_search.cpp
    src/kd_pr_search.cpp
    src/kd_fix_rad_search.cpp
    src/bd_tree.cpp
    src/bd_search.cpp
    src/bd_pr_search.cpp
    src/bd_fix_rad_search.cpp
    src/perf.cpp
)

# Set compile flags
# Use C++11 for this library since it uses the 'register' keyword which is deprecated in C++17
set_target_properties(ANN PROPERTIES 
    COMPILE_FLAGS "-O3"
    CXX_STANDARD 11
    CXX_STANDARD_REQUIRED YES
)

# Include directories
target_include_directories(ANN PUBLIC 
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)