# Tests for DynEarthSol
add_executable(run_tests tests.cxx)

# Link with the same libraries as the main executable
target_link_libraries(run_tests
    3x3
    ANN
    ${Boost_LIBRARIES}
)

# Add a test
add_test(NAME UnitTests COMMAND run_tests)

# Add a simple mesh test
add_test(
    NAME MeshTest 
    COMMAND $<TARGET_FILE:dynearthsol2d> ${CMAKE_SOURCE_DIR}/examples/rifting-2d.poly --sim.max_steps=5
)