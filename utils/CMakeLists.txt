# Utilities for DynEarthSol
# Python scripts don't need to be built, but we can install them

# Create a list of utility scripts
set(UTILITY_SCRIPTS
    2vtk.py
    plot_domain.py
    analyze_performance.py
)

# Copy scripts to build directory
foreach(SCRIPT ${UTILITY_SCRIPTS})
    configure_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/${SCRIPT}
        ${CMAKE_BINARY_DIR}/bin/${SCRIPT}
        COPYONLY
    )
endforeach()

# Install scripts
install(FILES ${UTILITY_SCRIPTS} DESTINATION bin
        PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE 
                    GROUP_READ GROUP_EXECUTE 
                    WORLD_READ WORLD_EXECUTE)