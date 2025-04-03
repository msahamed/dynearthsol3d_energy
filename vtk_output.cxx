#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <cstring>
#include <unistd.h>
#include <limits.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <ctime>
#include <chrono>

#include "constants.hpp"
#include "parameters.hpp"
#include "vtk_output.hpp"
#include "geometry.hpp"

namespace vtk_output {

void create_directory(const std::string& dirname) {
#if defined(_WIN32)
    int ret = mkdir(dirname.c_str());
#else
    int ret = mkdir(dirname.c_str(), 0755);
#endif
    if (ret != 0 && errno != EEXIST) {
        std::cerr << "Error creating directory: " << dirname << " - " << strerror(errno) << "\n";
    }
    else if (ret == 0) {
        std::cout << "  Created directory: " << dirname << std::endl;
    }
}

// Helper function to get timestamp string
std::string get_timestamp() {
    auto now = std::chrono::system_clock::now();
    std::time_t now_time = std::chrono::system_clock::to_time_t(now);
    
    std::stringstream ss;
    ss << std::put_time(std::localtime(&now_time), "%Y%m%d_%H%M%S");
    return ss.str();
}

// Global timestamp for this run (set once at first call)
static std::string run_timestamp = "";

// Helper function to get the output root path
std::string get_output_root() {
    // Use the PROJECT_ROOT_DIR definition from CMake if available
#ifdef PROJECT_ROOT_DIR
    return std::string(PROJECT_ROOT_DIR) + "/output";
#else
    // Fallback to using getcwd and manual detection
    char cwd[PATH_MAX];
    if (getcwd(cwd, sizeof(cwd)) == NULL) {
        std::cerr << "Error getting current working directory" << std::endl;
        return "./output"; // Fallback to relative path
    }
    
    std::string cwd_str(cwd);
    size_t build_pos = cwd_str.find("/build");
    if (build_pos != std::string::npos) {
        // We're in build directory, go up one level
        return cwd_str.substr(0, build_pos) + "/output";
    } else {
        // We're in project root
        return cwd_str + "/output";
    }
#endif
}

// Get the run-specific directory for a given model
std::string get_run_directory(const std::string& modelname) {
    // Initialize timestamp if it's the first call
    if (run_timestamp.empty()) {
        run_timestamp = get_timestamp();
    }
    
    // Create a directory name with model and timestamp
    return get_output_root() + "/" + modelname + "_" + run_timestamp;
}

void setup_output_directories(const std::string& modelname) {
    // Get the output root directory
    std::string output_root = get_output_root();
    
    // Create main output directory
    create_directory(output_root);
    
    // Get the run-specific directory
    std::string run_dir = get_run_directory(modelname);
    create_directory(run_dir);
    
    // Create subdirectories within the run directory
    std::string runs_subdir = run_dir + "/runs";
    std::string vtk_subdir = run_dir + "/vtk";
    std::string viz_subdir = run_dir + "/viz";
    
    create_directory(runs_subdir);
    create_directory(vtk_subdir);
    create_directory(viz_subdir);
    
    std::cout << "  Created output directories for model: " << modelname << " (timestamp: " << run_timestamp << ")" << std::endl;
}

void write_vtk_file(const Variables& var, int frame, double dt, const std::string& modelname) {
    // Ensure output directories exist
    setup_output_directories(modelname);
    
    // Construct filename with run-specific path
    std::string run_dir = get_run_directory(modelname);
    std::string vtk_dir = run_dir + "/vtk";
    std::stringstream filename;
    filename << vtk_dir << "/" << modelname << "_" << std::setw(6) << std::setfill('0') << frame << ".vtk";
    
    // Open output file
    std::ofstream vtk_file(filename.str().c_str());
    if (!vtk_file) {
        std::cerr << "Error: cannot open VTK file " << filename.str() << " for writing.\n";
        return;
    }
    
    // Write VTK file header
    vtk_file << "# vtk DataFile Version 3.0\n";
    vtk_file << "DynEarthSol output, frame " << frame << ", time " << var.time << " years\n";
    vtk_file << "ASCII\n";
    vtk_file << "DATASET UNSTRUCTURED_GRID\n";
    
    // Write node coordinates
    vtk_file << "POINTS " << var.nnode << " double\n";
    for (int i = 0; i < var.nnode; i++) {
        for (int j = 0; j < NDIMS; j++) {
            vtk_file << (*var.coord)[i][j] << " ";
        }
        if (NDIMS == 2) {
            vtk_file << "0.0";  // For 2D, add z=0 coordinate
        }
        vtk_file << "\n";
    }
    
    // Write element connectivity
    int points_per_elem = NODES_PER_ELEM;
    vtk_file << "CELLS " << var.nelem << " " << var.nelem * (points_per_elem + 1) << "\n";
    for (int i = 0; i < var.nelem; i++) {
        vtk_file << points_per_elem << " ";
        for (int j = 0; j < points_per_elem; j++) {
            vtk_file << (*var.connectivity)[i][j] << " ";
        }
        vtk_file << "\n";
    }
    
    // Write cell types
    vtk_file << "CELL_TYPES " << var.nelem << "\n";
    int vtk_cell_type;
    if (NDIMS == 2) {
        vtk_cell_type = 5;  // VTK_TRIANGLE
    } else {
        vtk_cell_type = 10; // VTK_TETRA
    }
    for (int i = 0; i < var.nelem; i++) {
        vtk_file << vtk_cell_type << "\n";
    }
    
    // Write point data (nodal fields)
    vtk_file << "POINT_DATA " << var.nnode << "\n";
    
    // Temperature
    vtk_file << "SCALARS temperature double 1\n";
    vtk_file << "LOOKUP_TABLE default\n";
    for (int i = 0; i < var.nnode; i++) {
        vtk_file << (*var.temperature)[i] << "\n";
    }
    
    // Velocity
    vtk_file << "VECTORS velocity double\n";
    for (int i = 0; i < var.nnode; i++) {
        for (int j = 0; j < NDIMS; j++) {
            vtk_file << (*var.vel)[i][j] << " ";
        }
        if (NDIMS == 2) {
            vtk_file << "0.0";  // For 2D, add vz=0 component
        }
        vtk_file << "\n";
    }
    
    // Element data
    vtk_file << "CELL_DATA " << var.nelem << "\n";
    
    // Strain rate
    vtk_file << "SCALARS strain_rate double 1\n";
    vtk_file << "LOOKUP_TABLE default\n";
    for (int i = 0; i < var.nelem; i++) {
        // Compute magnitude of strain rate tensor
        double sr_magnitude = 0;
        for (int j = 0; j < NSTR; j++) {
            sr_magnitude += (*var.strain_rate)[i][j] * (*var.strain_rate)[i][j];
        }
        sr_magnitude = std::sqrt(sr_magnitude);
        vtk_file << sr_magnitude << "\n";
    }
    
    // Material type
    vtk_file << "SCALARS material_type int 1\n";
    vtk_file << "LOOKUP_TABLE default\n";
    for (int i = 0; i < var.nelem; i++) {
        // Find the most abundant marker mattype in this element
        int_vec &a = (*var.elemmarkers)[i];
        int mattype = std::distance(a.begin(), std::max_element(a.begin(), a.end()));
        vtk_file << mattype << "\n";
    }
    
    // Stress (for 2D, output as 3-component vector)
    if (NDIMS == 2) {
        vtk_file << "VECTORS stress double\n";
        for (int i = 0; i < var.nelem; i++) {
            vtk_file << (*var.stress)[i][0] << " " // xx
                     << (*var.stress)[i][1] << " " // zz
                     << (*var.stress)[i][2] << "\n"; // xz
        }
    } else {
        // For 3D, output stress as tensor
        vtk_file << "TENSORS stress double\n";
        for (int i = 0; i < var.nelem; i++) {
            // Row 1: xx, xy, xz
            vtk_file << (*var.stress)[i][0] << " " 
                     << (*var.stress)[i][3] << " " 
                     << (*var.stress)[i][5] << "\n";
            // Row 2: xy, yy, yz
            vtk_file << (*var.stress)[i][3] << " " 
                     << (*var.stress)[i][1] << " " 
                     << (*var.stress)[i][4] << "\n";
            // Row 3: xz, yz, zz
            vtk_file << (*var.stress)[i][5] << " " 
                     << (*var.stress)[i][4] << " " 
                     << (*var.stress)[i][2] << "\n";
        }
    }
    
    // Plastic strain
    vtk_file << "SCALARS plastic_strain double 1\n";
    vtk_file << "LOOKUP_TABLE default\n";
    for (int i = 0; i < var.nelem; i++) {
        vtk_file << (*var.plstrain)[i] << "\n";
    }
    
    vtk_file.close();
    
    std::cout << "  Wrote VTK file: " << filename.str() << "\n";
}

} // namespace vtk_output