#ifndef DYNEARTHSOL3D_VTK_OUTPUT_HPP
#define DYNEARTHSOL3D_VTK_OUTPUT_HPP

#include <string>
#include "parameters.hpp"

class Variables;

namespace vtk_output {

    /**
     * Write simulation state to VTK file
     *
     * @param var Variables containing simulation state
     * @param frame Frame number for the output
     * @param dt Current time step
     * @param modelname Base name for output files
     */
    void write_vtk_file(const Variables& var, int frame, double dt, const std::string& modelname);

    /**
     * Create output directory structure
     *
     * @param modelname Base name for output files
     */
    void setup_output_directories(const std::string& modelname);
    
    /**
     * Get the absolute path to the output root directory
     * 
     * @return Absolute path to the output directory
     */
    std::string get_output_root();
    
    /**
     * Get the run-specific directory for a model
     * 
     * @param modelname Base name for the model
     * @return Path to the run directory (includes timestamp)
     */
    std::string get_run_directory(const std::string& modelname);
    
    /**
     * Get a timestamp string for the current time
     * 
     * @return Formatted timestamp string (YYYYMMDD_HHMMSS)
     */
    std::string get_timestamp();

} // namespace vtk_output

#endif // DYNEARTHSOL3D_VTK_OUTPUT_HPP