#include <iostream>

#ifdef USE_OMP
#include <omp.h>
#endif

#include "constants.hpp"
#include "parameters.hpp"
#include "benchmark.hpp"
#include "bc.hpp"
#include "binaryio.hpp"
#include "fields.hpp"
#include "geometry.hpp"
#include "ic.hpp"
#include "input.hpp"
#include "matprops.hpp"
#include "markerset.hpp"
#include "mesh.hpp"
#include "output.hpp"
#include "phasechanges.hpp"
#include "remeshing.hpp"
#include "rheology.hpp"
#ifdef WITH_OPENCL
#include "rheology-opencl.hpp"
#endif

#ifdef WIN32
#ifdef _MSC_VER
#define snprintf _snprintf
#endif // _MSC_VER
namespace std { using ::snprintf; }
#endif // WIN32

void init_var(const Param& param, Variables& var)
{
    var.time = 0;
    var.steps = 0;

    if (param.control.characteristic_speed == 0)
        var.max_vbc_val = find_max_vbc(param.bc);
    else
        var.max_vbc_val = param.control.characteristic_speed;

    // XXX: Hard coded boundary flag. If the order of ibound?? is changed
    //      in the future, the following lines have to be updated as well.
    var.vbc_types[0] = param.bc.vbc_x0;
    var.vbc_types[1] = param.bc.vbc_x1;
    var.vbc_types[2] = param.bc.vbc_y0;
    var.vbc_types[3] = param.bc.vbc_y1;
    var.vbc_types[4] = param.bc.vbc_z0;
    var.vbc_types[5] = param.bc.vbc_z1;
    var.vbc_types[6] = param.bc.vbc_n0;
    var.vbc_types[7] = param.bc.vbc_n1;
    var.vbc_types[8] = param.bc.vbc_n2;
    var.vbc_types[9] = param.bc.vbc_n3;

    var.vbc_values[0] = param.bc.vbc_val_x0;
    var.vbc_values[1] = param.bc.vbc_val_x1;
    var.vbc_values[2] = param.bc.vbc_val_y0;
    var.vbc_values[3] = param.bc.vbc_val_y1;
    var.vbc_values[4] = param.bc.vbc_val_z0;
    var.vbc_values[5] = param.bc.vbc_val_z1;
    var.vbc_values[6] = param.bc.vbc_val_n0;
    var.vbc_values[7] = param.bc.vbc_val_n1;
    var.vbc_values[8] = param.bc.vbc_val_n2;
    var.vbc_values[9] = param.bc.vbc_val_n3;
}


void init(const Param& param, Variables& var)
{
    BENCHMARK_START("Initialization");
    std::cout << "Initializing mesh and field data...\n";

    create_new_mesh(param, var);
    std::cout << "Mesh created\n";
    create_boundary_flags(var);
    std::cout << "Boundary flags created\n";
    create_boundary_nodes(var);
    std::cout << "Boundary nodes created\n";
    create_boundary_facets(var);
    std::cout << "Boundary facets created\n";
    create_support(var);
    std::cout << "Support created\n";
    create_elem_groups(var);
    std::cout << "Elem groups created\n";
    create_elemmarkers(param, var);
    std::cout << "Elem markers created\n";
    create_markers(param, var);
    std::cout << "Markers created\n";
    allocate_variables(param, var);
    std::cout << "Variables allocated\n";

    for(int i=0; i<var.nnode; i++)
        for(int d=0; d<NDIMS; d++)
            (*var.coord0)[i][d] = (*var.coord)[i][d];

    compute_volume(*var.coord, *var.connectivity, *var.volume);
    *var.volume_old = *var.volume;
    initial_material_properties(var, *var.rho);
    std::cout << "Material properties initialized\n";
    compute_mass(param, var.egroups, *var.connectivity, *var.volume, *var.mat,
                 var.max_vbc_val, *var.volume_n, *var.stressyy, *var.mass, *var.tmass, var);
    std::cout << "Mass computed\n";
    compute_shape_fn(*var.coord, *var.connectivity, *var.volume, var.egroups,
                     *var.shpdx, *var.shpdy, *var.shpdz);
    std::cout << "Shape fn computed\n";

    create_boundary_normals(var, var.bnormals, var.edge_vectors);
    apply_vbcs(param, var, *var.vel);
    std::cout << "VBCs applied\n";
    // temperature should be init'd before stress and strain
    initial_temperature(param, var, *var.temperature);
    std::cout << "Temperature initialized\n";
    initial_stress_state(param, var, *var.stress, *var.stressyy, *var.dP, *var.strain, var.compensation_pressure);
    std::cout << "Stress initialized\n";
    initial_weak_zone(param, var, *var.plstrain);
    std::cout << "Weak zone initialized\n";
    BENCHMARK_END();
}


void restart(const Param& param, Variables& var)
{
    std::cout << "Initializing mesh and field data from checkpoints...\n";

    /* Reading info file */
    {
        char filename[256];
        std::snprintf(filename, 255, "%s.info", param.sim.restarting_from_modelname.c_str());
        std::FILE *f = std::fopen(filename, "r");
        int frame, steps, nnode, nelem, nseg;

        while (1) {
            int n = std::fscanf(f, "%d %d %*f %*f %*f %d %d %d\n",
                                &frame, &steps, &nnode, &nelem, &nseg);

            if (n != 5) {
                std::cerr << "Error: reading info file: " << filename << '\n';
                std::exit(2);
            }
            if (frame == param.sim.restarting_from_frame)
                break;
        }



        var.steps = steps;
        var.nnode = nnode;
        var.nelem = nelem;
        var.nseg = nseg;

        std::fclose(f);
    }

    char filename_save[256];
    std::snprintf(filename_save, 255, "%s.save.%06d",
                  param.sim.restarting_from_modelname.c_str(), param.sim.restarting_from_frame);
    BinaryInput bin_save(filename_save);
    std::cout << "  Reading " << filename_save << "...\n";

    char filename_chkpt[256];
    std::snprintf(filename_chkpt, 255, "%s.chkpt.%06d",
                  param.sim.restarting_from_modelname.c_str(), param.sim.restarting_from_frame);
    BinaryInput bin_chkpt(filename_chkpt);
    std::cout << "  Reading " << filename_chkpt << "...\n";

    //
    // Following the same procedure in init()
    //

    // Reading mesh, replacing create_new_mesh()
    {
        var.coord = new array_t(var.nnode);
        bin_save.read_array(*var.coord, "coordinate");
        var.connectivity = new conn_t(var.nelem);
        bin_save.read_array(*var.connectivity, "connectivity");

        var.segment = new segment_t(var.nseg);
        bin_chkpt.read_array(*var.segment, "segment");
        var.segflag = new segflag_t(var.nseg);
        bin_chkpt.read_array(*var.segflag, "segflag");
        // Note: regattr is not needed for restarting
        // var.regattr = new regattr_t(var.nelem);
        // bin_chkpt.read_array(*var.regattr, "regattr", var.nelem);
    }

    create_boundary_flags(var);
    create_boundary_nodes(var);
    create_boundary_facets(var);
    create_support(var);
    create_elem_groups(var);
    create_elemmarkers(param, var);

    // Replacing create_markers()
    var.markersets.push_back(new MarkerSet(param, var, bin_chkpt, std::string("markerset")));
    if (param.control.has_hydration_processes) {
        var.hydrous_marker_index = var.markersets.size();
        var.markersets.push_back(new MarkerSet(param, var, bin_chkpt, std::string("hydrous-markerset")));
    }

    allocate_variables(param, var);

    // Initializing field variables
    bin_save.read_array(*var.coord0, "coord0");

    compute_volume(*var.coord, *var.connectivity, *var.volume);
    bin_chkpt.read_array(*var.volume_old, "volume_old");
    compute_mass(param, var.egroups, *var.connectivity, *var.volume, *var.mat,
                 var.max_vbc_val, *var.volume_n, *var.stressyy, *var.mass, *var.tmass, var);
    compute_shape_fn(*var.coord, *var.connectivity, *var.volume, var.egroups,
                     *var.shpdx, *var.shpdy, *var.shpdz);

    create_boundary_normals(var, var.bnormals, var.edge_vectors);
    apply_vbcs(param, var, *var.vel);

    {
        bin_save.read_array(*var.vel, "velocity");
        bin_save.read_array(*var.temperature, "temperature");
        bin_save.read_array(*var.strain_rate, "strain-rate");
        bin_save.read_array(*var.strain, "strain");
        bin_save.read_array(*var.stress, "stress");
        bin_save.read_array(*var.plstrain, "plastic strain");

        //======Energy Balance reletd parameters ================
        bin_chkpt.read_array(*var.tenergy, "total_energy");
        bin_chkpt.read_array(*var.venergy, "volumetric_energy");
        bin_chkpt.read_array(*var.denergy, "deviatoric_energy");
        bin_chkpt.read_array(*var.dP, "dP");
        bin_chkpt.read_array(*var.rho, "rho");
        bin_chkpt.read_array(*var.drho, "drho");
        //=======================================================

        if (param.mat.is_plane_strain)
            bin_chkpt.read_array(*var.stressyy, "stressyy");
    }

    compute_volume(*var.coord, *var.connectivity, *var.volume);
    bin_chkpt.read_array(*var.volume_old, "volume_old");
    compute_mass(param, var.egroups, *var.connectivity, *var.volume, *var.mat,
                 var.max_vbc_val, *var.volume_n, *var.stressyy, *var.mass, *var.tmass, var);
    compute_shape_fn(*var.coord, *var.connectivity, *var.volume, var.egroups,
                     *var.shpdx, *var.shpdy, *var.shpdz);

    apply_vbcs(param, var, *var.vel);

    // Misc. items
    {
        double_vec tmp(2);
        bin_chkpt.read_array(tmp, "time compensation_pressure");
        var.time = tmp[0];
        var.compensation_pressure = tmp[1];

        // the following fields are not required for restarting
        bin_save.read_array(*var.force, "force");
    }
}


void update_mesh(const Param& param, Variables& var)
{
    update_coordinate(var, *var.coord);
    surface_processes(param, var, *var.coord);

    var.volume->swap(*var.volume_old);
    compute_volume(*var.coord, *var.connectivity, *var.volume);
    compute_mass(param, var.egroups, *var.connectivity, *var.volume, *var.mat,
                 var.max_vbc_val, *var.volume_n, *var.stressyy, *var.mass, *var.tmass, var);
    compute_shape_fn(*var.coord, *var.connectivity, *var.volume, var.egroups,
                     *var.shpdx, *var.shpdy, *var.shpdz);
}


void isostasy_adjustment(const Param &param, Variables &var)
{
    std::cout << "Adjusting isostasy for " << param.ic.isostasy_adjustment_time_in_yr << " yrs...\n";

    var.dt = compute_dt(param, var);
    int iso_steps = param.ic.isostasy_adjustment_time_in_yr*YEAR2SEC / var.dt;

    for (int i=0; i<iso_steps; i++) {
        update_strain_rate(var, *var.strain_rate);
        compute_dvoldt(var, *var.ntmp);
        compute_edvoldt(var, *var.ntmp, *var.edvoldt);
        update_stress(var, *var.stress, *var.stressyy, *var.thermal_stress, *var.dP, *var.strain, *var.elastic_strain,
                      *var.plstrain, *var.delta_plstrain, *var.dtemp, *var.strain_rate, *var.power,
                      *var.tenergy, *var.venergy, *var.denergy);
        apply_NMD_to_Stress(var, *var.stress, *var.stressyy, *var.ediffStress, *var.ndiffStress, *var.dP);
        update_force(param, var, *var.force);
        update_velocity(var, *var.vel);

        // do not apply vbc to allow free boundary

        // displacment is vertical only
        #pragma omp parallel for default(none)          \
            shared(var, param)
        for (int i=0; i<var.nnode; ++i) {
            for (int j=0; j<NDIMS-1; ++j) {
                (*var.vel)[i][j] = 0;
            }
            if (param.bc.has_wrinkler_foundation == false &&
                (*var.bcflag)[i] & BOUNDZ0) {
                // holding bottom surface fixed
                (*var.vel)[i][NDIMS-1] = 0;
            }
        }

        update_mesh(param, var);

    }
    std::cout << "Adjusted isostasy for " << iso_steps << " steps.\n";
}



bool check_mesh_quality(const Param& param, Variables& var)
{
    int index;
    int quality = bad_mesh_quality(param, var, index);
    if (quality != 0) {
        remesh(param, var, quality);
        return true;
    }
    return false;
}

int main(int argc, char *argv[])
{
    std::ios::sync_with_stdio(false);

    // Start total timing
    BENCHMARK_START("Total Runtime");

    // OpenMP info
#ifdef USE_OMP
    int num_threads = omp_get_max_threads();
    std::cout << "=== OpenMP ENABLED: " << num_threads << " threads ===" << std::endl;
    
    // Set number of threads explicitly
    omp_set_num_threads(num_threads);
    
    // Verify in parallel region
    #pragma omp parallel
    {
        #pragma omp master
        {
            std::cout << "Parallel region active with " << omp_get_num_threads() << " threads" << std::endl;
        }
    }
#else
    std::cout << "=== OpenMP IS NOT ENABLED ===" << std::endl;
#endif

    // Parsing arguments
    Param param;
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " config_file\n";
        return 1;
    }
    get_input_parameters(argv[1], param);
    
    std::cout << "DEBUG: Loaded max_steps = " << param.sim.max_steps << std::endl;
    std::cout << "DEBUG: Loaded output_step_interval = " << param.sim.output_step_interval << std::endl;
    std::cout << "DEBUG: Loaded output_time_interval_in_yr = " << param.sim.output_time_interval_in_yr << std::endl;

    // Initialize GPU acceleration if available
#ifdef WITH_OPENCL
    if (param.control.use_gpu) {
        std::string device_type = param.control.gpu_device_type;
        if (initializeGPUAcceleration(device_type)) {
            std::cout << "GPU acceleration enabled: " << gpu_stress_update->getDeviceInfo() << std::endl;
        } else {
            std::cout << "Failed to initialize GPU acceleration, using CPU-only mode." << std::endl;
        }
    } else {
        std::cout << "GPU acceleration disabled by configuration." << std::endl;
    }
#endif

    // Variables related to FE mesh, defined in "mesh.hpp"
    Variables var;
    init_var(param, var);

    // Generating/retrieving mesh, allocating variables, setting up initial conditions
    if (param.sim.is_restarting) {
        restart(param, var);
    }
    else {
        init(param, var);
    }

    // Computing dt for the first time step
    double dt = compute_dt(param, var);
    std::cout << "Initial time step: " << dt << " seconds.\n";
    var.dt = dt;

    // Setting up the output manager
    double start_time = var.time;
    int start_frame = param.sim.is_restarting ? param.sim.restarting_from_frame : 0;
    Output output(param, start_time, start_frame);

    if (! param.sim.is_restarting) {
        // new simulation, save the initial condition
        output.write_checkpoint(param, var);
    }
    else {
        // copy time from checkpoint
        {
            char filename_save[256];
            std::snprintf(filename_save, 255, "%s.save.%06d",
                          param.sim.restarting_from_modelname.c_str(),
                          param.sim.restarting_from_frame);
            BinaryInput bin_save(filename_save);
            std::vector<double> tmp(1);
            bin_save.read_array(tmp, "time");
            var.time = tmp[0];
        }
        std::cout << "Restart from frame #" << param.sim.restarting_from_frame
                  << " of " << param.sim.restarting_from_modelname
                  << ". Time = " << var.time << " seconds.\n";
    }

    // Advancing the solution with explicit time-stepping scheme
    double next_output_time = var.time + param.sim.output_time_interval_in_yr * YEAR2SEC;
    while (var.time < param.sim.max_time_in_yr * YEAR2SEC) {
        BENCHMARK_START("Time Step");
        
        std::cout << "STEP: " << var.steps << "   TIME: " << var.time << "   dt: " << var.dt;
        if (param.sim.is_restarting)
            std::cout << "   *** restart ***";
        std::cout << '\n';

        ++var.steps;
        var.time += var.dt;

        try {
            BENCHMARK_START("Physics Update");
            update_temperature(param, var, *var.temperature,
                               *var.temp_power, *var.temp_pressure, *var.temp_density,
                               *var.dtemp, *var.dP, *var.ntmp, *var.stress,
                               *var.strain_rate, *var.stressyy, *var.drho, *var.rho,
                               *var.power, *var.powerTerm, *var.pressureTerm, *var.densityTerm);
            update_strain_rate(var, *var.strain_rate);
            update_stress(var, *var.stress, *var.stressyy, *var.thermal_stress, *var.dP,
                    *var.strain, *var.elastic_strain, *var.plstrain, *var.delta_plstrain, *var.dtemp, *var.strain_rate,
                    *var.power, *var.tenergy, *var.venergy, *var.denergy);
            phase_changes(param, var);
            BENCHMARK_END(); // End Physics Update

            BENCHMARK_START("Force & Velocity");
            update_force(param, var, *var.force);
            update_velocity(var, *var.vel);
            apply_vbcs(param, var, *var.vel);
            BENCHMARK_END(); // End Force & Velocity

            bool remesh = check_mesh_quality(param, var);
            if (remesh) {
                dt = compute_dt(param, var);
                var.dt = dt;
            }
            else {
                // No remeshing is done, update coordinate, volume, etc.
                update_coordinate(var, *var.coord);
                compute_volume(*var.coord, *var.connectivity, *var.volume);
                compute_shape_fn(*var.coord, *var.connectivity, *var.volume, var.egroups,
                                *var.shpdx, *var.shpdy, *var.shpdz);
            }

            // compute_mass() needs to be excuted before update_mesh()
            compute_mass(param, var.egroups, *var.connectivity, *var.volume, *var.mat,
                        var.max_vbc_val, *var.volume_n, *var.stressyy, *var.mass, *var.tmass, var);
        }
        catch (std::exception& e) {
            std::cout << "Exception at step " << var.steps
                      << ": " << e.what() << '\n';
            break;
        }

        param.sim.is_restarting = false;
        
        if (var.steps >= param.sim.max_steps) break;

        if (var.steps % param.sim.output_step_interval == 0 ||
            var.time >= next_output_time) {
            BENCHMARK_START("Output");
            output.write(var);
            BENCHMARK_END(); // End Output
            next_output_time += param.sim.output_time_interval_in_yr * YEAR2SEC;
        }

        // compute_dt needs to be excuted after mass, coord are updated
        dt = compute_dt(param, var);
        var.dt = dt;
        
        BENCHMARK_END(); // End Time Step
    }
    
    BENCHMARK_END(); // End Total Runtime
    
    // Print performance summary
    std::cout << "\n";
    
    // Write CSV for analysis
    std::string csv_filename = param.sim.modelname + "_timing.csv";
    g_benchmark.write_csv(csv_filename);
    
    // Cleanup GPU resources
#ifdef WITH_OPENCL
    finalizeGPUAcceleration();
#endif

    return 0;
}
