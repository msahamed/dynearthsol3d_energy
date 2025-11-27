// OpenCL kernel for stress update computation
// This kernel is designed to run one work item per mesh element

// Constants matching constants.hpp
#ifndef NDIMS
#ifdef THREED
#define NDIMS 3
#else
#define NDIMS 2
#endif
#endif

#ifndef NSTR
#define NSTR (NDIMS * (NDIMS + 1) / 2)
#endif

// Rheology types (must match MatProps definitions)
#define RH_ELASTIC 0
#define RH_VISCOUS 1
#define RH_MAXWELL 2
#define RH_EP 3
#define RH_EVP 4

// Input structure that matches StressUpdateParams in C++ code
typedef struct {
    // Material properties
    double bulkm;
    double shearm;
    double viscosity;
    double alpha;
    double amc;
    double anphi; 
    double anpsi;
    double hardn;
    double ten_max;
    
    // Element properties
    int rheol_type;
    double dt;
    double dT;
    double dv;
    
    // Input arrays
    double de[NSTR];
    double edot[NSTR];
    double s_in[NSTR];
    double es_in[NSTR];
    double syy_in;
    double plstrain_in;
} StressParams;

// Output structure
typedef struct {
    // Output arrays
    double s_out[NSTR];
    double es_out[NSTR];
    double syy_out;
    double depls_out;
    double power[3]; // thermal, viscous, dissipative
    double pressure_change;
    int failure_mode;
} StressResults;

// Helper functions
inline double trace(const double* a) {
    double result = 0.0;
    for (int i = 0; i < NDIMS; i++) {
        result += a[i];
    }
    return result;
}

inline double second_invariant(const double* s) {
#ifdef THREED
    double s11 = s[0];
    double s22 = s[1];
    double s33 = s[2];
    double s12 = s[3];
    double s13 = s[4];
    double s23 = s[5];
    
    return sqrt(0.5 * ((s11-s22)*(s11-s22) + (s22-s33)*(s22-s33) + (s33-s11)*(s33-s11) +
                        6.0 * (s12*s12 + s23*s23 + s13*s13)));
#else
    double s11 = s[0];
    double s22 = s[1];
    double s12 = s[2];
    
    return sqrt(0.25 * (s11-s22)*(s11-s22) + s12*s12);
#endif
}

// Apply elastic rheology model
void elastic_kernel(double bulkm, double shearm, const double* de, double* s) {
    double lambda = bulkm - 2.0 / 3.0 * shearm;
    double dev = trace(de);

    for (int i = 0; i < NDIMS; ++i)
        s[i] += 2.0 * shearm * de[i] + lambda * dev;
    for (int i = NDIMS; i < NSTR; ++i)
        s[i] += 2.0 * shearm * de[i];
}

// Apply viscous rheology model
void viscous_kernel(double bulkm, double viscosity, double total_dv, const double* edot, double* s) {
    double dev = trace(edot) / NDIMS;

    for (int i = 0; i < NDIMS; ++i)
        s[i] = 2.0 * viscosity * (edot[i] - dev) + bulkm * total_dv;
    for (int i = NDIMS; i < NSTR; ++i)
        s[i] = 2.0 * viscosity * edot[i];
}

// Apply Maxwell rheology model
void maxwell_kernel(double bulkm, double shearm, double viscosity, double dt, double dv, const double* de, double* s) {
    // non-dimensional parameter: dt/ relaxation time
    double tmp = 0.5 * dt * shearm / viscosity;
    double f1 = 1.0 - tmp;
    double f2 = 1.0 / (1.0 + tmp);

    double dev = trace(de) / NDIMS;
    double s0 = trace(s) / NDIMS;

    // convert back to total stress
    for (int i = 0; i < NDIMS; ++i)
        s[i] = ((s[i] - s0) * f1 + 2.0 * shearm * (de[i] - dev)) * f2 + s0 + bulkm * dv;
    for (int i = NDIMS; i < NSTR; ++i)
        s[i] = (s[i] * f1 + 2.0 * shearm * de[i]) * f2;
}

// Main kernel for stress update
__kernel void stress_update(
    __global const StressParams* input,
    __global StressResults* output
) {
    // Get global element index
    int e = get_global_id(0);
    
    // Create local copies for better performance
    StressParams params = input[e];
    StressResults results;
    
    // Initialize output arrays from input
    for (int i = 0; i < NSTR; i++) {
        results.s_out[i] = params.s_in[i];
        results.es_out[i] = params.es_in[i];
    }
    results.syy_out = params.syy_in;
    results.depls_out = 0.0;
    results.failure_mode = 0;
    results.power[0] = results.power[1] = results.power[2] = 0.0;
    
    // Calculate thermal stress
    double tstress = -params.bulkm * params.alpha * params.dT;
    
    // Compute pressure change before stress update
#ifdef THREED
    double pressure_old = -(results.s_out[0] + results.s_out[1] + results.s_out[2]) / NDIMS;
#else
    double pressure_old = -(results.s_out[0] + results.s_out[1] + results.syy_out) / 3.0;
#endif
    
    // Local variables for convenience
    double* s = results.s_out;
    double* de = params.de;
    
    // Apply the rheology model based on element type
    switch (params.rheol_type) {
    case RH_ELASTIC:
        elastic_kernel(params.bulkm, params.shearm, de, s);
        break;
    case RH_VISCOUS:
        {
            double total_dv = trace(de);
            viscous_kernel(params.bulkm, params.viscosity, total_dv, params.edot, s);
        }
        break;
    case RH_MAXWELL:
        maxwell_kernel(params.bulkm, params.shearm, params.viscosity, params.dt, params.dv, de, s);
        break;
    case RH_EP:
        // Elasto-plastic model is complex and will be implemented
        // as a CPU fallback for now
        break;
    default:
        // Unknown rheology type - no update
        break;
    }
    
    // Compute final pressure and pressure change
#ifdef THREED
    double pressure_new = -(s[0] + s[1] + s[2]) / NDIMS;
#else
    double pressure_new = -(s[0] + s[1] + results.syy_out) / 3.0;
#endif
    results.pressure_change = pressure_new - pressure_old;
    
    // Write results back to global memory
    output[e] = results;
} 