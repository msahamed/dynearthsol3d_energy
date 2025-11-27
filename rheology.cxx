#include <cmath>
#include <iostream>

#include "3x3-C/dsyevh3.h"

#include "constants.hpp"
#include "parameters.hpp"
#include "matprops.hpp"
#include "rheology.hpp"
#include "utils.hpp"


static void principal_stresses3(const double* s, double p[3], double v[3][3])
{
    /* s is a flattened stress vector, with the components {XX, YY, ZZ, XY, XZ, YZ}.
     * Returns the eigenvalues p and eignvectors v.
     * The eigenvalues are ordered such that p[0] <= p[1] <= p[2].
     */

    // unflatten s to a 3x3 tensor, only the upper part is needed.
    double a[3][3];
    a[0][0] = s[0];
    a[1][1] = s[1];
    a[2][2] = s[2];
    a[0][1] = s[3];
    a[0][2] = s[4];
    a[1][2] = s[5];

    dsyevh3(a, v, p);

    // reorder p and v
    if (p[0] > p[1]) {
        double tmp, b[3];
        tmp = p[0];
        p[0] = p[1];
        p[1] = tmp;
        for (int i=0; i<3; ++i)
            b[i] = v[i][0];
        for (int i=0; i<3; ++i)
            v[i][0] = v[i][1];
        for (int i=0; i<3; ++i)
            v[i][1] = b[i];
    }
    if (p[1] > p[2]) {
        double tmp, b[3];
        tmp = p[1];
        p[1] = p[2];
        p[2] = tmp;
        for (int i=0; i<3; ++i)
            b[i] = v[i][1];
        for (int i=0; i<3; ++i)
            v[i][1] = v[i][2];
        for (int i=0; i<3; ++i)
            v[i][2] = b[i];
    }
    if (p[0] > p[1]) {
        double tmp, b[3];
        tmp = p[0];
        p[0] = p[1];
        p[1] = tmp;
        for (int i=0; i<3; ++i)
            b[i] = v[i][0];
        for (int i=0; i<3; ++i)
            v[i][0] = v[i][1];
        for (int i=0; i<3; ++i)
            v[i][1] = b[i];
    }
}


static void principal_stresses2(const double* s, double p[2],
                                double& cos2t, double& sin2t)
{
    /* 's' is a flattened stress vector, with the components {XX, ZZ, XZ}.
     * Returns the eigenvalues 'p', and the direction cosine of the
     * eigenvectors in the X-Z plane.
     * The eigenvalues are ordered such that p[0] <= p[1].
     */

    // center and radius of Mohr circle
    double s0 = 0.5 * (s[0] + s[1]);
    double rad = second_invariant(s);

    // principal stresses in the X-Z plane
    p[0] = s0 - rad;
    p[1] = s0 + rad;

    {
        // direction cosine and sine of 2*theta
        const double eps = 1e-15;
        double a = 0.5 * (s[0] - s[1]);
        double b = - rad; // always negative
        if (b < -eps) {
            cos2t = a / b;
            sin2t = s[2] / b;
        }
        else {
            cos2t = 1;
            sin2t = 0;
        }
    }
}


static void elastic(double bulkm, double shearm, const double* de, double* s)
{
    /* increment the stress s according to the incremental strain de */
    double lambda = bulkm - 2. /3 * shearm;
    double dev = trace(de);

    for (int i=0; i<NDIMS; ++i)
        s[i] += 2 * shearm * de[i] + lambda * dev;
    for (int i=NDIMS; i<NSTR; ++i)
        s[i] += 2 * shearm * de[i];
}


static void maxwell(double bulkm, double shearm, double viscosity, double dt,
                    double dv, const double* de, double* s)
{
    // non-dimensional parameter: dt/ relaxation time
    double tmp = 0.5 * dt * shearm / viscosity;
    double f1 = 1 - tmp;
    double f2 = 1 / (1  + tmp);

    double dev = trace(de) / NDIMS;
    double s0 = trace(s) / NDIMS;

    // convert back to total stress
    for (int i=0; i<NDIMS; ++i)
        s[i] = ((s[i] - s0) * f1 + 2 * shearm * (de[i] - dev)) * f2 + s0 + bulkm * dv;
    for (int i=NDIMS; i<NSTR; ++i)
        s[i] = (s[i] * f1 + 2 * shearm * de[i]) * f2;
}


static void viscous(double bulkm, double viscosity, double total_dv,
                    const double* edot, double* s)
{
    /* Viscous Model + incompressibility enforced by bulk modulus */

    double dev = trace(edot) / NDIMS;

    for (int i=0; i<NDIMS; ++i)
        s[i] = 2 * viscosity * (edot[i] - dev) + bulkm * total_dv;
    for (int i=NDIMS; i<NSTR; ++i)
        s[i] = 2 * viscosity * edot[i];
}

static void compute_elastic_strain(const double* total_dstrain, const double depls1, const double depls2, double* estrain){
    double cos2t, sin2t;
    double p[2], des[2];
    principal_stresses2(total_dstrain, p, cos2t, sin2t);
    des[0] = p[0] - depls1;
    des[1] = p[1] - depls2;

    // principal strain to tensor
    double dc2 = (des[0] - des[1]) * cos2t;
    double dss = des[0] + des[1];
    estrain[0] += 0.5 * (dss + dc2);
    estrain[1] += 0.5 * (dss - dc2);
    estrain[2] += 0.5 * (des[0] - des[1]) * sin2t;
}

static void elasto_plastic(double bulkm, double shearm,
                           double amc, double anphi, double anpsi,
                           double hardn, double ten_max,
                           const double* de, double& depls, double* s,
                           int &failure_mode)
{
    /* Elasto-plasticity (Mohr-Coulomb criterion)
     *
     * failure_mode --
     *   0: no failure
     *   1: tensile failure
     *  10: shear failure
     */

    // elastic trial stress
    elastic(bulkm, shearm, de, s);
    depls = 0;
    failure_mode = 0;

    //
    // transform to principal stress coordinate system
    //
    // eigenvalues (principal stresses)
    double p[NDIMS];
#ifdef THREED
    // eigenvectors
    double v[3][3];
    principal_stresses3(s, p, v);
#else
    // In 2D, we only construct the eigenvectors from
    // cos(2*theta) and sin(2*theta) of Mohr circle
    double cos2t, sin2t;
    principal_stresses2(s, p, cos2t, sin2t);
#endif

    // composite (shear and tensile) yield criterion
    double fs = p[0] - p[NDIMS-1] * anphi + amc;
    double ft = p[NDIMS-1] - ten_max;

    if (fs > 0 && ft < 0) {
        // no failure
        return;
    }

    // yield, shear or tensile?
    double pa = std::sqrt(1 + anphi*anphi) + anphi;
    double ps = ten_max * anphi - amc;
    double h  = p[NDIMS-1] - ten_max + pa * (p[0] - ps);
    double a1 = bulkm + 4. / 3 * shearm;
    double a2 = bulkm - 2. / 3 * shearm;

    double alam;
    if (h < 0) {
        // shear failure
        failure_mode = 10;

        alam = fs / (a1 - a2*anpsi + a1*anphi*anpsi - a2*anphi + 2*std::sqrt(anphi)*hardn);
        p[0] -= alam * (a1 - a2 * anpsi);
#ifdef THREED
        p[1] -= alam * (a2 - a2 * anpsi);
#endif
        p[NDIMS-1] -= alam * (a2 - a1 * anpsi);

        // 2nd invariant of plastic strain
#ifdef THREED
        /* // plastic strain in the principle directions, depls2 is always 0
         * double depls1 = alam;
         * double depls3 = -alam * anpsi;
         * double deplsm = (depls1 + depls3) / 3;
         * depls = std::sqrt(((depls1-deplsm)*(depls1-deplsm) +
         *                    (-deplsm)*(-deplsm) +
         *                    (depls3-deplsm)*(depls3-deplsm) +
         *                    deplsm*deplsm) / 2);
         */
        // the equations above can be reduce to:
        depls = std::fabs(alam) * std::sqrt((7 + 4*anpsi + 7*anpsi*anpsi) / 18);
#else
        /* // plastic strain in the principle directions
         * double depls1 = alam;
         * double depls2 = -alam * anpsi;
         * double deplsm = (depls1 + depls2) / 2;
         * depls = std::sqrt(((depls1-deplsm)*(depls1-deplsm) +
         *                    (depls2-deplsm)*(depls2-deplsm) +
         *                    deplsm*deplsm) / 2);
         */
        // the equations above can be reduce to:
        depls = std::fabs(alam) * std::sqrt((3 + 2*anpsi + 3*anpsi*anpsi) / 8);
#endif
    }
    else {
        // tensile failure
        failure_mode = 1;

        alam = ft / a1;
        p[0] -= alam * a2;
#ifdef THREED
        p[1] -= alam * a2;
#endif
        p[NDIMS-1] -= alam * a1;

        // 2nd invariant of plastic strain
#ifdef THREED
        /* double depls1 = 0;
         * double depls3 = alam;
         * double deplsm = (depls1 + depls3) / 3;
         * depls = std::sqrt(((depls1-deplsm)*(depls1-deplsm) +
         *                    (-deplsm)*(-deplsm) +
         *                    (depls3-deplsm)*(depls3-deplsm) +
         *                    deplsm*deplsm) / 2);
         */
        depls = std::fabs(alam) * std::sqrt(7. / 18);
#else
        /* double depls1 = 0;
         * double depls3 = alam;
         * double deplsm = (depls1 + depls3) / 2;
         * depls = std::sqrt(((depls1-deplsm)*(depls1-deplsm) +
         *                    (depls2-deplsm)*(depls2-deplsm) +
         *                    deplsm*deplsm) / 2);
         */
        depls = std::fabs(alam) * std::sqrt(3. / 8);
#endif
    }

    // rotate the principal stresses back to global axes
    {
#ifdef THREED
        double ss[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
        for(int m=0; m<3; m++) {
            for(int n=m; n<3; n++) {
                for(int k=0; k<3; k++) {
                    ss[m][n] += v[m][k] * v[n][k] * p[k];
                }
            }
        }
        s[0] = ss[0][0];
        s[1] = ss[1][1];
        s[2] = ss[2][2];
        s[3] = ss[0][1];
        s[4] = ss[0][2];
        s[5] = ss[1][2];
#else
        double dc2 = (p[0] - p[1]) * cos2t;
        double dss = p[0] + p[1];
        s[0] = 0.5 * (dss + dc2);
        s[1] = 0.5 * (dss - dc2);
        s[2] = 0.5 * (p[0] - p[1]) * sin2t;
#endif
    }
}


static void elasto_plastic2d(double bulkm, double shearm, double& t_power,
                             double& v_power, double& d_power, double amc,
                             double anphi, double anpsi, double hardn, double ten_max,
                             const double* de, double& depls, double* s, 
                             double &syy, int &failure_mode, double thermal_stress, double* estrain)
{
    /* Elasto-plasticity (Mohr-Coulomb criterion) */

    /* This function is derived from geoFLAC.
     * The original code in geoFLAC assumes 2D plain strain formulation,
     * i.e. there are 3 principal stresses (PSs) and only 2 principal strains
     * (Strain_yy, Strain_xy, Strain_yz must be 0).
     * Here, the code is adopted to pure 2D or 3D plane strain.
     *
     * failure_mode --
     *   0: no failure
     *   1: tensile failure, all PSs exceed tensile limit
     *   2: tensile failure, 2 PSs exceed tensile limit
     *   3: tensile failure, 1 PS exceeds tensile limit
     *  10: pure shear failure
     *  11, 12, 13: tensile + shear failure
     *  20, 21, 22, 23: shear + tensile failure
     */

    depls = 0;
    failure_mode = 0;

    // elastic trial stress
    // double thermalStress = (1  * bulkm * alpha * dT);
    double a1 = bulkm + 4. / 3 * shearm;
    double a2 = bulkm - 2. / 3 * shearm;
    double sxx = s[0] + (de[1] * a2) + (de[0] * a1) + thermal_stress;
    double szz = s[1] + (de[0] * a2) + (de[1] * a1) + thermal_stress;
    double sxz = s[2] + de[2]  * 2 * shearm;
    syy += ((de[0] + de[1]) * a2) + thermal_stress; // Stress YY component, plane strain


    // transform to principal stress coordinate system
    //
    // eigenvalues (principal stresses)
    double p[3];

    // In 2D, we only construct the eigenvectors from
    // cos(2*theta) and sin(2*theta) of Mohr circle
    double cos2t, sin2t;
    int n1, n2, n3;

    {
        // center and radius of Mohr circle
        double s0 = 0.5 * (sxx + szz);
        double rad = 0.5 * std::sqrt((sxx-szz)*(sxx-szz) + 4*sxz*sxz);

        // principal stresses in the X-Z plane
        double si = s0 - rad;
        double sii = s0 + rad;

        // direction cosine and sine of 2*theta
        const double eps = 1e-15;
        if (rad > eps) {
            cos2t = 0.5 * (szz - sxx) / rad;
            sin2t = -sxz / rad;
        }
        else {
            cos2t = 1;
            sin2t = 0;
        }

        // sort p.s.
#if 1
        //
        // 3d plane strain
        //
        if (syy > sii) {
            // syy is minor p.s.
            n1 = 0;
            n2 = 1;
            n3 = 2;
            p[0] = si;
            p[1] = sii;
            p[2] = syy;
        }
        else if (syy < si) {
            // syy is major p.s.
            n1 = 1;
            n2 = 2;
            n3 = 0;
            p[0] = syy;
            p[1] = si;
            p[2] = sii;
        }
        else {
            // syy is intermediate
            n1 = 0;
            n2 = 2;
            n3 = 1;
            p[0] = si;
            p[1] = syy;
            p[2] = sii;
        }
#else
        /* XXX: This case gives unreasonable result. Don't know why... */

        //
        // pure 2d case
        //
        n1 = 0;
        n2 = 2;
        p[0] = si;
        p[1] = syy;
        p[2] = sii;
#endif
    }

    // Possible tensile failure scenarios
    // 1. S1 (least tensional or greatest compressional principal stress) > ten_max:
    //    all three principal stresses must be greater than ten_max.
    //    Assign ten_max to all three and don't do anything further.
    if( p[0] >= ten_max ) {
        s[0] = s[1] = syy = ten_max;
        s[2] = 0.0;
        failure_mode = 1;
        return;
    }

    // 2. S2 (intermediate principal stress) > ten_max:
    //    S2 and S3 must be greater than ten_max.
    //    Assign ten_max to these two and continue to the shear failure block.
    if( p[1] >= ten_max ) {
        p[1] = p[2] = ten_max;
        failure_mode = 2;
    }

    // 3. S3 (most tensional or least compressional principal stress) > ten_max:
    //    Only this must be greater than ten_max.
    //    Assign ten_max to S3 and continue to the shear failure block.
    else if( p[2] >= ten_max ) {
        p[2] = ten_max;
        failure_mode = 3;
    }


    // shear yield criterion
    double fs = p[0] - p[2] * anphi + amc;
    if (fs >= 0.0) {
        // Tensile failure case S2 or S3 could have happened!!
        // XXX: Need to rationalize why exit without doing anything further.
        s[0] = sxx;
        s[1] = szz;
        s[2] = sxz;
        return;
    }

    failure_mode += 10;

    // shear failure
    const double alams = fs / (a1 - a2*anpsi + a1*anphi*anpsi - a2*anphi + hardn);
    p[0] -= alams * (a1 - a2 * anpsi);
    p[1] -= alams * (a2 - a2 * anpsi);
    p[2] -= alams * (a2 - a1 * anpsi);


    //Principal Plastic strain
    double depls1 = alams;
    double depls3 = -alams * anpsi;

    double deplsm = (depls1 + depls3) / 3;
    double vstress = (p[0]+p[1]+p[2])/3;

    compute_elastic_strain(de, depls1, depls3, estrain);

    // 2nd invariant of plastic strain
    depls = 0.5 * std::fabs(alams + alams * anpsi);
    // plasticpower = (p[0] * depls1) + (p[2] * depls3);

    t_power = (p[0] * depls1) + (p[2] * depls3);
    v_power = (vstress * deplsm);
    d_power = ((p[0]-vstress) * (depls1-deplsm)) + ((p[1]-vstress)*(-deplsm)) + ((p[2]-vstress) * (depls3-deplsm));
    //r_power = t_power - v_power - d_power;

    //***********************************
    // The following seems redundant but... this is how it goes in geoFLAC.
    //
    // Possible tensile failure scenarios
    // 1. S1 (least tensional or greatest compressional principal stress) > ten_max:
    //    all three principal stresses must be greater than ten_max.
    //    Assign ten_max to all three and don't do anything further.
    if( p[0] >= ten_max ) {
        s[0] = s[1] = syy = ten_max;
        s[2] = 0.0;
        failure_mode += 20;
        return;
    }

    // 2. S2 (intermediate principal stress) > ten_max:
    //    S2 and S3 must be greater than ten_max.
    //    Assign ten_max to these two and continue to the shear failure block.
    if( p[1] >= ten_max ) {
        p[1] = p[2] = ten_max;
        failure_mode += 20;
    }

    // 3. S3 (most tensional or least compressional principal stress) > ten_max:
    //    Only this must be greater than ten_max.
    //    Assign ten_max to S3 and continue to the shear failure block.
    else if( p[2] >= ten_max ) {
        p[2] = ten_max;
        failure_mode += 20;
    }
    //***********************************


    // rotate the principal stresses back to global axes
    {
      double dc2 = (p[n1] - p[n2]) * cos2t;
      double dss = p[n1] + p[n2];
      s[0] = 0.5 * (dss + dc2);
      s[1] = 0.5 * (dss - dc2);
      s[2] = 0.5 * (p[n1] - p[n2]) * sin2t;
      syy = p[n3];
    }
}

// GPU-friendly data structures for stress update
struct StressUpdateParams {
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
    
    // Input arrays (read-only)
    double de[NSTR];
    double edot[NSTR];
    double s_in[NSTR];
    double es_in[NSTR];
    double syy_in;
    double plstrain_in;
    
    // Output arrays (write-only)
    double s_out[NSTR];
    double es_out[NSTR];
    double syy_out;
    double depls_out;
    double power[3]; // thermal, viscous, dissipative
    double pressure_change;
    int failure_mode;
};

// Compute kernel suitable for GPU offloading - handles a single element's stress update
static void compute_stress_update_kernel(StressUpdateParams& params) 
{
    // Local copies for better register usage on GPU
    double* s = params.s_out;
    double& syy = params.syy_out;
    double* es = params.es_out;
    double* de = params.de;
    
    // Initialize output with input values
    for (int i = 0; i < NSTR; i++) {
        s[i] = params.s_in[i];
        es[i] = params.es_in[i];
    }
    syy = params.syy_in;
    
    // Calculate thermal stress
    double tstress = -params.bulkm * params.alpha * params.dT;
    double t_power = 0;
    double v_power = 0;
    double d_power = 0;
    
    // Compute pressure change before stress update
#ifdef THREED
    double pressure_old = -(s[0] + s[1] + s[2]) / NDIMS;
#else
    double pressure_old = -(s[0] + s[1] + syy) / 3.0;
#endif
    
    // Apply the rheology model based on element type
    switch (params.rheol_type) {
    case MatProps::rh_elastic:
        elastic(params.bulkm, params.shearm, de, s);
        break;
    case MatProps::rh_viscous:
        {
            double total_dv = trace(de);
            viscous(params.bulkm, params.viscosity, total_dv, params.edot, s);
        }
        break;
    case MatProps::rh_maxwell:
        maxwell(params.bulkm, params.shearm, params.viscosity, params.dt, params.dv, de, s);
        break;
    case MatProps::rh_ep:
        {
            if (params.dT != 0) {
                tstress = -params.bulkm * params.alpha * params.dT;
                // Store thermal_stress update in t_power temporarily
                t_power = tstress;
            }
            
            if (true) { // Always apply plane strain for now 
                elasto_plastic2d(params.bulkm, params.shearm, t_power, v_power, d_power, 
                                 params.amc, params.anphi, params.anpsi, params.hardn, 
                                 params.ten_max, de, params.depls_out, s, syy, 
                                 params.failure_mode, tstress, es);
            }
            else {
                elasto_plastic(params.bulkm, params.shearm, params.amc, params.anphi, 
                               params.anpsi, params.hardn, params.ten_max,
                               de, params.depls_out, s, params.failure_mode);
            }
        }
        break;
    default:
        // Unknown rheology type - no update
        break;
    }
    
    // Compute final pressure and pressure change
#ifdef THREED
    double pressure_new = -(s[0] + s[1] + s[2]) / NDIMS;
#else
    double pressure_new = -(s[0] + s[1] + syy) / 3.0;
#endif
    params.pressure_change = pressure_new - pressure_old;
    
    // Store power values
    params.power[0] = t_power; // thermal
    params.power[1] = v_power; // viscous
    params.power[2] = d_power; // dissipative
}

void update_stress(const Variables& var, tensor_t& stress, double_vec& stressyy, double_vec& thermal_stress, double_vec& dP,
                   tensor_t& strain, tensor_t& elastic_strain, double_vec& plstrain, double_vec& delta_plstrain, double_vec& dtemp,
                   tensor_t& strain_rate, double_vec& power, double_vec& tenergy, double_vec& venergy, double_vec& denergy)
{
    const int rheol_type = var.mat->rheol_type;
    const size_t nelem = var.nelem;
    
    // First pass: prepare input data for all elements
    std::vector<StressUpdateParams> elem_params(nelem);
    
    #pragma omp parallel for default(none) \
        shared(var, stress, stressyy, elem_params, strain_rate, strain, dtemp, elastic_strain, plstrain, rheol_type, nelem)
    for (int e=0; e<nelem; ++e) {
        StressUpdateParams& params = elem_params[e];
        
        // Compute average temperature change
        double dT = 0;
        const int *conn = (*var.connectivity)[e];
        for (int i = 0; i < NODES_PER_ELEM; ++i) {
            dT += dtemp[conn[i]];
        }
        params.dT = dT / NODES_PER_ELEM;
        
        // Initialize material parameters
        params.rheol_type = rheol_type;
        params.dt = var.dt;
        params.bulkm = var.mat->bulkm(e);
        params.shearm = var.mat->shearm(e);
        params.viscosity = var.mat->visc(e);
        params.alpha = var.mat->alpha(e);
        params.dv = (*var.volume)[e] / (*var.volume_old)[e] - 1;
        params.plstrain_in = plstrain[e];
        
        // If plastic model, get plastic properties
        if (rheol_type == MatProps::rh_ep) {
            var.mat->plastic_props(e, plstrain[e],
                                  params.amc, params.anphi, params.anpsi, 
                                  params.hardn, params.ten_max);
        }
        
        // Copy input stress/strain arrays
        double* s = stress[e];
        double& syy = stressyy[e];
        double* es = strain[e];
        double* edot = strain_rate[e];
        double* estrain = elastic_strain[e];
        
        for (int i = 0; i < NSTR; i++) {
            params.s_in[i] = s[i];
            params.es_in[i] = estrain[i];
            
            // Compute strain increment
            params.edot[i] = edot[i];
            params.de[i] = edot[i] * var.dt;
            
            // Apply anti-mesh locking correction if i < NDIMS
            if (i < NDIMS && rheol_type != MatProps::rh_viscous) {
                double div = trace(edot);
                params.edot[i] += ((*var.edvoldt)[e] - div) / NDIMS;
                params.de[i] = params.edot[i] * var.dt;
            }
        }
        params.syy_in = syy;
    }
    
    // This section would be offloaded to GPU in final implementation
    // Here we're just isolating the computation kernel
    #pragma omp parallel for default(none) \
        shared(elem_params, nelem)
    for (int e=0; e<nelem; ++e) {
        compute_stress_update_kernel(elem_params[e]);
    }
    
    // Update global arrays with results
    #pragma omp parallel for default(none) \
        shared(var, stress, stressyy, thermal_stress, dP, strain, elastic_strain, plstrain, delta_plstrain, \
               strain_rate, power, tenergy, venergy, denergy, elem_params, nelem)
    for (int e=0; e<nelem; ++e) {
        const StressUpdateParams& params = elem_params[e];
        
        // Copy results back
        double* s = stress[e];
        double& syy = stressyy[e];
        double* estrain = elastic_strain[e];
        
        for (int i = 0; i < NSTR; i++) {
            s[i] = params.s_out[i];
            estrain[i] = params.es_out[i];
        }
        syy = params.syy_out;
        
        // Update other result fields
        dP[e] = params.pressure_change;
        
        if (params.rheol_type == MatProps::rh_ep) {
            delta_plstrain[e] = params.depls_out;
            plstrain[e] += params.depls_out;
            
            // Update energy components
            tenergy[e] = params.power[0]; 
            venergy[e] = params.power[1];
            denergy[e] = params.power[2];
            power[e] = params.power[0] + params.power[1] + params.power[2];
            
            // Update thermal stress
            if (params.dT != 0) {
                thermal_stress[e] += -params.bulkm * params.alpha * params.dT;
            }
        }
    }
}
