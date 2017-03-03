#include <iostream>
#include "constants.hpp"
#include "parameters.hpp"
#include "bc.hpp"
#include "matprops.hpp"
#include "utils.hpp"
#include "fields.hpp"


void allocate_variables(const Param &param, Variables& var)
{
    const int n = var.nnode;
    const int e = var.nelem;
    const int tn = (var.bnodes[5]).size(); // Top surface node numbers

    var.volume = new double_vec(e);
    var.volume_old = new double_vec(e);
    var.volume_n = new double_vec(n);

    var.mass = new double_vec(n);
    var.tmass = new double_vec(n);

    var.edvoldt = new double_vec(e);


    {
        // these fields are reallocated during remeshing interpolation
        var.temperature = new double_vec(n);
        var.temp_power=new double_vec(n);
        var.temp_pressure = new double_vec(n);
        var.temp_density = new double_vec(n);
        var.coord0 = new array_t(n);
        var.plstrain = new double_vec(e);
        var.delta_plstrain = new double_vec(e);
        var.vel = new array_t(n, 0);
        var.topography = new array_t(tn, 0);
        var.strain = new tensor_t(e, 0);
        var.stress = new tensor_t(e, 0);
        var.stressyy = new double_vec(e, 0);
        var.thermal_stress = new double_vec(e, 0);
        var.ediffStress = new double_vec(e);
        var.ndiffStress = new double_vec(n);

        // Energy balance equation related variables :
        var.dtemp = new double_vec(n); // Temperature difference
        var.dP    = new double_vec(e);    // Pressure difference
        var.rho   = new double_vec(e);   // density
        var.drho  = new double_vec(e);   // density difference
        var.power    = new double_vec(e);
        var.tenergy    = new double_vec(e);
        var.venergy    = new double_vec(e);
        var.denergy    = new double_vec(e);
        var.powerTerm  = new double_vec(n);
        var.pressureTerm     = new double_vec(n);
        var.densityTerm      = new double_vec(n);

    }

    var.ntmp= new double_vec(n);

    var.force = new array_t(n, 0);

    var.strain_rate = new tensor_t(e, 0);

    var.shpdx = new shapefn(e);
    if (NDIMS == 3) var.shpdy = new shapefn(e);
    var.shpdz = new shapefn(e);

    var.mat = new MatProps(param, var);
}


void reallocate_variables(const Param& param, Variables& var)
{
    const int n = var.nnode;
    const int e = var.nelem;

    delete var.volume;
    delete var.volume_old;
    delete var.volume_n;
    var.volume     = new double_vec(e);
    var.volume_old = new double_vec(e);
    var.volume_n   = new double_vec(n);

    delete var.mass;
    delete var.tmass;
    var.mass  = new double_vec(n);
    var.tmass = new double_vec(n);

    delete var.edvoldt;
    var.edvoldt = new double_vec(e);

    delete var.ediffStress;
    var.ediffStress = new double_vec(e);

    delete var.ndiffStress;
    var.ndiffStress = new double_vec(n);

    delete var.ntmp;
    var.ntmp = new double_vec(n);

    delete var.force;
    var.force = new array_t(n, 0);

    delete var.strain_rate;
    var.strain_rate = new tensor_t(e, 0);

    delete var.shpdx;
    delete var.shpdz;
    var.shpdx = new shapefn(e);
    if (NDIMS == 3) {
        delete var.shpdy;
        var.shpdy = new shapefn(e);
    }
    var.shpdz = new shapefn(e);

    delete var.mat;
    var.mat = new MatProps(param, var);
}


void update_temperature(const Param &param, const Variables &var, double_vec &temperature,
                        double_vec &temp_power, double_vec &temp_pressure, double_vec &temp_density,
                        double_vec &dtemp,double_vec &dP, double_vec &tdot, tensor_t &stress,
                        tensor_t &strain_rate, double_vec &stressyy, double_vec &drho, double_vec &rho,
                        double_vec &power, double_vec &powerTerm, double_vec &pressureTerm, double_vec &densityTerm)
{
    tdot.assign(var.nnode, 0);
    powerTerm.assign(powerTerm.size(),0);
    pressureTerm.assign(pressureTerm.size(),0);
    densityTerm.assign(densityTerm.size(), 0);
    class ElemFunc_temperature : public ElemFunc
    {
    private:
        const Variables  &var;
        const double_vec &temperature;
              double_vec &temp_power;
              double_vec &temp_pressure;
              double_vec &temp_density;
              double_vec &dtemp;
              double_vec &dP;
              double_vec &drho;
              double_vec &rho;
              tensor_t   &stress;
              tensor_t   &strain_rate;
              double_vec &stressyy;
              double_vec &power;
              double_vec &powerTerm;
              double_vec &pressureTerm;
              double_vec &densityTerm;
              double_vec &tdot;
    public:
        ElemFunc_temperature(const Variables &var, const double_vec &temperature, double_vec &temp_power,
                             double_vec &temp_pressure, double_vec &temp_density, double_vec &dtemp,
                             double_vec &dP, tensor_t &stress, tensor_t   &strain_rate, double_vec &stressyy,
                             double_vec &drho, double_vec &rho, double_vec &power, double_vec &powerTerm,
                             double_vec &pressureTerm, double_vec &densityTerm, double_vec &tdot) :
                             var(var), temperature(temperature), temp_power(temp_power), temp_pressure(temp_pressure),
                             temp_density(temp_density), dtemp(dtemp), dP(dP), stress(stress), strain_rate(strain_rate),
                             stressyy(stressyy), drho(drho), rho(rho), power(power),pressureTerm(pressureTerm),
                             powerTerm(powerTerm), densityTerm(densityTerm),tdot(tdot) {};
        void operator()(int e)
        {
            // diffusion matrix
            double D[NODES_PER_ELEM][NODES_PER_ELEM];

            const int *conn = (*var.connectivity)[e];
            double kv = var.mat->k(e) *  (*var.volume)[e]; // thermal conductivity * volumn
            double* s = stress[e];
            double* edot = strain_rate[e];
            const double *shpdx = (*var.shpdx)[e];
    #ifdef THREED
            const double *shpdy = (*var.shpdy)[e];
    #endif
            const double *shpdz = (*var.shpdz)[e];
            for (int i=0; i<NODES_PER_ELEM; ++i) {
                for (int j=0; j<NODES_PER_ELEM; ++j) {
    #ifdef THREED
                    D[i][j] = (shpdx[i] * shpdx[j] +
                               shpdy[i] * shpdy[j] +
                               shpdz[i] * shpdz[j]);
    #else
                    D[i][j] = (shpdx[i] * shpdx[j] +
                               shpdz[i] * shpdz[j]);
    #endif
                }
            }
            double temp = 0;
            for (int i = 0; i < NODES_PER_ELEM; ++i) {
                temp += temperature[conn[i]];
            }

            double& syy = stressyy[e];

    #ifdef THREED
            double P = -(s[0] + s[1] + s[2]) / NDIMS;
            double vedot = edot[0]+ edot[1] + edot[2];
    #else
            double P = -(s[0] + s[1] + syy) / 3;
            double vedot = edot[0]+ edot[1];
    #endif
            double alpha    = var.mat->alpha(e);
            double T = temp / NODES_PER_ELEM;
            double plastic  = (*var.power)[e] * (*var.volume)[e]/NODES_PER_ELEM;
            double pressure = T * alpha * (*var.dP)[e] * (*var.volume)[e] / NODES_PER_ELEM;
            double den      = P * T * alpha * vedot * (*var.volume)[e] / NODES_PER_ELEM;

            for (int i = 0; i < NODES_PER_ELEM; ++i) {
                double diffusion = 0;
                powerTerm[conn[i]]        += plastic;
                pressureTerm[conn[i]]     += pressure;
                densityTerm[conn[i]]      += den;

                for (int j = 0; j<NODES_PER_ELEM; ++j)
                    diffusion += D[i][j] * temperature[conn[j]];

                tdot[conn[i]] += diffusion * kv;
            }
        }
    } elemf(var, temperature, temp_power, temp_pressure, temp_density, dtemp, dP, stress, strain_rate,
            stressyy, drho, rho, power, powerTerm, pressureTerm, densityTerm, tdot);

    loop_all_elem(var.egroups, elemf);

    // Combining temperature update and bc in the same loop for efficiency,
    // since only the top boundary has Dirichlet bc, and all the other boundaries
    // have no heat flux bc.
    #pragma omp parallel for default(none)      \
    shared(var, param, tdot, temperature, temp_power, temp_pressure, temp_density, \
          dtemp, powerTerm, pressureTerm, densityTerm, std::cout);

    for (int n=0; n<var.nnode; ++n) {

      double temp_old        = temperature[n];
      double diffusion_term  = -(tdot[n] * var.dt) / (*var.tmass)[n];
      double power_term      = powerTerm[n] / (*var.tmass)[n];
      double pressure_term   = pressureTerm[n] / (*var.tmass)[n];
      double density_term    = densityTerm[n] / (*var.tmass)[n];
      double surface_temp    = param.bc.surface_temperature;

      if ((*var.bcflag)[n] & BOUNDZ1){
        temperature[n]    = surface_temp;
        temp_power[n]     = surface_temp;
        temp_pressure[n]  = surface_temp;
        temp_density[n]   = surface_temp;
        dtemp[n]          = temp_old - temperature[n];
      }else{
        temperature[n]    += (diffusion_term) + (power_term) + (pressure_term) + (density_term);
        temp_power[n]     += power_term;
        temp_pressure[n]  += pressure_term;
        temp_density[n]   += density_term;
        dtemp[n]          = temperature[n] - temp_old;
      }
    }
}


void initial_material_properties(const Variables& var, double_vec& rho) {

    for (int e = 0; e<var.nelem; ++e) {
        rho[e] = var.mat->rho(e);
    }
}

void update_density(const Variables& var, double_vec& rho,
                    double_vec& drho, tensor_t& strain_rate)
{
   for (int e = 0; e<var.nelem; ++e) {
      double rho_old  = 0;
      rho_old = rho[e];
      double *edot = (*var.strain_rate)[e];

#ifdef THREED
      double vedot = edot[0]+ edot[1] + edot[2];
#else
      double vedot = edot[0]+ edot[1];
#endif

      rho[e]  = rho_old * (1 - var.dt * vedot);
      drho[e] = rho[e] - rho_old;
   }
}


void update_strain_rate(const Variables& var, tensor_t& strain_rate)
{
    double *v[NODES_PER_ELEM];

    #pragma omp parallel for default(none) \
        shared(var, strain_rate) private(v)
    for (int e=0; e<var.nelem; ++e) {
        const int *conn = (*var.connectivity)[e];
        const double *shpdx = (*var.shpdx)[e];
        const double *shpdz = (*var.shpdz)[e];
        double *s = (*var.strain_rate)[e];

        for (int i=0; i<NODES_PER_ELEM; ++i)
            v[i] = (*var.vel)[conn[i]];

        // XX component
        int n = 0;
        s[n] = 0;
        for (int i=0; i<NODES_PER_ELEM; ++i)
            s[n] += v[i][0] * shpdx[i];

#ifdef THREED
        const double *shpdy = (*var.shpdy)[e];
        // YY component
        n = 1;
        s[n] = 0;
        for (int i=0; i<NODES_PER_ELEM; ++i)
            s[n] += v[i][1] * shpdy[i];
#endif

        // ZZ component
#ifdef THREED
        n = 2;
#else
        n = 1;
#endif
        s[n] = 0;
        for (int i=0; i<NODES_PER_ELEM; ++i)
            s[n] += v[i][NDIMS-1] * shpdz[i];

#ifdef THREED
        // XY component
        n = 3;
        s[n] = 0;
        for (int i=0; i<NODES_PER_ELEM; ++i)
            s[n] += 0.5 * (v[i][0] * shpdy[i] + v[i][1] * shpdx[i]);
#endif

        // XZ component
#ifdef THREED
        n = 4;
#else
        n = 2;
#endif
        s[n] = 0;
        for (int i=0; i<NODES_PER_ELEM; ++i)
            s[n] += 0.5 * (v[i][0] * shpdz[i] + v[i][NDIMS-1] * shpdx[i]);

#ifdef THREED
        // YZ component
        n = 5;
        s[n] = 0;
        for (int i=0; i<NODES_PER_ELEM; ++i)
            s[n] += 0.5 * (v[i][1] * shpdz[i] + v[i][2] * shpdy[i]);
#endif
    }
}


static void apply_damping(const Param& param, const Variables& var, array_t& force)
{
    // flatten 2d arrays to simplify indexing
    double* ff = force.data();
    const double* v = var.vel->data();
    const double small_vel = 1e-13;
    #pragma omp parallel for default(none)          \
        shared(var, param, ff, v)
    for (int i=0; i<var.nnode*NDIMS; ++i) {
        if (std::fabs(v[i]) > small_vel) {
            ff[i] -= param.control.damping_factor * std::copysign(ff[i], v[i]);
        }
    }
}


void update_force(const Param& param, const Variables& var, array_t& force)
{
    std::fill_n(force.data(), var.nnode*NDIMS, 0);

    class ElemFunc_force : public ElemFunc
    {
    private:
        const Variables &var;
        array_t &force;
        const double gravity;
    public:
        ElemFunc_force(const Variables &var, array_t &force, double gravity) :
            var(var), force(force), gravity(gravity) {};
        void operator()(int e)
        {
            const int *conn = (*var.connectivity)[e];
            const double *shpdx = (*var.shpdx)[e];
#ifdef THREED
            const double *shpdy = (*var.shpdy)[e];
#endif
            const double *shpdz = (*var.shpdz)[e];
            double *s = (*var.stress)[e];
            double vol = (*var.volume)[e];

            double buoy = 0;
            if (gravity != 0)
                buoy = (*var.rho)[e] * gravity / NODES_PER_ELEM;

            for (int i=0; i<NODES_PER_ELEM; ++i) {
                double *f = force[conn[i]];
#ifdef THREED
                f[0] -= (s[0]*shpdx[i] + s[3]*shpdy[i] + s[4]*shpdz[i]) * vol;
                f[1] -= (s[3]*shpdx[i] + s[1]*shpdy[i] + s[5]*shpdz[i]) * vol;
                f[2] -= (s[4]*shpdx[i] + s[5]*shpdy[i] + s[2]*shpdz[i] + buoy) * vol;
#else
                f[0] -= (s[0]*shpdx[i] + s[2]*shpdz[i]) * vol;
                f[1] -= (s[2]*shpdx[i] + s[1]*shpdz[i] + buoy) * vol;
#endif
            }
        }
    } elemf(var, force, param.control.gravity);

    loop_all_elem(var.egroups, elemf);

    apply_stress_bcs(param, var, force);

    if (param.control.is_quasi_static) {
        apply_damping(param, var, force);
    }
}


void update_velocity(const Variables& var, array_t& vel)
{
    const double* m = &(*var.mass)[0];
    // flatten 2d arrays to simplify indexing
    const double* f = var.force->data();
    double* v = vel.data();
    #pragma omp parallel for default(none) \
        shared(var, m, f, v)
    for (int i=0; i<var.nnode*NDIMS; ++i) {
        int n = i / NDIMS;
        v[i] += var.dt * f[i] / m[n];
    }
}


void update_coordinate(const Variables& var, array_t& coord)
{
    double* x = var.coord->data();
    const double* v = var.vel->data();

    #pragma omp parallel for default(none) \
        shared(var, x, v)
    for (int i=0; i<var.nnode*NDIMS; ++i) {
        x[i] += v[i] * var.dt;
    }
}


namespace {

    void jaumann_rate_3d(double *s, double dt, double w3, double w4, double w5)
    {
        double s_inc[NSTR];

        s_inc[0] =-2.0 * s[3] * w3 - 2.0 * s[4] * w4;
        s_inc[1] = 2.0 * s[3] * w3 - 2.0 * s[5] * w5;
        s_inc[2] = 2.0 * s[4] * w4 + 2.0 * s[5] * w5;
        s_inc[3] = s[0] * w3 - s[1] * w3 - s[4] * w5 - s[5] * w4;
        s_inc[4] = s[0] * w4 - s[2] * w4 + s[3] * w5 - s[5] * w3;
        s_inc[5] = s[1] * w5 - s[2] * w5 + s[3] * w4 + s[4] * w3;

        for(int i=0; i<NSTR; ++i)  {
            s[i] += dt * s_inc[i];
        }
    }


    void jaumann_rate_2d(double *s, double dt, double w2)
    {
        double s_inc[NSTR];

        s_inc[0] =-2.0 * s[2] * w2;
        s_inc[1] = 2.0 * s[2] * w2;
        s_inc[2] = s[0] * w2 - s[1] * w2;

        for(int i=0; i<NSTR; ++i)  {
            s[i] += dt * s_inc[i];
        }
    }

}


void rotate_stress(const Variables &var, tensor_t &stress, tensor_t &strain)
{
    // The spin rate tensor, W, and the Cauchy stress tensor, S, are
    //     [  0  w3  w4]     [s0 s3 s4]
    // W = [-w3   0  w5],  S=[s3 s1 s5].
    //     [-w4 -w5   0]     [s4 s5 s2]
    //
    // Stress (and strain) increment based on the Jaumann rate is
    // dt*(S*W-W*S). S*W-W*S is also symmetric.
    // So, following the indexing of stress tensor,
    // we get
    // sj[0] = dt * ( -2 * s3 * w3 - 2 * s4 * w4)
    // sj[1] = dt * (  2 * s3 * w3 - 2 * s5 * w5)
    // sj[2] = dt * (  2 * s4 * w4 + 2 * s5 * w5)
    // sj[3] = dt * ( s0 * w3 - s1 * w3 - s4 * w5 - s5 * w4)
    // sj[4] = dt * ( s0 * w4 - s2 * w4 + s3 * w5 - s5 * w3)
    // sj[5] = dt * ( s1 * w5 - s2 * w5 + s3 * w4 + s4 * w3)

    #pragma omp parallel for default(none) \
        shared(var, stress, strain)
    for (int e=0; e<var.nelem; ++e) {
        const int *conn = (*var.connectivity)[e];

#ifdef THREED

        double w3, w4, w5;
        {
            const double *shpdx = (*var.shpdx)[e];
            const double *shpdy = (*var.shpdy)[e];
            const double *shpdz = (*var.shpdz)[e];

            double *v[NODES_PER_ELEM];
            for (int i=0; i<NODES_PER_ELEM; ++i)
                v[i] = (*var.vel)[conn[i]];

            w3 = 0;
            for (int i=0; i<NODES_PER_ELEM; ++i)
                w3 += 0.5 * (v[i][0] * shpdy[i] - v[i][1] * shpdx[i]);

            w4 = 0;
            for (int i=0; i<NODES_PER_ELEM; ++i)
                w4 += 0.5 * (v[i][0] * shpdz[i] - v[i][NDIMS-1] * shpdx[i]);

            w5 = 0;
            for (int i=0; i<NODES_PER_ELEM; ++i)
                w5 += 0.5 * (v[i][1] * shpdz[i] - v[i][NDIMS-1] * shpdy[i]);
        }

        jaumann_rate_3d(stress[e], var.dt, w3, w4, w5);
        jaumann_rate_3d(strain[e], var.dt, w3, w4, w5);

#else

        double w2;
        {
            const double *shpdx = (*var.shpdx)[e];
            const double *shpdz = (*var.shpdz)[e];

            double *v[NODES_PER_ELEM];
            for (int i=0; i<NODES_PER_ELEM; ++i)
                v[i] = (*var.vel)[conn[i]];

            w2 = 0;
            for (int i=0; i<NODES_PER_ELEM; ++i)
                w2 += 0.5 * (v[i][NDIMS-1] * shpdx[i] - v[i][0] * shpdz[i]);
        }

        jaumann_rate_2d(stress[e], var.dt, w2);
        jaumann_rate_2d(strain[e], var.dt, w2);

#endif
    }
}
