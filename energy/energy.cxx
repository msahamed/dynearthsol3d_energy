#include "energy.hpp"
#include "../parameters.hpp"
#include "../matprops.hpp"

void Energy::allocate_energy_variables(const Param &param, Variables& var){
    const int n = var.nnode;
    const int e = var.nelem;

    // Energy balance equation related variables :
    dtemp = new double_vec(n); // Temperature difference
    dP    = new double_vec(e); // Pressure difference
    rho   = new double_vec(e); // density
    drho  = new double_vec(e); // density difference
    power    = new double_vec(e);
    tenergy    = new double_vec(e);
    venergy    = new double_vec(e);
    denergy    = new double_vec(e);
    powerTerm  = new double_vec(n);
    pressureTerm     = new double_vec(n);
    densityTerm      = new double_vec(n);
    thermal_energy = new double_vec(e);
    elastic_energy = new double_vec(e);
}

void Energy::initial_material_properties(const Variables& var, double_vec& rho) {
    for (int e = 0; e<var.nelem; ++e) {
        rho[e] = var.mat->rho(e);
    }
}

void Energy::update_density(const Variables& var, double_vec& rho,
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


void Energy::update_thermal_energy(const Variables& var, 
                            double_vec& thermal_energy){

  // #pragma omp parallel for default(none)      \
  // shared(var, connectivity, temperature. volume, thermal_energy, std::cout);
  for (int e = 0; e<var.nelem; ++e) {
      const int *conn = (*var.connectivity)[e];
      double temp =0.0;
      for (int i = 0; i < NODES_PER_ELEM; ++i) {
          temp += (*var.temperature)[conn[i]];
      }
      double T = temp / NODES_PER_ELEM;
      double vol = (*var.volume)[e];
      thermal_energy[e] = (*this->rho)[e] * var.mat->cp(e) * T * vol;
  }
}

void Energy::update_elastic_energy(const Variables& var, double_vec& elastic_energy){
    // #pragma omp parallel for default(none)      \
    // shared(var, stress, elastic_strain. volume, elastic_energy, std::cout);  
  for (int e = 0; e<var.nelem; ++e) {
      double* s = (*var.stress)[e];
      double* es = (*var.elastic_strain)[e];
      double energy = 0.0;
      for (int el = 0; el< NSTR; el++) energy += s[el] * es[el];
      elastic_energy[e] = energy * (*var.volume)[e];
  }
}

void update_energy_balance(const Param &param, const Variables &var, double_vec &temperature,
                        double_vec &temp_power, double_vec &temp_pressure, double_vec &temp_density,
                        double_vec &dtemp,double_vec &dP, double_vec &tdot, tensor_t &stress,
                        tensor_t &strain_rate, double_vec &stressyy, double_vec &drho, double_vec &rho,
                        double_vec &power, double_vec &powerTerm, double_vec &pressureTerm, 
                        double_vec &densityTerm)
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
                             var(var), temperature(temperature), temp_power(temp_power), temp_pressure(temp_pressure), temp_density(temp_density), dtemp(dtemp), dP(dP), stress(stress),strain_rate(strain_rate), stressyy(stressyy), drho(drho), rho(rho), power(power),pressureTerm(pressureTerm), powerTerm(powerTerm), densityTerm(densityTerm),tdot(tdot) {};
        void operator()(int e)
        {
            // diffusion matrix
            double D[NODES_PER_ELEM][NODES_PER_ELEM];

            const int *conn = (*var.connectivity)[e];
            double kv = var.mat->k(e) *  (*var.volume)[e]; // thermal conductivity * volume
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