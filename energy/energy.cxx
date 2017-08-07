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
