#ifndef DYNEARTHSOL3D_FIELDS_HPP
#define DYNEARTHSOL3D_FIELDS_HPP

void allocate_variables(const Param &param, Variables& var);
void reallocate_variables(const Param &param, Variables& var);
void update_temperature(const Param &param, const Variables &var, double_vec &temperature,
                        double_vec &temp_power, double_vec &temp_pressure, double_vec &temp_density,
                        double_vec &dtemp,double_vec &dP, double_vec &tdot, tensor_t &stress,
                        tensor_t &strain_rate, double_vec &stressyy, double_vec &drho, double_vec &rho,
                        double_vec &power, double_vec &powerTerm, double_vec &pressureTerm, double_vec &densityTerm);

void update_strain_rate(const Variables& var, tensor_t& strain_rate);
void update_force(const Param& param, const Variables& var, array_t& force);
void update_velocity(const Variables& var, array_t& vel);
void update_coordinate(const Variables& var, array_t& coord);
void rotate_stress(const Variables &var, tensor_t &stress, tensor_t &strain);
// void update_density(const Variables& var, double_vec& rho,
//                     double_vec& drho, tensor_t& strain_rate);
// void initial_material_properties(const Variables& var, double_vec& rho);
// void update_thermal_energy(const Variables& var, double_vec& thermal_energy);
// void update_elastic_energy(const Variables& var, double_vec& elastic_energy);
#endif
