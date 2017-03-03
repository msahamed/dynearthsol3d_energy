#ifndef DYNEARTHSOL3D_RHEOLOGY_HPP
#define DYNEARTHSOL3D_RHEOLOGY_HPP

void update_stress(const Variables& var, tensor_t& stress, double_vec& stressyy, double_vec& thermal_stress, double_vec& dP,
                   tensor_t& strain, double_vec& plstrain, double_vec& delta_plstrain, double_vec& dtemp,
                   tensor_t& strain_rate, double_vec& power, double_vec& tenergy, double_vec& venergy,double_vec& denergy);

#endif
