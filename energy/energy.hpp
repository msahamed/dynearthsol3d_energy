#ifndef DYNEARTHSOL3D_ENERGY_HPP
#define DYNEARTHSOL3D_ENERGY_HPP

#include "../parameters.hpp"

class Energy{
    
    public:
        // energy balance equation related variables
        double_vec *dtemp;
        double_vec *dP;
        double_vec *rho;
        double_vec *drho;
        double_vec *power;
        double_vec *tenergy;
        double_vec *venergy;
        double_vec *denergy;
        double_vec *powerTerm;
        double_vec *pressureTerm;
        double_vec *densityTerm;
        double_vec *thermal_energy;
        double_vec *elastic_energy;
        double_vec *temp_power, *temp_pressure, *temp_density;
        double_vec *thermal_stress, *ediffStress, *ndiffStress;
        tensor_t *elastic_strain;
        const Variables  &var;
        const Param &param;

        Energy(const Param &param, Variables& var):param(param), var(var){};
        void allocate_energy_variables(const Param &param, Variables& var);
        void initial_material_properties(Variables& var, double_vec &rho);

        void update_stress(const Variables& var, tensor_t& stress, double_vec& stressyy, double_vec& thermal_stress, double_vec& dP,
            tensor_t& strain, tensor_t& elastic_strain, double_vec& plstrain, double_vec& delta_plstrain, double_vec& dtemp,
            tensor_t& strain_rate, double_vec& power, double_vec& tenergy, double_vec& venergy,double_vec& denergy);
        
        void update_temperature(const Param &param, const Variables &var, double_vec &temperature,
            double_vec &temp_power, double_vec &temp_pressure, double_vec &temp_density,
            double_vec &dtemp,double_vec &dP, double_vec &tdot, tensor_t &stress,
            tensor_t &strain_rate, double_vec &stressyy, double_vec &drho, double_vec &rho,
            double_vec &power, double_vec &powerTerm, double_vec &pressureTerm, double_vec &densityTerm);

        void apply_NMD_to_Stress(const Variables &var, tensor_t& stress, double_vec &stressyy,
            double_vec & ediffStress, double_vec & ndiffStress, double_vec& dP);

        void update_density(const Variables& var, double_vec& rho,
                            double_vec& drho, tensor_t& strain_rate);
        void initial_material_properties(const Variables& var, double_vec& rho);
        void update_thermal_energy(const Variables& var, double_vec& thermal_energy);
        void update_elastic_energy(const Variables& var, double_vec& elastic_energy);


};

#endif