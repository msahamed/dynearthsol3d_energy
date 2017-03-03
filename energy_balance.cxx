void update_stress_energy(const Variables& var, tensor_t& stress, double_vec& stressyy, double_vec& thermal_stress, double_vec& dP,
                   tensor_t& strain, double_vec& plstrain, double_vec& delta_plstrain, double_vec& dtemp,
                   tensor_t& strain_rate, double_vec& power, double_vec& tenergy, double_vec& venergy, double_vec& denergy)
{
    const int rheol_type = var.mat->rheol_type;

    #pragma omp parallel for default(none)                           \
         shared(var, stress, stressyy, power, dP, dtemp, strain, plstrain, delta_plstrain, strain_rate, std::cerr)
    for (int e=0; e<var.nelem; ++e) {

        // stress, strain and strain_rate of this element
        double* s    = stress[e];
        double& syy  = stressyy[e];
        double* es   = strain[e];
        double* edot = strain_rate[e];

        // anti-mesh locking correction on strain rate
        if(1){
            double div = trace(edot);
            //double div2 = ((*var.volume)[e] / (*var.volume_old)[e] - 1) / var.dt;
            for (int i=0; i<NDIMS; ++i) {
                edot[i] += ((*var.edvoldt)[e] - div) / NDIMS;  // XXX: should NDIMS -> 3 in plane strain?
            }
        }


        double dT = 0;
        const int *conn = (*var.connectivity)[e];
        for (int i = 0; i < NODES_PER_ELEM; ++i) {
          dT += dtemp[conn[i]];
        }
        dT /= NODES_PER_ELEM;

        double de[NSTR];
        double alpha = var.mat->alpha(e);
        // double thermal_strain = (alpha * dT)/3;

        // normal components
        for (int i = 0; i<NDIMS; ++i){
           es[i] += (edot[i] * var.dt); // + thermal_strain;
           de[i] =  (edot[i] * var.dt); // + thermal_strain;
        }

        // Shear components
        for (int i = NDIMS; i<NSTR; ++i) {
            es[i] += edot[i] * var.dt;
            de[i] = edot[i] * var.dt;
        }

        switch (rheol_type) {
        case MatProps::rh_elastic:
            {
                double bulkm = var.mat->bulkm(e);
                double shearm = var.mat->shearm(e);
#ifdef THREED
                double pressure_old=-(s[0] + s[1] + s[2]) / NDIMS;
#else
                double pressure_old=-(s[0] + s[1]) / NDIMS;
#endif
                elastic(bulkm, shearm, de, s);
#ifdef THREED
                double pressure_new = (-(s[0] + s[1] + s[2]) / NDIMS);
#else
                double pressure_new=-(s[0] + s[1] + syy) / 3;
#endif
                dP[e]=(pressure_new-pressure_old);
            }
            break;
        case MatProps::rh_viscous:
            {
                double bulkm = var.mat->bulkm(e);
                double viscosity = var.mat->visc(e);
                double total_dv = trace(es);
                viscous(bulkm, viscosity, total_dv, edot, s);
            }
            break;
        case MatProps::rh_maxwell:
            {
                double bulkm = var.mat->bulkm(e);
                double shearm = var.mat->shearm(e);
                double viscosity = var.mat->visc(e);
                double dv = (*var.volume)[e] / (*var.volume_old)[e] - 1;
                maxwell(bulkm, shearm, viscosity, var.dt, dv, de, s);
            }
            break;
        case MatProps::rh_ep:
            {
                double t_power = 0;
                double v_power = 0;
                double d_power = 0;
                double depls = 0;
                double bulkm = var.mat->bulkm(e);
                double shearm = var.mat->shearm(e);
                double amc, anphi, anpsi, hardn, ten_max;
                var.mat->plastic_props(e, plstrain[e],
                                       amc, anphi, anpsi, hardn, ten_max);
                int failure_mode;

#ifdef THREED
                double pressure_old = -(s[0]+s[1]+s[2])/NDIMS;
#else
                double pressure_old = -(s[0]+s[1]+ syy)/3.0; // plane strain
#endif

                double thermal_stress = -bulkm * alpha * dT;;
                thermal_stress[e] += thermal_stress;

                if (var.mat->is_plane_strain) {
                    elasto_plastic2d(bulkm, shearm, t_power,v_power, d_power, amc, anphi, anpsi, hardn, ten_max, de, depls, s, syy, failure_mode, thermal_stress);
                }
                else {
                    elasto_plastic(bulkm, shearm, amc, anphi, anpsi, hardn, ten_max,
                                   de, depls, s, failure_mode);
                }

#ifdef THREED
                double pressure_new = -(s[0]+s[1]+s[2])/NDIMS;
#else
                double pressure_new = -(s[0]+s[1]+syy)/3;
#endif

                dP[e] = (pressure_new-pressure_old);
                plstrain[e] += depls;
                delta_plstrain[e] = depls;
                power[e]  = t_power;

                tenergy[e] += t_power;
                venergy[e] += v_power;
                denergy[e] += d_power;


            }
            break;
        case MatProps::rh_evp:
            {
                double t_power = 0;
                double v_power = 0;
                double d_power = 0;
                double depls        = 0;
                double bulkm        = var.mat->bulkm(e);
                double shearm       = var.mat->shearm(e);
                double viscosity    = var.mat->visc(e);
                double dv           = (*var.volume)[e] / (*var.volume_old)[e] - 1;

                // stress due to maxwell rheology
                double sv[NSTR];
                for (int i=0; i<NSTR; ++i) sv[i] = s[i];
                maxwell(bulkm, shearm, viscosity, var.dt, dv, de, sv);
                double svII = second_invariant2(sv);

                double amc, anphi, anpsi, hardn, ten_max;
                var.mat->plastic_props(e, plstrain[e],
                                       amc, anphi, anpsi, hardn, ten_max);

                // stress due to elasto-plastic rheology
                double sp[NSTR], spyy;
                for (int i=0; i<NSTR; ++i) sp[i] = s[i];
                int failure_mode;

                double dT = 0;
                const int *conn = (*var.connectivity)[e];
                for (int i = 0; i < NODES_PER_ELEM; ++i) {
                     dT += dtemp[conn[i]];
                }
                dT /= NODES_PER_ELEM;


                if (var.mat->is_plane_strain) {
                    spyy = syy;
                    elasto_plastic2d(bulkm, shearm, t_power,v_power, d_power, amc, anphi, anpsi, hardn, ten_max,
                                     de, depls, s, syy, failure_mode, thermal_stress);
                }
                else {
                    elasto_plastic(bulkm, shearm, amc, anphi, anpsi, hardn, ten_max,
                                   de, depls, sp, failure_mode);
                }
                double spII = second_invariant2(sp);

                // use the smaller as the final stress
                if (svII < spII)
                    for (int i=0; i<NSTR; ++i) s[i] = sv[i];
                else {
                    for (int i=0; i<NSTR; ++i) s[i] = sp[i];
                    plstrain[e] += depls;
                    delta_plstrain[e] = depls;
                    syy = spyy;
                }
            }
            break;
        default:
            std::cerr << "Error: unknown rheology type: " << rheol_type << "\n";
            std::exit(1);
            break;
        }
        // std::cerr << "stress " << e << ": ";
        // print(std::cerr, s, NSTR);
        // std::cerr << '\n';
    }
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
