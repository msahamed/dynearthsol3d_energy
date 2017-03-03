#ifndef DYNEARTHSOL3D_IC_READ_TEMP_HPP
#define DYNEARTHSOL3D_IC_READ_TEMP_HPP

void read_external_temperature_from_comsol(const Param &param,
                                           const Variables &var,
                                           double_vec &temperature);

#endif
