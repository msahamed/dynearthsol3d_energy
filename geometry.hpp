#ifndef DYNEARTHSOL3D_GEOMETRY_HPP
#define DYNEARTHSOL3D_GEOMETRY_HPP

double dist2(const double* a, const double* b);
void compute_volume(const array_t &coord, const conn_t &connectivity,
                    double_vec &volume);

void compute_dvoldt(const Variables &var, double_vec &dvoldt);

void compute_edvoldt(const Variables &var, double_vec &dvoldt,
                     double_vec &edvoldt);
void apply_NMD_to_Stress(const Variables &var, tensor_t& stress, double_vec &stressyy,
                        double_vec & ediffStress, double_vec & ndiffStress, double_vec& dP);

double compute_dt(const Param& param, const Variables& var);

void compute_mass(const Param &param,
                  const int_vec &egroups, const conn_t &connectivity,
                  const double_vec &volume, const MatProps &mat,
                  double max_vbc_val, double_vec &volume_n, double_vec& stressyy,
                  double_vec &mass, double_vec &tmass, const Variables& var);

void compute_shape_fn(const array_t &coord, const conn_t &connectivity,
                      const double_vec &volume,
                      const int_vec &egroups,
                      shapefn &shpdx, shapefn &shpdy, shapefn &shpdz);

double elem_quality(const array_t &coord, const conn_t &connectivity,
                    const double_vec &volume, int e);

double worst_elem_quality(const array_t &coord, const conn_t &connectivity,
                          const double_vec &volume, int &worst_elem);

#endif
