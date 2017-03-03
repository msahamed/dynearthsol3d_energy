#ifndef DYNEARTHSOL3D_BC_HPP
#define DYNEARTHSOL3D_BC_HPP

bool is_on_boundary(const Variables &var, int node);
double find_max_vbc(const BC &bc);
void create_boundary_normals(const Variables &var, double bnormals[nbdrytypes][NDIMS],
                             std::map<std::pair<int,int>, double*>  &edge_vectors);
void apply_vbcs(const Param &param, const Variables &var, array_t &vel);
void apply_stress_bcs(const Param& param, const Variables& var, array_t& force);
void surface_processes(const Param& param, const Variables& var, array_t& coord);

#endif
