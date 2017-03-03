#include <iostream>

#include "constants.hpp"
#include "parameters.hpp"
#include "matprops.hpp"

#include "ic-read-temp.hpp"
#include "ic.hpp"

namespace {

    class Zone
    {
    public:
        virtual ~Zone() {};
        virtual bool contains(const double x[NDIMS]) const = 0;
    };

    class Empty_zone : public Zone
    {
    public:
        bool contains(const double x[NDIMS]) const {return false;}
    };


    class Planar_zone : public Zone
    {
    private:
        const double az, incl;
        const double halfwidth; // in meter
#ifdef THREED
        const double ymin, ymax; // in meter
#endif
        const double zmin, zmax; // in meter
        const double *x0;

    public:
        Planar_zone(const double center[NDIMS], double azimuth, double inclination, double halfwidth_,
#ifdef THREED
                    double ymin_, double ymax_,
#endif
                    double zmin_, double zmax_) :
            az(std::tan(azimuth * DEG2RAD)), incl(1/std::tan(inclination * DEG2RAD)), halfwidth(halfwidth_),
#ifdef THREED
            ymin(ymin_), ymax(ymax_),
#endif
            zmin(zmin_), zmax(zmax_),
            x0(center) // Copy the pointer only, not the data. The caller needs to keep center alive.
        {}

        bool contains(const double x[NDIMS]) const
        {
            // Is x within halfwidth distance to a plane cutting through x0?
            return (x[NDIMS-1] > zmin &&
                    x[NDIMS-1] < zmax &&
#ifdef THREED
                    x[1] > ymin &&
                    x[1] < ymax &&
#endif
                    std::fabs( (x[0] - x0[0])
#ifdef THREED
                               - az * (x[1] - x0[1])
#endif
                               + incl * (x[NDIMS-1] - x0[NDIMS-1]) ) < halfwidth );
        }
    };


    class Ellipsoidal_zone : public Zone
    {
    private:
        const double *x0;
        double semi_axis2[NDIMS];

    public:
        Ellipsoidal_zone(const double center[NDIMS], const double semi_axis[NDIMS]) :
            x0(center) // Copy the pointer only, not the data. The caller needs to keep center alive.
        {
            for(int i=0; i<NDIMS; i++)
                semi_axis2[i] =  semi_axis[i] * semi_axis[i];
        }

        bool contains(const double x[NDIMS]) const
        {
            return ( (x[0] - x0[0])*(x[0] - x0[0])/semi_axis2[0]
#ifdef THREED
                     + (x[1] - x0[1])*(x[1] - x0[1])/semi_axis2[1]
#endif
                     + (x[NDIMS-1] - x0[NDIMS-1])*(x[NDIMS-1] - x0[NDIMS-1])/semi_axis2[NDIMS-1] < 1 );
        }
    };

} // anonymous namespace

double get_maximum_depth(const Variables& var, double zlength){
  int nnodes = (var.bnodes[5]).size();
  int index = (var.bnodes[5])[0];
  double first_node = (*var.coord)[index][0];
  int depth = std::abs(zlength) + first_node;

  for (int i = 0; i < nnodes; i++){
    index = (var.bnodes[5])[i];
    if (depth > std::abs(zlength) + (*var.coord)[index][NDIMS-1]){
      depth = depth;
    }else{
      depth = std::abs(zlength) + (*var.coord)[index][NDIMS-1];
    }
  }
}

double get_actual_depth(const Variables& var, double zcenter, double xcenter){

  int nnodes = (var.bnodes[5]).size();
  int index = (var.bnodes[5])[0];
  double first_node = (*var.coord)[index][0];
  int diff = abs(xcenter - first_node);
  int xvalue = 0;
  int zvalue = 0;

  for (int i = 0; i < nnodes; i++){
    index = (var.bnodes[5])[i];
    if (diff > abs(xcenter - (*var.coord)[index][0])){
        diff = abs(xcenter - (*var.coord)[index][0]);
        xvalue = (*var.coord)[index][0];
        zvalue = (*var.coord)[index][NDIMS-1];
    }
  }

  return zvalue + std::abs(zcenter);
}

void initial_stress_state(const Param &param, const Variables &var,
                          tensor_t &stress, double_vec &stressyy, double_vec &dP, tensor_t &strain,
                          double &compensation_pressure)
{
    if (param.control.gravity == 0) {
        compensation_pressure = 0;
        return;
    }

    // lithostatic condition for stress and strain
    double rho = var.mat->rho(0);
    double ks = var.mat->bulkm(0);

    for (int e=0; e<var.nelem; ++e) {
        const int *conn = (*var.connectivity)[e];
        double zcenter = 0;
        double xcenter = 0;
        for (int i=0; i<NODES_PER_ELEM; ++i) {
            xcenter += (*var.coord)[conn[i]][0];
            zcenter += (*var.coord)[conn[i]][NDIMS-1];
        }
        xcenter /= NODES_PER_ELEM;
        zcenter /= NODES_PER_ELEM;

        double p;
        if (param.ic.is_topo_considered){
          double depth = get_actual_depth(var, zcenter, xcenter);
          p = var.mat->rho(e)* param.control.gravity * depth;
        }else{
          p = ref_pressure(param, zcenter);
        }


        if (param.control.ref_pressure_option == 1 ||
            param.control.ref_pressure_option == 2) {
            ks = var.mat->bulkm(e);
        }

        for (int i=0; i<NDIMS; ++i) {
            stress[e][i] = -p;
            strain[e][i] = -p / ks / NDIMS;
        }
        if (param.mat.is_plane_strain)
            stressyy[e] = -p;
    }

    if (param.ic.is_topo_considered){
      double depth = get_maximum_depth(var, param.mesh.zlength);
      compensation_pressure = ref_pressure(param, depth);
    }else{
      compensation_pressure = ref_pressure(param, -param.mesh.zlength);
    }

}


void initial_weak_zone(const Param &param, const Variables &var,
                       double_vec &plstrain)
{
    Zone *weakzone;

    // TODO: adding different types of weak zone
    double plane_center[NDIMS]; // this variable must outlive weakzone
    switch (param.ic.weakzone_option) {
    case 0:
        weakzone = new Empty_zone();
        break;
    case 1:
        // a planar weak zone, cut through top center
        plane_center[0] = param.ic.weakzone_xcenter * param.mesh.xlength;
#ifdef THREED
        plane_center[1] = param.ic.weakzone_ycenter * param.mesh.ylength;
#endif
        plane_center[NDIMS-1] = -param.ic.weakzone_zcenter * param.mesh.zlength;
        weakzone = new Planar_zone(plane_center,
                                   param.ic.weakzone_azimuth,
                                   param.ic.weakzone_inclination,
                                   param.ic.weakzone_halfwidth * param.mesh.resolution,
#ifdef THREED
                                   param.ic.weakzone_y_min * param.mesh.ylength,
                                   param.ic.weakzone_y_max * param.mesh.ylength,
#endif
                                   -param.ic.weakzone_depth_max * param.mesh.zlength,
                                   -param.ic.weakzone_depth_min * param.mesh.zlength);
        break;
    case 2:
        // a ellipsoidal weak zone
        double semi_axis[NDIMS];
        plane_center[0] = param.ic.weakzone_xcenter * param.mesh.xlength;
        semi_axis[0] = param.ic.weakzone_xsemi_axis;
#ifdef THREED
        plane_center[1] = param.ic.weakzone_ycenter * param.mesh.ylength;
        semi_axis[1] = param.ic.weakzone_ysemi_axis;
#endif
        plane_center[NDIMS-1] = -param.ic.weakzone_zcenter * param.mesh.zlength;
        semi_axis[NDIMS-1] = param.ic.weakzone_zsemi_axis;
        weakzone = new Ellipsoidal_zone(plane_center, semi_axis);
        break;
    default:
        std::cerr << "Error: unknown weakzone_option: " << param.ic.weakzone_option << '\n';
        std::exit(1);
    }

    for (int e=0; e<var.nelem; ++e) {
        const int *conn = (*var.connectivity)[e];
        // the coordinate of the center of this element
        double center[NDIMS] = {0};
        for (int i=0; i<NODES_PER_ELEM; ++i) {
            for (int d=0; d<NDIMS; ++d) {
                center[d] += (*var.coord)[conn[i]][d];
            }
        }
        for (int d=0; d<NDIMS; ++d) {
            center[d] /= NODES_PER_ELEM;
        }

        if (weakzone->contains(center))
            plstrain[e] = param.ic.weakzone_plstrain;

        // Find the most abundant marker mattype in this element
        // int_vec &a = (*var.elemmarkers)[e];
        // int material = std::distance(a.begin(), std::max_element(a.begin(), a.end()));
    }

    delete weakzone;
}


void initial_temperature(const Param &param, const Variables &var,
                         double_vec &temperature)
{
    switch(param.ic.temperature_option) {
    case 0:
        {
            const double age = param.ic.oceanic_plate_age_in_yr * YEAR2SEC;
            const MatProps &mat = *var.mat;
            const double diffusivity = mat.k(0) / mat.rho(0) / mat.cp(0); // thermal diffusivity of 0th element
            double depth;

            for (int i=0; i<var.nnode; ++i) {

              if (param.ic.is_topo_considered){
                depth = get_actual_depth(var, (*var.coord)[i][NDIMS-1], (*var.coord)[i][0]);
              }else{
                depth = (*var.coord)[i][NDIMS-1];
              }
                double w = - depth/ std::sqrt(4 * diffusivity * age);
                temperature[i] = param.bc.surface_temperature +
                    (param.bc.mantle_temperature - param.bc.surface_temperature) * std::erf(w);
            }
            break;
        }
    case 90:
        read_external_temperature_from_comsol(param, var, *var.temperature);
        break;
    default:
        std::cout << "Error: unknown ic.temperature option: " << param.ic.temperature_option << '\n';
        std::exit(1);
    }
}
