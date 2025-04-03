#include "algorithm"
#include "iostream"

#include "ANN/ANN.h"

#include "constants.hpp"
#include "parameters.hpp"

#include "barycentric-fn.hpp"
#include "utils.hpp"
#include "brc-interpolation.hpp"

namespace { // anonymous namespace

typedef Array2D<double,NODES_PER_ELEM> brc_t;


void interpolate_field(const brc_t &brc, const int_vec &el, const conn_t &connectivity,
                       const double_vec &source, double_vec &target)
{
    const size_t target_size = target.size();
    const size_t chunk_size = 64; // Cache-friendly chunk size
    
    #pragma omp parallel for default(none)          \
        shared(brc, el, connectivity, source, target, target_size, chunk_size) \
        schedule(static, chunk_size)
    for (std::size_t i=0; i<target_size; i++) {
        int e = el[i];
        const int *conn = connectivity[e];
        const double *brc_row = brc[i];
        
        // Prefetch next iteration's data
        if (i < target_size - 1) {
            const int next_e = el[i+1];
            const double *next_source = &source[next_e];
            const double *next_brc = brc[i+1];
            __builtin_prefetch(next_source, 0, 3);
            __builtin_prefetch(next_brc, 0, 3);
        }
        
        double result = 0;
        #pragma omp simd reduction(+:result)
        for (int j=0; j<NODES_PER_ELEM; j++) {
            result += source[conn[j]] * brc_row[j];
        }
        target[i] = result;
    }
}


void interpolate_field(const brc_t &brc, const int_vec &el, const conn_t &connectivity,
                       const array_t &source, array_t &target)
{
    const size_t target_size = target.size();
    const size_t chunk_size = 64; // Cache-friendly chunk size
    
    #pragma omp parallel for default(none)          \
        shared(brc, el, connectivity, source, target, target_size, chunk_size) \
        schedule(static, chunk_size)
    for (std::size_t i=0; i<target_size; i++) {
        int e = el[i];
        const int *conn = connectivity[e];
        const double *brc_row = brc[i];
        double *target_row = target[i];
        
        // Prefetch next iteration's data
        if (i < target_size - 1) {
            const int next_e = el[i+1];
            const double *next_source = source[next_e];
            const double *next_brc = brc[i+1];
            __builtin_prefetch(next_source, 0, 3);
            __builtin_prefetch(next_brc, 0, 3);
        }
        
        // Compute interpolation for all dimensions
        #pragma omp simd
        for (int d=0; d<NDIMS; d++) {
            target_row[d] = 0;
            for (int j=0; j<NODES_PER_ELEM; j++) {
                target_row[d] += source[conn[j]][d] * brc_row[j];
            }
        }
    }
}


void prepare_interpolation(const Variables &var,
                           const Barycentric_transformation &bary,
                           const array_t &old_coord,
                           const conn_t &old_connectivity,
                           const std::vector<int_vec> &old_support,
                           brc_t &brc, int_vec &el)
{
    // for each new coord point, find the enclosing old element

    // ANN requires double** as input
    double **points = new double*[old_coord.size()];
    for (std::size_t i=0; i<old_coord.size(); i++) {
        points[i] = const_cast<double*>(old_coord[i]);
    }
    ANNkd_tree kdtree(points, old_coord.size(), NDIMS);

    const int k = 1;
    const double eps = 0;
    int nn_idx[k];
    double dd[k];

    // Note: kdtree.annkSearch() is not thread-safe, cannot use openmp in this loop
    for (int i=0; i<var.nnode; i++) {
        double *q = (*var.coord)[i];

        // find the nearest point nn in old_coord
        kdtree.annkSearch(q, k, nn_idx, dd, eps);
        int nn = nn_idx[0];

        // elements surrounding nn
        const int_vec &nn_elem = old_support[nn];

        // std::cout << i << " ";
        // print(std::cout, q, NDIMS);
        // std::cout << " " << nn << " " << dd[0] << '\n';

        double r[NDIMS];
        int e;

        // shortcut: q is exactly the same as nn
        if (dd[0] == 0) {
            e = nn_elem[0];
            bary.transform(q, e, r);
            // r should be a permutation of [1, 0, 0]
            // normalize r to remove round-off error
            for (int d=0; d<NDIMS; d++) {
                if (r[d] > 0.9)
                    r[d] = 1;
                else
                    r[d] = 0;
            }
            goto found;
        }

        // loop over (old) elements surrounding nn to find
        // the element that is enclosing q
        for (std::size_t j=0; j<nn_elem.size(); j++) {
            e = nn_elem[j];
            bary.transform(q, e, r);
            if (bary.is_inside(r)) {
                // std::cout << e << " ";
                // print(std::cout, r, NDIMS);
                // std::cout << '\n';
                goto found;
            }
        }

        /* not_found */

        {
            /* Situation: q is in the upper element, but its nearest point is o!
             * we won't find the enclosing element with the method above
             *     x
             *    / \   <-- this is a large triangle
             *   / q                            \
             *  x---- x
             *   \-o-/   <-- this is a small triangle
             */

            // this array contains the elements that have been searched so far
            int_vec searched(nn_elem);

            // search through elements that are neighbors of nn_elem
            for (std::size_t j=0; j<nn_elem.size(); j++) {
                int ee = nn_elem[j];
                const int *conn = old_connectivity[ee];
                for (int m=0; m<NODES_PER_ELEM; m++) {
                    // np is a node close to q
                    int np = conn[m];
                    const int_vec &np_elem = old_support[np];
                    for (std::size_t j=0; j<np_elem.size(); j++) {
                        e = np_elem[j];
                        auto it = std::find(searched.begin(), searched.end(), e);
                        if (it != searched.end()) {
                            // this element has been searched before
                            continue;
                        }
                        searched.push_back(e);
                        bary.transform(q, e, r);
                        // std::cout << e << " ";
                        // print(std::cout, r, NDIMS);
                        // std::cout << " ... \n";
                        if (bary.is_inside(r)) {
                            goto found;
                        }
                    }
                }
            }
        }
        {
            //std::cout << "New node is outside of the old domain. \n";

            // Situation: q must be outside the old domain
            // using nearest old_coord instead
            e = nn_elem[0];
            bary.transform(points[nn], e, r);
        }
    found:
        el[i] = e;
        double sum = 0;
        for (int d=0; d<NDIMS; d++) {
            brc[i][d] = r[d];
            sum += r[d];
        }
        brc[i][NODES_PER_ELEM-1] = 1 - sum;
    }

    delete [] points;

    // print(std::cout, *var.coord);
    // std::cout << '\n';
    // print(std::cout, el);
    // std::cout << '\n';
    // print(std::cout, bar);
    // std::cout << '\n';
}

} // anonymous namespace


void barycentric_node_interpolation(Variables &var,
                                    const Barycentric_transformation &bary,
                                    const array_t &old_coord,
                                    const conn_t &old_connectivity)
{
    int_vec el(var.nnode);
    brc_t brc(var.nnode);
    prepare_interpolation(var, bary, old_coord, old_connectivity, *var.support, brc, el);
    double_vec *a;

    a = new double_vec(var.nnode);
    interpolate_field(brc, el, old_connectivity, *var.temperature, *a);
    delete var.temperature;
    var.temperature = a;

    a = new double_vec(var.nnode);
    interpolate_field(brc, el, old_connectivity, *var.temp_power, *a);
    delete var.temp_power;
    var.temp_power = a;

    a = new double_vec(var.nnode);
    interpolate_field(brc, el, old_connectivity, *var.temp_pressure, *a);
    delete var.temp_pressure;
    var.temp_pressure = a;

    a = new double_vec(var.nnode);
    interpolate_field(brc, el, old_connectivity, *var.temp_density, *a);
    delete var.temp_density;
    var.temp_density = a;

    a = new double_vec(var.nnode);
    interpolate_field(brc, el, old_connectivity, *var.dtemp, *a);
    delete var.dtemp;
    var.dtemp = a;

    a = new double_vec(var.nnode);
    interpolate_field(brc, el, old_connectivity, *var.powerTerm, *a);
    delete var.powerTerm;
    var.powerTerm = a;

    a = new double_vec(var.nnode);
    interpolate_field(brc, el, old_connectivity, *var.pressureTerm, *a);
    delete var.pressureTerm;
    var.pressureTerm = a;

    a = new double_vec(var.nnode);
    interpolate_field(brc, el, old_connectivity, *var.densityTerm, *a);
    delete var.densityTerm;
    var.densityTerm = a;

    a = new double_vec(var.nnode);
    interpolate_field(brc, el, old_connectivity, *var.ndiffStress, *a);
    delete var.ndiffStress;
    var.ndiffStress = a;

    array_t *b = new array_t(var.nnode);
    interpolate_field(brc, el, old_connectivity, *var.vel, *b);
    delete var.vel;
    var.vel = b;

    b = new array_t(var.nnode);
    interpolate_field(brc, el, old_connectivity, *var.coord0, *b);
    delete var.coord0;
    var.coord0 = b;

}


void barycentric_node_interpolation_forT(const Variables &var,
                                         const Barycentric_transformation &bary,
                                         const array_t &input_coord,
                                         const conn_t &input_connectivity,
                                         const std::vector<int_vec> &input_support,
					 const double_vec &inputtemperature,
					 double_vec &outputtemperature)
{
    int_vec el(var.nnode);
    brc_t brc(var.nnode);
    prepare_interpolation(var, bary, input_coord, input_connectivity, input_support, brc, el);

    interpolate_field(brc, el, input_connectivity, inputtemperature, outputtemperature);
}
