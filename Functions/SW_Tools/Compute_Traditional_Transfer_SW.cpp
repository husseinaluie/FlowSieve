#include <math.h>
#include <algorithm>
#include <vector>
#include <omp.h>
#include "../../differentiation_tools.hpp"
#include "../../functions.hpp"
#include "../../functions_sw.hpp"

void Compute_Traditional_Transfer_SW(
        std::vector<double> & out_array, 
        const std::vector<double> & u_bar, 
        const std::vector<double> & v_bar, 
        const std::vector<double> & PE,
        const std::vector<double> & time,
        const std::vector<double> & depth,
        const std::vector<double> & latitude,
        const std::vector<double> & longitude,
        const std::vector<double> & mask,
        const int & Ntime,
        const int & Ndepth,
        const int & Nlat,
        const int & Nlon
        ){

    const int num_pts = Ntime * Ndepth * Nlat * Nlon;

    int Itime, Idepth, Ilat, Ilon, index;

    //
    std::vector<const std::vector<double>*> deriv_fields;
    deriv_fields.push_back(&u_bar);
    deriv_fields.push_back(&v_bar);

    std::vector<double*> x_deriv_vals, y_deriv_vals, z_deriv_vals;
    double dudx, dvdy;
    
    //
    #pragma omp parallel \
    default(none) \
    shared(u_bar, v_bar, PE, deriv_fields, out_array, \
            mask, longitude, latitude, Ntime, Ndepth, Nlat, Nlon ) \
    private(Itime, Idepth, Ilat, Ilon, index,\
            x_deriv_vals, y_deriv_vals, z_deriv_vals,\
            dudx, dvdy)
    {

        x_deriv_vals.clear();
        x_deriv_vals.push_back(&dudx);
        x_deriv_vals.push_back(NULL);

        y_deriv_vals.clear();
        y_deriv_vals.push_back(NULL);
        y_deriv_vals.push_back(&dvdy);

        z_deriv_vals.clear();
        z_deriv_vals.push_back(NULL);
        z_deriv_vals.push_back(NULL);

        #pragma omp for collapse(1) schedule(static)
        for (index = 0; index < num_pts; index++) {

            Index1to4(index, Itime, Idepth, Ilat, Ilon,
                             Ntime, Ndepth, Nlat, Nlon);

            // We'll need some derivatives
            Cart_derivatives_at_point(
                    x_deriv_vals, y_deriv_vals,
                    z_deriv_vals, deriv_fields,
                    latitude, longitude,
                    Itime, Idepth, Ilat, Ilon,
                    Ntime, Ndepth, Nlat, Nlon,
                    mask);

            out_array.at(index) = (dudx + dvdy) * PE.at(index);

        }
    }
}

