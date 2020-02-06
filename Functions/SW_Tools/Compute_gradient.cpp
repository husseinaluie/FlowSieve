#include <math.h>
#include <algorithm>
#include <vector>
#include <omp.h>
#include "../../differentiation_tools.hpp"
#include "../../functions.hpp"
#include "../../functions_sw.hpp"

void Compute_gradient(
        std::vector<double> & field_x, 
        std::vector<double> & field_y, 
        const std::vector<double> & field,
        const std::vector<double> & latitude,
        const std::vector<double> & longitude,
        const int & Ntime,
        const int & Ndepth,
        const int & Nlat,
        const int & Nlon,
        const std::vector<double> & mask
        ){

    const int num_pts = Ntime * Ndepth * Nlat * Nlon;

    int Itime, Idepth, Ilat, Ilon, index;

    //
    std::vector<const std::vector<double>*> deriv_fields;
    deriv_fields.push_back(&field);

    std::vector<double*> x_deriv_vals, y_deriv_vals, z_deriv_vals;
    double dfdx, dfdy;
    
    //
    #pragma omp parallel \
    default(none) \
    shared(field, deriv_fields, field_x, field_y, \
            mask, longitude, latitude, Ntime, Ndepth, Nlat, Nlon ) \
    private(Itime, Idepth, Ilat, Ilon, index,\
            x_deriv_vals, y_deriv_vals, z_deriv_vals,\
            dfdx, dfdy)
    {

        x_deriv_vals.clear();
        x_deriv_vals.push_back(&dfdx);

        y_deriv_vals.clear();
        y_deriv_vals.push_back(&dfdy);

        z_deriv_vals.clear();
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

            // store values
            field_x.at(index) = dfdx;
            field_y.at(index) = dfdy;

        }
    }
}

