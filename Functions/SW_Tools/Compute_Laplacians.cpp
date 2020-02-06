#include <math.h>
#include <algorithm>
#include <vector>
#include <omp.h>
#include "../../differentiation_tools.hpp"
#include "../../functions.hpp"
#include "../../functions_sw.hpp"

void Compute_Laplacians(
        std::vector<double> & h_lap_u, 
        std::vector<double> & h_lap_v, 
        std::vector<double> &   lap_h, 
        const std::vector<double> & full_u, 
        const std::vector<double> & full_v, 
        const std::vector<double> & full_h, 
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
    deriv_fields.push_back(&full_u);
    deriv_fields.push_back(&full_v);
    deriv_fields.push_back(&full_h);

    std::vector<double*> x_deriv_vals, y_deriv_vals, z_deriv_vals;
    double dudxx, dudyy, dvdxx, dvdyy, dhdxx, dhdyy;
    
    //
    #pragma omp parallel \
    default(none) \
    shared(full_u, full_v, full_h, deriv_fields, \
            h_lap_u, h_lap_v, lap_h, mask, longitude, latitude,\
            Ntime, Ndepth, Nlat, Nlon) \
    private(Itime, Idepth, Ilat, Ilon, index,\
            x_deriv_vals, y_deriv_vals, z_deriv_vals,\
            dudxx, dudyy, dvdxx, dvdyy, dhdxx, dhdyy)
    {

        x_deriv_vals.clear();
        x_deriv_vals.push_back(&dudxx);
        x_deriv_vals.push_back(&dvdxx);
        x_deriv_vals.push_back(&dhdxx);

        y_deriv_vals.clear();
        y_deriv_vals.push_back(&dudyy);
        y_deriv_vals.push_back(&dvdyy);
        y_deriv_vals.push_back(&dhdyy);

        z_deriv_vals.clear();
        z_deriv_vals.push_back(NULL);
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
                    mask, 2);

            h_lap_u.at(index) = full_h.at(index) * (dudxx + dudyy);
            h_lap_v.at(index) = full_h.at(index) * (dvdxx + dvdyy);
              lap_h.at(index) =                    (dhdxx + dhdyy);

        }
    }
}
