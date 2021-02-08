#include "../constants.hpp"
#include "../functions.hpp"
#include "../differentiation_tools.hpp"
#include <algorithm>
#include <vector>
#include <omp.h>
#include <math.h>

void compute_div_vel(
        std::vector<double> & div,
        std::vector<double> & OkuboWeiss,
        const std::vector<double> & u_x,
        const std::vector<double> & u_y,
        const std::vector<double> & u_z,
        const std::vector<double> & longitude,
        const std::vector<double> & latitude,
        const int Ntime,
        const int Ndepth,
        const int Nlat,
        const int Nlon,
        const std::vector<bool> & mask
        ) {

    const int OMP_chunksize = get_omp_chunksize(Nlat,Nlon);

    double div_tmp, OkuboWeiss_tmp;

    int Itime, Idepth, Ilat, Ilon;
    size_t index;
    const size_t Npts = u_x.size();

    double ux_x, ux_y, ux_z,
           uy_x, uy_y, uy_z,
           uz_x, uz_y, uz_z;
    std::vector<double*> x_deriv_vals, y_deriv_vals, z_deriv_vals;
    std::vector<const std::vector<double>*> deriv_fields;

    deriv_fields.push_back(&u_x);
    deriv_fields.push_back(&u_y);
    deriv_fields.push_back(&u_z);
    
    #pragma omp parallel \
    default(none) \
    shared( div, OkuboWeiss, latitude, longitude, mask,\
            u_x, u_y, u_z, deriv_fields)\
    private(Itime, Idepth, Ilat, Ilon, index,\
            ux_x, ux_y, ux_z, uy_x, uy_y, uy_z, uz_x, uz_y, uz_z, \
            x_deriv_vals, y_deriv_vals, z_deriv_vals,\
            div_tmp, OkuboWeiss_tmp)
    {

        x_deriv_vals.push_back(&ux_x);
        x_deriv_vals.push_back(&uy_x);
        x_deriv_vals.push_back(&uz_x);

        y_deriv_vals.push_back(&ux_y);
        y_deriv_vals.push_back(&uy_y);
        y_deriv_vals.push_back(&uz_y);

        z_deriv_vals.push_back(&ux_z);
        z_deriv_vals.push_back(&uy_z);
        z_deriv_vals.push_back(&uz_z);

        #pragma omp for collapse(1) schedule(guided, OMP_chunksize)
        for (index = 0; index < Npts; index++) {

            div_tmp = constants::fill_value;
            OkuboWeiss_tmp = constants::fill_value;

            if ( mask.at(index) ) { // Skip land areas

                Index1to4(index, Itime, Idepth, Ilat, Ilon,
                                 Ntime, Ndepth, Nlat, Nlon);

                Cart_derivatives_at_point(
                        x_deriv_vals, y_deriv_vals,
                        z_deriv_vals, deriv_fields,
                        latitude, longitude,
                        Itime, Idepth, Ilat, Ilon,
                        Ntime, Ndepth, Nlat, Nlon,
                        mask);

                // u_i,i
                div_tmp = ux_x + uy_y + uz_z;

                // 2D Okubo-Weiss
                // (u_1,1 + u_2,2)^2 + 4 u_1,2 u_2,1
                OkuboWeiss_tmp = pow(ux_x - uy_y, 2) + 4 * ux_y * uy_x;

            }

            div.at(index) = div_tmp;
            OkuboWeiss.at(index) = OkuboWeiss_tmp;

        } // end index
    } // end pragma
} // end function
