#include "../constants.hpp"
#include "../functions.hpp"
#include "../differentiation_tools.hpp"
#include <algorithm>
#include <vector>
#include <omp.h>
#include <math.h>

void compute_div_vel(
        std::vector<double> & div,
        const std::vector<double> & u_x,
        const std::vector<double> & u_y,
        const std::vector<double> & u_z,
        const std::vector<double> & longitude,
        const std::vector<double> & latitude,
        const int Ntime,
        const int Ndepth,
        const int Nlat,
        const int Nlon,
        const std::vector<double> & mask
        ) {

    double div_tmp;

    int Itime, Idepth, Ilat, Ilon, index, mask_index;

    double ux_x, uy_y, uz_z;
    std::vector<double*> x_deriv_vals, y_deriv_vals, z_deriv_vals;
    std::vector<const std::vector<double>*> deriv_fields;

    deriv_fields.push_back(&u_x);
    deriv_fields.push_back(&u_y);
    deriv_fields.push_back(&u_z);
    
    #pragma omp parallel \
    default(none) \
    shared( div, Idepth, Itime, stdout,\
            latitude, longitude, mask,\
            u_x, u_y, u_z, deriv_fields)\
    private(Ilat, Ilon, index, mask_index, \
            ux_x, uy_y, uz_z,\
            x_deriv_vals, y_deriv_vals, z_deriv_vals,\
            div_tmp)
    {

        x_deriv_vals.push_back(&ux_x);
        x_deriv_vals.push_back(NULL);
        x_deriv_vals.push_back(NULL);

        y_deriv_vals.push_back(NULL);
        y_deriv_vals.push_back(&uy_y);
        y_deriv_vals.push_back(NULL);

        z_deriv_vals.push_back(NULL);
        z_deriv_vals.push_back(NULL);
        z_deriv_vals.push_back(&uz_z);

        #pragma omp for collapse(2) schedule(dynamic)
        for (Ilat = 0; Ilat < Nlat; Ilat++) {
            for (Ilon = 0; Ilon < Nlon; Ilon++) {

                mask_index = Index(
                        0,     0,      Ilat, Ilon,
                        Ntime, Ndepth, Nlat, Nlon);

                for (Itime = 0; Itime < Ntime; Itime++) {
                    for (Idepth = 0; Idepth < Ndepth; Idepth++) {

                        index = Index(Itime, Idepth, Ilat, Ilon,
                                      Ntime, Ndepth, Nlat, Nlon);

                        if (mask.at(mask_index) == 1) { // Skip land areas

                            Cart_derivatives_at_point(
                                    x_deriv_vals, y_deriv_vals,
                                    z_deriv_vals, deriv_fields,
                                    latitude, longitude,
                                    Itime, Idepth, Ilat, Ilon,
                                    Ntime, Ndepth, Nlat, Nlon,
                                    mask);

                            // u_i,i
                            div_tmp = ux_x + uy_y + uz_z;

                        } // end if(water) block
                        else { // if(land)
                            div_tmp = constants::fill_value;
                        } // end if(land)

                        div.at(index) = div_tmp;

                    } // end Idepth loop
                } // end Itime loop
            } // end Ilon loop
        } // end Ilat loop
    } // end pragma
} // end function
