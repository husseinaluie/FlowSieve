#include "../constants.hpp"
#include "../functions.hpp"
#include "../differentiation_tools.hpp"
#include <algorithm>
#include <vector>
#include <omp.h>
#include <math.h>

void compute_div_vel(
        std::vector<double> & div,
        const std::vector<double> * full_div,
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

    int Ilat, Ilon, index, mask_index;

    double ux_x, uy_y, uz_z;
    
    for (int Itime = 0; Itime < Ntime; Itime++) {
        for (int Idepth = 0; Idepth < Ndepth; Idepth++) {
            #pragma omp parallel \
            default(none) \
            shared( div, Idepth, Itime, stdout,\
                    latitude, longitude, mask,\
                    u_x, u_y, u_z, full_div)\
            private(Ilat, Ilon, index, mask_index, \
                    ux_x, uy_y, uz_z,\
                    div_tmp)
            {
                #pragma omp for collapse(2) schedule(dynamic)
                for (Ilat = 0; Ilat < Nlat; Ilat++) {
                    for (Ilon = 0; Ilon < Nlon; Ilon++) {

                        div_tmp = 0.;

                        // Convert our four-index to a one-index
                        index = Index(
                                Itime, Idepth, Ilat, Ilon,
                                Ntime, Ndepth, Nlat, Nlon);
                        mask_index = Index(
                                0,     0,      Ilat, Ilon,
                                Ntime, Ndepth, Nlat, Nlon);

                        if (mask.at(mask_index) == 1) { // Skip land areas

                            ux_x = Cart_derivative_at_point(u_x, latitude, longitude, "x",
                                    Itime, Idepth, Ilat, Ilon,
                                    Ntime, Ndepth, Nlat, Nlon,
                                    mask);
                            uy_y = Cart_derivative_at_point(u_y, latitude, longitude, "y",
                                    Itime, Idepth, Ilat, Ilon,
                                    Ntime, Ndepth, Nlat, Nlon,
                                    mask);
                            uz_z = Cart_derivative_at_point(u_z, latitude, longitude, "z",
                                    Itime, Idepth, Ilat, Ilon,
                                    Ntime, Ndepth, Nlat, Nlon,
                                    mask);

                            // u_i,i
                            div_tmp = ux_x + uy_y + uz_z;
                        }

                        div.at(index) = div_tmp;
                        if (full_div != NULL) {
                            div.at(index) = full_div->at(index) - div.at(index);
                        }

                    } // end Ilon loop
                } // end Ilat loop
            } // end pragma
        } // end Idepth loop
    } // end Itime loop
}
