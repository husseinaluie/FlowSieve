#include "../constants.hpp"
#include "../functions.hpp"
#include "../differentiation_tools.hpp"
#include "../preprocess.hpp"
#include <algorithm>
#include <vector>
#include <omp.h>
#include <math.h>

void toroidal_Lap_F(
        std::vector<double> & out_arr,
        const std::vector<double> & F,
        const std::vector<double> & longitude,
        const std::vector<double> & latitude,
        const int Ntime,
        const int Ndepth,
        const int Nlat,
        const int Nlon,
        const std::vector<bool> & mask
        ) {

    // ret = ddlon(vel_lat) / cos_lat - ddlat( u_lon * cos_lat ) / cos_lat 
    //     = ddlon(vel_lat) / cos_lat - ddlat( u_lon ) + u_lon * tan_lat

    int Ilat, Ilon, index;
    const int Itime  = 0;
    const int Idepth = 0;
    double d1Fdlat1, d2Fdlat2, d2Fdlon2, tmp;
    std::vector<double*> lon2_deriv_vals, lat2_deriv_vals, lat1_deriv_vals;
    std::vector<const std::vector<double>*> deriv_fields;
    bool is_pole;

    deriv_fields.push_back(&F);
    
    #pragma omp parallel \
    default(none) \
    shared( out_arr, latitude, longitude, mask,\
            F, deriv_fields)\
    private(Ilat, Ilon, index, tmp, is_pole, \
            d2Fdlon2, d1Fdlat1, d2Fdlat2, \
            lon2_deriv_vals, lat1_deriv_vals, lat2_deriv_vals)
    {

        lon2_deriv_vals.push_back(&d2Fdlon2);
        lat1_deriv_vals.push_back(&d1Fdlat1);
        lat2_deriv_vals.push_back(&d2Fdlat2);

        #pragma omp for collapse(2) schedule(guided)
        for (Ilat = 0; Ilat < Nlat; ++Ilat) {
            for (Ilon = 0; Ilon < Nlon; ++Ilon) {

                index = Index(0, 0, Ilat, Ilon,
                              1, 1, Nlat, Nlon);
                tmp = constants::fill_value;

                if (mask.at(index)) { // Skip land areas

                    // Second lon derivative
                    spher_derivative_at_point(
                            lon2_deriv_vals, deriv_fields,
                            longitude, "lon",
                            Itime, Idepth, Ilat, Ilon,
                            Ntime, Ndepth, Nlat, Nlon,
                            mask, 2);

                    // First lat derivative
                    spher_derivative_at_point(
                            lat1_deriv_vals, deriv_fields,
                            latitude, "lat",
                            Itime, Idepth, Ilat, Ilon,
                            Ntime, Ndepth, Nlat, Nlon,
                            mask, 1);

                    // Second lat derivative
                    spher_derivative_at_point(
                            lat2_deriv_vals, deriv_fields,
                            latitude, "lat",
                            Itime, Idepth, Ilat, Ilon,
                            Ntime, Ndepth, Nlat, Nlon,
                            mask, 2);

                    // If we're too close to the pole (less than 0.01 degrees), bad things happen
                    is_pole = std::fabs( std::fabs( latitude.at(Ilat) * 180.0 / M_PI ) - 90 ) < 0.01;

                    if ( not(is_pole) ) {
                        tmp =     d2Fdlon2 / pow(cos(latitude.at(Ilat)), 2)
                                - d1Fdlat1 * tan(latitude.at(Ilat));
                                + d2Fdlat2;
                    } else {
                        tmp = 0.;
                    }

                    tmp *= 1. / pow(constants::R_earth, 2.);

                }
                out_arr.at(index) = tmp;

            } // end Ilon
        } // end Ilat
    } // end pragma
}
