#include "../../constants.hpp"
#include "../../functions.hpp"
#include "../../differentiation_tools.hpp"
#include <algorithm>
#include <vector>
#include <omp.h>
#include <math.h>

void toroidal_curl_u_dot_er(
        std::vector<double> & out_arr,
        const std::vector<double> & u_lon,
        const std::vector<double> & u_lat,
        const std::vector<double> & longitude,
        const std::vector<double> & latitude,
        const int Ntime,
        const int Ndepth,
        const int Nlat,
        const int Nlon,
        const std::vector<double> & mask
        ) {

    // ret = ddlon(vel_lat) / cos_lat - ddlat( u_lon * cos_lat ) / cos_lat 
    //     = ddlon(vel_lat) / cos_lat - ddlat( u_lon ) + u_lon * tan_lat

    int Itime, Idepth, Ilat, Ilon, index;
    double dulat_dlon, dulon_dlat, tmp;
    std::vector<double*> lon_deriv_vals, lat_deriv_vals;
    std::vector<const std::vector<double>*> deriv_fields;

    deriv_fields.push_back(&u_lon);
    deriv_fields.push_back(&u_lat);
    
    #pragma omp parallel \
    default(none) \
    shared( out_arr, latitude, longitude, mask,\
            u_lon, u_lat, deriv_fields)\
    private(Itime, Idepth, Ilat, Ilon, index, tmp, \
            dulat_dlon, dulon_dlat, lon_deriv_vals, lat_deriv_vals)
    {

        lon_deriv_vals.push_back(NULL);
        lon_deriv_vals.push_back(&dulat_dlon);

        lat_deriv_vals.push_back(&dulon_dlat);
        lat_deriv_vals.push_back(NULL);

        #pragma omp for collapse(1) schedule(guided)
        for (index = 0; index < (int)u_lon.size(); index++) {

            tmp = constants::fill_value;

            if (mask.at(index) == 1) { // Skip land areas

                Index1to4(index, Itime, Idepth, Ilat, Ilon,
                                 Ntime, Ndepth, Nlat, Nlon);

                spher_derivative_at_point(
                        lon_deriv_vals, deriv_fields,
                        longitude, "lon",
                        Itime, Idepth, Ilat, Ilon,
                        Ntime, Ndepth, Nlat, Nlon,
                        mask);

                spher_derivative_at_point(
                        lat_deriv_vals, deriv_fields,
                        latitude, "lat",
                        Itime, Idepth, Ilat, Ilon,
                        Ntime, Ndepth, Nlat, Nlon,
                        mask);

                //  = ddlon(vel_lat) / cos_lat - ddlat( u_lon ) + u_lon * tan_lat
                tmp =   dulat_dlon / cos(latitude.at(Ilat))
                      - dulon_dlat 
                      + u_lon.at(index) * tan(latitude.at(Ilat));

                tmp *= 1. / constants::R_earth;

            }
            out_arr.at(index) = tmp;

        } // end index
    } // end pragma
}
