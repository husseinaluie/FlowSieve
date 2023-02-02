#include "../constants.hpp"
#include "../functions.hpp"
#include "../differentiation_tools.hpp"
#include <algorithm>
#include <vector>
#include <omp.h>
#include <math.h>

void toroidal_vel_div(  
        std::vector<double> & div,
        const std::vector<double> & vel_lon,
        const std::vector<double> & vel_lat,
        const std::vector<double> & longitude,
        const std::vector<double> & latitude,
        const int Ntime,
        const int Ndepth,
        const int Nlat,
        const int Nlon,
        const std::vector<bool> & mask
    ) {

    size_t index;
    int Itime, Idepth, Ilat, Ilon;
    double dulon_dlon, dulat_dlat, ulat, cos_lat, sin_lat, tmp_val;
    std::vector<double*> lon_deriv_vals, lat_deriv_vals;
    std::vector<const std::vector<double>*> lon_deriv_fields, lat_deriv_fields;
    bool is_pole;

    lon_deriv_fields.push_back(&vel_lon);
    lat_deriv_fields.push_back(&vel_lat);

    const size_t Npts = Ntime * Ndepth * Nlat * Nlon;

    #pragma omp parallel \
    default(none) \
    shared( latitude, longitude, mask, div, vel_lon, vel_lat, \
            lon_deriv_fields, lat_deriv_fields )\
    private(Itime, Idepth, Ilat, Ilon, index, cos_lat, sin_lat, tmp_val, ulat, \
            dulon_dlon, dulat_dlat, lon_deriv_vals, lat_deriv_vals, is_pole ) \
    firstprivate( Nlon, Nlat, Ndepth, Ntime, Npts )
    {

        lon_deriv_vals.push_back(&dulon_dlon);
        lat_deriv_vals.push_back(&dulat_dlat);

        #pragma omp for collapse(1) schedule(guided)
        for (index = 0; index < Npts; ++index) {

            tmp_val = constants::fill_value;

            if (mask.at(index)) { // Skip land areas

                Index1to4(index, Itime, Idepth, Ilat, Ilon,
                                 Ntime, Ndepth, Nlat, Nlon);

                spher_derivative_at_point(
                        lon_deriv_vals, lon_deriv_fields,
                        longitude, "lon",
                        Itime, Idepth, Ilat, Ilon,
                        Ntime, Ndepth, Nlat, Nlon,
                        mask);

                spher_derivative_at_point(
                        lat_deriv_vals, lat_deriv_fields,
                        latitude, "lat",
                        Itime, Idepth, Ilat, Ilon,
                        Ntime, Ndepth, Nlat, Nlon,
                        mask);

                cos_lat = cos(latitude.at(Ilat));
                sin_lat = sin(latitude.at(Ilat));

                ulat = vel_lat.at(index);

                // If we're too close to the pole (less than 0.01 degrees), bad things happen
                is_pole = std::fabs( std::fabs( latitude.at(Ilat) * 180.0 / M_PI ) - 90 ) < 0.01;

                if ( is_pole ) {
                    tmp_val = 0.;
                } else {
                    tmp_val = dulon_dlon  +  cos_lat * dulat_dlat  -  ulat * sin_lat ;
                    tmp_val *= 1. / ( cos_lat * constants::R_earth );
                }

            }
            div.at(index) = tmp_val;
        }
    }
}

