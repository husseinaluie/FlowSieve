#include "../constants.hpp"
#include "../functions.hpp"
#include "../differentiation_tools.hpp"
#include <algorithm>
#include <vector>
#include <omp.h>
#include <math.h>

void toroidal_vel_from_F(  
        std::vector<double> & vel_lon,
        std::vector<double> & vel_lat,
        const std::vector<double> & F,
        const std::vector<double> & longitude,
        const std::vector<double> & latitude,
        const int Ntime,
        const int Ndepth,
        const int Nlat,
        const int Nlon,
        const std::vector<double> & mask
    ) {

    int Itime, Idepth, Ilat, Ilon, index;
    double dFdlon, dFdlat, cos_lat, tmp_lon, tmp_lat;
    std::vector<double*> lon_deriv_vals, lat_deriv_vals;
    std::vector<const std::vector<double>*> deriv_fields;

    deriv_fields.push_back(&F);

    #pragma omp parallel \
    default(none) \
    shared( latitude, longitude, mask, F, vel_lon, vel_lat, deriv_fields)\
    private(Itime, Idepth, Ilat, Ilon, index, cos_lat, tmp_lon, tmp_lat, \
            dFdlon, dFdlat, lon_deriv_vals, lat_deriv_vals)
    {

        lon_deriv_vals.push_back(&dFdlon);
        lat_deriv_vals.push_back(&dFdlat);

        #pragma omp for collapse(1) schedule(guided)
        for (index = 0; index < (int)F.size(); ++index) {

            tmp_lon = constants::fill_value;
            tmp_lat = constants::fill_value;

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

                cos_lat = cos(latitude.at(Ilat));

                tmp_lon = - dFdlat /  constants::R_earth;
                tmp_lat =   dFdlon / (constants::R_earth * cos_lat);

            }
            vel_lon.at(index) = tmp_lon;
            vel_lat.at(index) = tmp_lat;
        }
    }
}
