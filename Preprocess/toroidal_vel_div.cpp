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
        const dataset & source_data,
        const std::vector<bool> & mask
    ) {

    size_t index;
    int Itime, Idepth, Ilat, Ilon;
    double dulon_dlon, dulat_dlat, local_lat, ulat, cos_lat, sin_lat, tmp_val;
    std::vector<double*> lon_deriv_vals, lat_deriv_vals;
    std::vector<const std::vector<double>*> lon_deriv_fields, lat_deriv_fields;
    bool is_pole;

    lon_deriv_fields.push_back(&vel_lon);
    lat_deriv_fields.push_back(&vel_lat);

    const std::vector<double>   &latitude   = source_data.latitude,
                                &longitude  = source_data.longitude;
    const std::vector<int>  &myCounts = source_data.myCounts;

    const int   Ntime  = myCounts.at(0),
                Ndepth = myCounts.at(1),
                Nlat   = (constants::GRID_TYPE == constants::GridType::MeshGrid) ? myCounts.at(2) : 1,
                Nlon   = (constants::GRID_TYPE == constants::GridType::MeshGrid) ? myCounts.at(3) : 1;
    const size_t Npts  = vel_lon.size();

    #pragma omp parallel \
    default(none) \
    shared( latitude, longitude, mask, div, vel_lon, vel_lat, \
            lon_deriv_fields, lat_deriv_fields, source_data )\
    private(Itime, Idepth, Ilat, Ilon, index, cos_lat, sin_lat, tmp_val, ulat, \
            local_lat, dulon_dlon, dulat_dlat, lon_deriv_vals, lat_deriv_vals, is_pole ) \
    firstprivate( Nlon, Nlat, Ndepth, Ntime, Npts )
    {

        lon_deriv_vals.push_back(&dulon_dlon);
        lat_deriv_vals.push_back(&dulat_dlat);

        #pragma omp for collapse(1) schedule(guided)
        for (index = 0; index < Npts; ++index) {

            tmp_val = constants::fill_value;

            if (mask.at(index)) { // Skip land areas

                if (constants::GRID_TYPE == constants::GridType::MeshGrid) {
                    Index1to4(index, Itime, Idepth, Ilat, Ilon,
                                     Ntime, Ndepth, Nlat, Nlon);

                    spher_derivative_at_point(
                            lon_deriv_vals, lon_deriv_fields, longitude, "lon",
                            source_data, Itime, Idepth, Ilat, Ilon, mask);

                    spher_derivative_at_point(
                            lat_deriv_vals, lat_deriv_fields, latitude, "lat",
                            source_data, Itime, Idepth, Ilat, Ilon, mask);
                } else {
                    spher_derivative_at_point(
                            lon_deriv_vals, lon_deriv_fields, longitude, "lon",
                            source_data, 0, 0, index, index, mask);

                    spher_derivative_at_point(
                            lat_deriv_vals, lat_deriv_fields, latitude, "lat",
                            source_data, 0, 0, index, index, mask);
                }

                local_lat = (constants::GRID_TYPE == constants::GridType::MeshGrid) 
                                ? latitude.at(Ilat)
                                : latitude.at(index);

                cos_lat = cos(local_lat);
                sin_lat = sin(local_lat);

                ulat = vel_lat.at(index);

                // If we're too close to the pole (less than 0.01 degrees), bad things happen
                is_pole = std::fabs( std::fabs( local_lat * 180.0 / M_PI ) - 90 ) < 0.01;

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

