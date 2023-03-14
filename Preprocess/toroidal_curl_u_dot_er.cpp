#include "../constants.hpp"
#include "../functions.hpp"
#include "../differentiation_tools.hpp"
#include "../preprocess.hpp"
#include <algorithm>
#include <vector>
#include <omp.h>
#include <math.h>

void toroidal_curl_u_dot_er(
        std::vector<double> & out_arr,
        const std::vector<double> & u_lon,
        const std::vector<double> & u_lat,
        const dataset & source_data,
        const std::vector<bool>   & mask,
        const std::vector<double> * seed
        ) {

    // ret = ddlon(vel_lat) / cos_lat - ddlat( u_lon * cos_lat ) / cos_lat 
    //     = ddlon(vel_lat) / cos_lat - ddlat( u_lon ) + u_lon * tan_lat

    size_t index, index_sub;
    int Itime, Idepth, Ilat, Ilon;
    double dulat_dlon, dulon_dlat, tmp, local_lat;
    std::vector<double*> lon_deriv_vals, lat_deriv_vals;
    std::vector<const std::vector<double>*> deriv_fields;
    bool is_pole;

    deriv_fields.push_back(&u_lon);
    deriv_fields.push_back(&u_lat);

    const std::vector<double>   &latitude   = source_data.latitude,
                                &longitude  = source_data.longitude;
    const std::vector<int>  &myCounts = source_data.myCounts;

    const int   Ntime  = myCounts.at(0),
                Ndepth = myCounts.at(1),
                Nlat   = (constants::GRID_TYPE == constants::GridType::MeshGrid) ? myCounts.at(2) : 1,
                Nlon   = (constants::GRID_TYPE == constants::GridType::MeshGrid) ? myCounts.at(3) : 1;
    const size_t Npts  = u_lon.size();
    
    #pragma omp parallel \
    default(none) \
    shared( out_arr, latitude, longitude, mask, source_data, \
            u_lon, u_lat, deriv_fields, seed )\
    private(Itime, Idepth, Ilat, Ilon, index, index_sub, tmp, is_pole, \
            local_lat, dulat_dlon, dulon_dlat, lon_deriv_vals, lat_deriv_vals ) \
    firstprivate( Nlon, Nlat, Ndepth, Ntime )
    {

        lon_deriv_vals.push_back(NULL);
        lon_deriv_vals.push_back(&dulat_dlon);

        lat_deriv_vals.push_back(&dulon_dlat);
        lat_deriv_vals.push_back(NULL);

        #pragma omp for collapse(1) schedule(guided)
        for (index = 0; index < Npts; ++index) {

            if (mask.at(index)) { // Skip land areas
                if (constants::GRID_TYPE == constants::GridType::MeshGrid) {
                    Index1to4(index, Itime, Idepth, Ilat, Ilon,
                                     Ntime, Ndepth, Nlat, Nlon);

                    spher_derivative_at_point(
                            lon_deriv_vals, deriv_fields, longitude, "lon",
                            source_data, Itime, Idepth, Ilat, Ilon, mask);

                    spher_derivative_at_point(
                            lat_deriv_vals, deriv_fields, latitude, "lat",
                            source_data, Itime, Idepth, Ilat, Ilon, mask);
                } else {
                    spher_derivative_at_point(
                            lon_deriv_vals, deriv_fields, longitude, "lon",
                            source_data, 0, 0, index, index, mask);

                    spher_derivative_at_point(
                            lat_deriv_vals, deriv_fields, latitude, "lat",
                            source_data, 0, 0, index, index, mask);
                }

                // If we're too close to the pole (less than 0.01 degrees), bad things happen
                local_lat = (constants::GRID_TYPE == constants::GridType::MeshGrid) 
                                ? latitude.at(Ilat)
                                : latitude.at(index);
                is_pole = std::fabs( std::fabs( local_lat * 180.0 / M_PI ) - 90 ) < 0.01;

                if (is_pole) {
                    tmp = 0.;
                } else {
                    //  = ddlon(vel_lat) / cos_lat - ddlat( u_lon ) + u_lon * tan_lat
                    tmp =   dulat_dlon / cos(local_lat)
                          - dulon_dlat 
                          + u_lon.at(index) * tan(local_lat);
                    tmp *= 1. / constants::R_earth;
                }

                // If we have a seed, then subtract it off now
                if ( not(seed == NULL) ) { tmp = tmp - seed->at(index); }

            }

            out_arr.at(index) = tmp;

        }
    } // end pragma
}
