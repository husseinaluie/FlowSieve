#include "../../constants.hpp"
#include "../../functions.hpp"
#include "../../differentiation_tools.hpp"
#include <algorithm>
#include <vector>
#include <omp.h>
#include <math.h>

void curl_vec_dot_er(
        std::vector<double> & out_arr,
        const std::vector<double> & u_lon,
        const std::vector<double> & u_lat,
        const dataset & source_data,
        const std::vector<short int> & mask
        ) {

    // ret = ddlon(vel_lat) / cos_lat - ddlat( u_lon * cos_lat ) / cos_lat 

    const size_t Npts  = u_lon.size();
    size_t index;
    const std::vector<double>   &latitude   = source_data.latitude,
                                &longitude  = source_data.longitude;

    const std::vector<int>  &myCounts = source_data.myCounts;

    const int   Ntime  = myCounts.at(0),
                Ndepth = myCounts.at(1),
                Nlat   = (constants::GRID_TYPE == constants::GridType::MeshGrid) ? myCounts.at(2) : 1,
                Nlon   = (constants::GRID_TYPE == constants::GridType::MeshGrid) ? myCounts.at(3) : 1;

    // Store u_lon * cos_lat in an array
    std::vector<double> cos_lat_u_lon( Npts, 0. );

    #pragma omp parallel \
    default(none) \
    shared( latitude, mask, source_data, u_lon, cos_lat_u_lon )\
    private( index )\
    firstprivate( Nlon, Nlat, Ndepth, Ntime, Npts )
    {
        #pragma omp for collapse(1) schedule(static)
        for (index = 0; index < Npts; ++index) {

            int Itime, Idepth, Ilat, Ilon;

            if (mask.at(index)) { // Skip land areas
                if (constants::GRID_TYPE == constants::GridType::MeshGrid) {
                    Index1to4(index, Itime, Idepth, Ilat, Ilon,
                                     Ntime, Ndepth, Nlat, Nlon);
                }

                // check if we're too close to the pole
                double local_lat = (constants::GRID_TYPE == constants::GridType::MeshGrid) 
                                ? latitude.at(Ilat)
                                : latitude.at(index);

                cos_lat_u_lon[index] = cos(local_lat) * u_lon[index];
            }
        }
    } // end pragma

    std::vector<const std::vector<double>*> deriv_fields = { &cos_lat_u_lon, &u_lat };
    
    #pragma omp parallel \
    default(none) \
    shared( out_arr, latitude, longitude, mask, source_data, deriv_fields )\
    private( index ) \
    firstprivate( Nlon, Nlat, Ndepth, Ntime, Npts )
    {

        double dulat_dlon, duloncoslat_dlat;

        std::vector<double*> lon_deriv_vals = { NULL, &dulat_dlon }, 
                             lat_deriv_vals = { &duloncoslat_dlat, NULL };

        #pragma omp for collapse(1) schedule(guided)
        for (index = 0; index < Npts; ++index) {

            int Itime, Idepth, Ilat, Ilon;

            double tmp = 0.;
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

                // check if we're too close to the pole
                double local_lat = (constants::GRID_TYPE == constants::GridType::MeshGrid) 
                                ? latitude.at(Ilat)
                                : latitude.at(index);

                // ret = ddlon(vel_lat) / cos_lat - ddlat( u_lon * cos_lat ) / cos_lat
                tmp = dulat_dlon - duloncoslat_dlat;
                tmp *= 1. / ( cos(local_lat) * constants::R_earth );
            }

            out_arr.at(index) = tmp;

        }
    } // end pragma
}
