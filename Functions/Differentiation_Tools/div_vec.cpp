#include "../../constants.hpp"
#include "../../functions.hpp"
#include "../../differentiation_tools.hpp"
#include <algorithm>
#include <vector>
#include <omp.h>
#include <math.h>

void div_vec(
        std::vector<double> & out_arr,
        const std::vector<double> & u_lon,
        const std::vector<double> & u_lat,
        const dataset & source_data,
        const std::vector<short int> & mask
        ) {

    // ret = ddlon(vel_lon) / cos_lat + ddlat( u_lat * cos_lat ) / cos_lat 

    const size_t Npts  = u_lon.size();

    size_t index;
    const std::vector<double>   &latitude   = source_data.latitude,
                                &longitude  = source_data.longitude;
    const std::vector<int>  &myCounts = source_data.myCounts;

    const int   Ntime  = myCounts.at(0),
                Ndepth = myCounts.at(1),
                Nlat   = (constants::GRID_TYPE == constants::GridType::MeshGrid) ? myCounts.at(2) : 1,
                Nlon   = (constants::GRID_TYPE == constants::GridType::MeshGrid) ? myCounts.at(3) : 1;

    // Store u_lat * cos_lat in an array
    std::vector<double> cos_lat_u_lat( Npts, 0. );

    #pragma omp parallel \
    default(none) \
    shared( latitude, source_data, u_lat, cos_lat_u_lat )\
    private( index )\
    firstprivate( Nlon, Nlat, Ndepth, Ntime, Npts )
    {
        #pragma omp for collapse(1) schedule(guided)
        for (index = 0; index < Npts; ++index) {

            int Itime, Idepth, Ilat, Ilon;

            if (constants::GRID_TYPE == constants::GridType::MeshGrid) {
                source_data.index1to4_local(index, Itime, Idepth, Ilat, Ilon);
            }

            // check if we're too close to the pole
            double local_lat = (constants::GRID_TYPE == constants::GridType::MeshGrid) 
                            ? latitude.at(Ilat)
                            : latitude.at(index);

            cos_lat_u_lat[index] = cos(local_lat) * u_lat[index];
        }
    } // end pragma

    std::vector<const std::vector<double>*> deriv_fields = { &u_lon, &cos_lat_u_lat };
    
    #pragma omp parallel \
    default(none) \
    shared( out_arr, latitude, longitude, mask, source_data, deriv_fields )\
    private( index ) \
    firstprivate( Nlon, Nlat, Ndepth, Ntime, Npts )
    {

        double dulon_dlon, dulatcoslat_dlat;

        std::vector<double*> lon_deriv_vals = { &dulon_dlon, NULL }, 
                             lat_deriv_vals = { NULL, &dulatcoslat_dlat };

        #pragma omp for collapse(1) schedule(guided)
        for (index = 0; index < Npts; ++index) {

            int Itime, Idepth, Ilat, Ilon;

            if (constants::GRID_TYPE == constants::GridType::MeshGrid) {
                source_data.index1to4_local(index, Itime, Idepth, Ilat, Ilon);

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
            double tmp = dulon_dlon + dulatcoslat_dlat;
            if (constants::GRID_TYPE == constants::GridType::MeshGrid) {
                tmp *= 1. / ( cos(local_lat) * constants::R_earth );
            } else {
                // non-mesh grids have the R_earth baked in
                tmp *= 1. / cos(local_lat);
            }

            out_arr.at(index) = tmp;

        }
    } // end pragma
}
