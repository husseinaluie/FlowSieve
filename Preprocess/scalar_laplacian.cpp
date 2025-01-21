#include "../constants.hpp"
#include "../functions.hpp"
#include "../differentiation_tools.hpp"
#include <algorithm>
#include <vector>
#include <omp.h>
#include <math.h>

void scalar_laplacian(  
        std::vector<double> & Lap_scalar,
        const std::vector<double> & scalar,
        const dataset & source_data,
        const std::vector<short int> & mask
    ) {

    size_t index;
    std::vector<const std::vector<double>*> deriv_field;
    deriv_field.push_back(&scalar);

    const std::vector<double>   &latitude   = source_data.latitude,
                                &longitude  = source_data.longitude;
    const std::vector<int>  &myCounts = source_data.myCounts;

    const int   Ntime  = myCounts.at(0),
                Ndepth = myCounts.at(1),
                Nlat   = (constants::GRID_TYPE == constants::GridType::MeshGrid) ? myCounts.at(2) : 1,
                Nlon   = (constants::GRID_TYPE == constants::GridType::MeshGrid) ? myCounts.at(3) : 1;
    const size_t Npts  = scalar.size();

    #pragma omp parallel default(none) \
    shared( latitude, longitude, mask, Lap_scalar, scalar, deriv_field, source_data )\
    private( index ) firstprivate( Nlon, Nlat, Ndepth, Ntime, Npts )
    {
        double d2dlon2, d2dlat2, ddlat;
        std::vector<double*> lon2_deriv_vals, lat2_deriv_vals, lat1_deriv_vals;

        lon2_deriv_vals.push_back(&d2dlon2);
        lat2_deriv_vals.push_back(&d2dlat2);
        lat1_deriv_vals.push_back(&ddlat);

        #pragma omp for collapse(1) schedule(guided)
        for (index = 0; index < Npts; ++index) {

            double tmp_val = 0;
            if (mask.at(index)) { // Skip land areas

                int Ilat;
                if (constants::GRID_TYPE == constants::GridType::MeshGrid) {
                    int Itime, Idepth, Ilon;
                    Index1to4(index, Itime, Idepth, Ilat, Ilon,
                                     Ntime, Ndepth, Nlat, Nlon);

                    spher_derivative_at_point(
                            lon2_deriv_vals, deriv_field, longitude, "lon",
                            source_data, Itime, Idepth, Ilat, Ilon, mask, 2);

                    spher_derivative_at_point(
                            lat2_deriv_vals, deriv_field, latitude, "lat",
                            source_data, Itime, Idepth, Ilat, Ilon, mask, 2);

                    spher_derivative_at_point(
                            lat1_deriv_vals, deriv_field, latitude, "lat",
                            source_data, Itime, Idepth, Ilat, Ilon, mask, 1);
                } else {
                    spher_derivative_at_point(
                            lon2_deriv_vals, deriv_field, longitude, "lon",
                            source_data, 0, 0, index, index, mask, 2);

                    spher_derivative_at_point(
                            lat2_deriv_vals, deriv_field, latitude, "lat",
                            source_data, 0, 0, index, index, mask, 2);

                    spher_derivative_at_point(
                            lat1_deriv_vals, deriv_field, latitude, "lat",
                            source_data, 0, 0, index, index, mask, 1);
                }

                double local_lat = (constants::GRID_TYPE == constants::GridType::MeshGrid) 
                                ? latitude.at(Ilat)
                                : latitude.at(index);

                // If we're too close to the pole bad things happen
                bool is_pole = std::fabs( std::fabs( local_lat * 180.0 / M_PI ) - 90 ) < 1e-6;

                if ( is_pole ) {
                    tmp_val = 0.;
                } else {
                    double cos_lat = cos(local_lat);
                    double tan_lat = tan(local_lat);

                    tmp_val = d2dlon2 / pow( cos_lat, 2 ) + d2dlat2  - tan_lat * ddlat;
                    tmp_val *= 1. / pow(constants::R_earth, 2);
                }

            }
            Lap_scalar.at(index) = tmp_val;
        }
    }
}

