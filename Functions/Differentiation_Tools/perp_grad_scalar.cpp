#include "../../constants.hpp"
#include "../../functions.hpp"
#include "../../differentiation_tools.hpp"
#include <algorithm>
#include <vector>
#include <omp.h>
#include <math.h>

void perp_grad_scalar(  
        std::vector<double> & vel_lon,
        std::vector<double> & vel_lat,
        const std::vector<double> & F,
        const dataset & source_data,
        const std::vector<short int> & mask
    ) {

    size_t index;
    std::vector<const std::vector<double>*> deriv_fields = { &F };

    const std::vector<double>   &latitude   = source_data.latitude,
                                &longitude  = source_data.longitude;
    const std::vector<int>  &myCounts = source_data.myCounts;

    const int   Ntime  = myCounts.at(0),
                Ndepth = myCounts.at(1),
                Nlat   = (constants::GRID_TYPE == constants::GridType::MeshGrid) ? myCounts.at(2) : 1,
                Nlon   = (constants::GRID_TYPE == constants::GridType::MeshGrid) ? myCounts.at(3) : 1;
    const size_t Npts  = (constants::GRID_TYPE == constants::GridType::MeshGrid) ? Nlat*Nlon : latitude.size();

    #pragma omp parallel \
    default(none) \
    shared( latitude, longitude, mask, vel_lon, vel_lat, deriv_fields, source_data )\
    private( index ) \
    firstprivate( Nlon, Nlat, Ndepth, Ntime, Npts )
    {

        double dFdlon = 0, dFdlat = 0;
        std::vector<double*> lon_deriv_vals = { &dFdlon }, 
                             lat_deriv_vals = { &dFdlat };

        #pragma omp for collapse(1) schedule(guided)
        for (index = 0; index < Npts; ++index) {

            double tmp_lon = 0.;
            double tmp_lat = 0.;

            int Itime, Idepth, Ilat, Ilon;
            if (constants::GRID_TYPE == constants::GridType::MeshGrid) {
                source_data.index1to4_local(index, Itime, Idepth, Ilat, Ilon);

                spher_derivative_at_point(
                        lon_deriv_vals, deriv_fields, longitude, "lon", source_data,
                        Itime, Idepth, Ilat, Ilon, mask);

                spher_derivative_at_point(
                        lat_deriv_vals, deriv_fields, latitude, "lat", source_data,
                        Itime, Idepth, Ilat, Ilon, mask);
            } else {
                spher_derivative_at_point(
                        lon_deriv_vals, deriv_fields, longitude, "lon", source_data,
                        0, 0, index, index, mask);

                spher_derivative_at_point(
                        lat_deriv_vals, deriv_fields, latitude, "lat", source_data,
                        0, 0, index, index, mask);
            }

            if (constants::CARTESIAN) {
                tmp_lon = - dFdlat;
                tmp_lat =   dFdlon;
            } else {
                double local_lat = (constants::GRID_TYPE == constants::GridType::MeshGrid) 
                                ? latitude.at(Ilat)
                                : latitude.at(index);
                bool is_pole = std::fabs( std::fabs( local_lat * 180.0 / M_PI ) - 90 ) < 1e-10;
                double cos_lat = cos(local_lat);

                if (constants::GRID_TYPE == constants::GridType::MeshGrid) {
                    tmp_lon = is_pole ? 0. : - dFdlat /  constants::R_earth;
                    tmp_lat = is_pole ? 0. :   dFdlon / (constants::R_earth * cos_lat);
                } else {
                    // non-mesh grids have the R_earth baked in
                    tmp_lon = is_pole ? 0. : - dFdlat;
                    tmp_lat = is_pole ? 0. :   dFdlon / cos_lat;
                }
            }

            vel_lon.at(index) = tmp_lon;
            vel_lat.at(index) = tmp_lat;
        }
    }
}
