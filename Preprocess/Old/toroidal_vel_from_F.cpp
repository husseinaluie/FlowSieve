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
        const dataset & source_data,
        const std::vector<short int> & mask
    ) {

    int Itime, Idepth, Ilat, Ilon;
    size_t index;
    double dFdlon, dFdlat, local_lat, cos_lat, tmp_lon, tmp_lat;
    std::vector<double*> lon_deriv_vals, lat_deriv_vals;
    std::vector<const std::vector<double>*> deriv_fields;
    bool is_pole;

    const std::vector<double>   &latitude   = source_data.latitude,
                                &longitude  = source_data.longitude;
    const std::vector<int>  &myCounts = source_data.myCounts;

    const int   Ntime  = myCounts.at(0),
                Ndepth = myCounts.at(1),
                Nlat   = (constants::GRID_TYPE == constants::GridType::MeshGrid) ? source_data.Nlat : 1,
                Nlon   = (constants::GRID_TYPE == constants::GridType::MeshGrid) ? source_data.Nlon : 1;
    const size_t Npts  = mask.size();
    //const size_t Npts  = (constants::GRID_TYPE == constants::GridType::MeshGrid) ? Nlat*Nlon : latitude.size();

    deriv_fields.push_back(&F);
    fprintf( stdout, "%d, %d, %d, %d\n", Ntime, Ndepth, Nlat, Nlon );
    fprintf( stdout, "%'zu, %'zu, %'zu, %'zu\n", vel_lon.size(), vel_lat.size(), F.size(), mask.size() );

    #pragma omp parallel \
    default(none) \
    shared( latitude, longitude, mask, F, vel_lon, vel_lat, deriv_fields, source_data)\
    private(Itime, Idepth, Ilat, Ilon, index, local_lat, cos_lat, tmp_lon, tmp_lat, \
            dFdlon, dFdlat, lon_deriv_vals, lat_deriv_vals, is_pole) \
    firstprivate( Nlon, Nlat, Ndepth, Ntime, Npts )
    {

        lon_deriv_vals.push_back(&dFdlon);
        lat_deriv_vals.push_back(&dFdlat);

        #pragma omp for collapse(1) schedule(static)
        for (index = 0; index < Npts; ++index) {

            tmp_lon = 0.;
            tmp_lat = 0.;

            if (mask.at(index)) { // Skip land areas

                if (constants::GRID_TYPE == constants::GridType::MeshGrid) {
                    Index1to4(index, Itime, Idepth, Ilat, Ilon,
                                     Ntime, Ndepth, Nlat, Nlon);

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
                    // If we're too close to the pole bad things happen
                    local_lat = (constants::GRID_TYPE == constants::GridType::MeshGrid) 
                                    ? latitude.at(Ilat)
                                    : latitude.at(index);
                    is_pole = std::fabs( std::fabs( local_lat * 180.0 / M_PI ) - 90 ) < 1e-4;
                    cos_lat = cos(local_lat);

                    tmp_lon = is_pole ? 0. : - dFdlat /  constants::R_earth;
                    tmp_lat = is_pole ? 0. :   dFdlon / (constants::R_earth * cos_lat);
                }

            }
            vel_lon.at(index) = tmp_lon;
            vel_lat.at(index) = tmp_lat;
        }
    }
}
