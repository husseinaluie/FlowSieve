#include <math.h>
#include "../functions.hpp"
#include "../constants.hpp"

/*!
 * \brief Wrapper that applies vel_Spher_to_Cart_at_point to every point in the domain.
 *
 * @param[in,out]   u_x,u_y,u_z         Computed Cartesian velocities
 * @param[in]       u_r,u_lon,u_lat     Spherical velocities to convert
 * @param[in]       source_data         dataset class instance containing data (Psi, Phi, etc)
 *
 */
void vel_Spher_to_Cart(
            std::vector<double> & u_x,
            std::vector<double> & u_y,
            std::vector<double> & u_z,
            const std::vector<double> & u_r,
            const std::vector<double> & u_lon,
            const std::vector<double> & u_lat,
            const dataset & source_data
        ) {

    const std::vector<double>   &latitude   = source_data.latitude,
                                &longitude  = source_data.longitude;

    const std::vector<bool> &mask = source_data.mask;

    const int   Ntime   = source_data.Ntime,
                Ndepth  = source_data.Ndepth,
                Nlat    = source_data.Nlat,
                Nlon    = source_data.Nlon;

    size_t index;
    int Itime, Idepth, Ilat, Ilon;

    if (constants::CARTESIAN) {
        u_x = u_lon;
        u_y = u_lat;
        u_z = u_r;
    } else {
        #pragma omp parallel default(none) \
        private( Itime, Idepth, Ilat, Ilon, index ) \
        shared( u_x, u_y, u_z, u_r, u_lon, u_lat, longitude, latitude, mask, source_data )
        {
            #pragma omp for collapse(1) schedule(guided)
            for (index = 0; index < u_lon.size(); ++index) {

                if ( mask.at(index) ) { // Skip land areas
                    Index1to4( index, Itime, Idepth, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon );

                    vel_Spher_to_Cart_at_point(     
                            u_x.at(index), u_y.at(  index), u_z.at(  index),
                            u_r.at(index), u_lon.at(index), u_lat.at(index),
                            longitude.at(Ilon), latitude.at(Ilat));
                }
            }
        }
    }
}
