#include <math.h>
#include "../functions.hpp"
#include "../constants.hpp"

/*!
 * \brief Wrapper that applies vel_Spher_to_Cart_at_point to every point in the domain.
 *
 * @param[in,out]   u_r,u_lon,u_lat     Computed Spherical velocities
 * @param[in]       u_x,u_y,u_z         Cartesian velocities to convert
 * @param[in]       source_data         dataset class instance containing data (Psi, Phi, etc)
 *
 */
void vel_Cart_to_Spher(
            std::vector<double> & u_r,
            std::vector<double> & u_lon,
            std::vector<double> & u_lat,
            const std::vector<double> & u_x,
            const std::vector<double> & u_y,
            const std::vector<double> & u_z,
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
        u_lon = u_x;
        u_lat = u_y;
        u_r   = u_z;
    } else {
        #pragma omp parallel default(none) \
        private( Itime, Idepth, Ilat, Ilon, index ) \
        shared( u_x, u_y, u_z, u_r, u_lon, u_lat, longitude, latitude, mask, source_data ) \
        firstprivate( Nlon, Nlat, Ndepth, Ntime )
        {
            #pragma omp for collapse(1) schedule(guided)
            for (index = 0; index < u_lon.size(); ++index) {

                if ( mask.at(index) ) { // Skip land areas
                    Index1to4( index, Itime, Idepth, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon );

                    vel_Cart_to_Spher_at_point(     
                            u_r.at(index), u_lon.at(index), u_lat.at(index),
                            u_x.at(index), u_y.at(  index), u_z.at(  index),
                            longitude.at(Ilon), latitude.at(Ilat) );
                }
            }
        }
    }
}
