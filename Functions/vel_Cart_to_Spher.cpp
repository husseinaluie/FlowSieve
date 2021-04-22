#include <math.h>
#include "../functions.hpp"
#include "../constants.hpp"

/*!
 * \brief Wrapper that applies vel_Spher_to_Cart_at_point to every point in the domain.
 *
 * @param[in,out]   u_r,u_lon,u_lat                     Computed Spherical velocities
 * @param[in]       u_x,u_y,u_z                         Cartesian velocities to convert
 * @param[in]       mask                                differentiate land from water 
 * @param[in]       time, depth, latitude, longitude    grid vectors (1D)
 *
 */
void vel_Cart_to_Spher(
            std::vector<double> & u_r,
            std::vector<double> & u_lon,
            std::vector<double> & u_lat,
            const std::vector<double> & u_x,
            const std::vector<double> & u_y,
            const std::vector<double> & u_z,
            const std::vector<bool> & mask,
            const std::vector<double> & time,
            const std::vector<double> & depth,
            const std::vector<double> & latitude,
            const std::vector<double> & longitude
        ) {

    size_t index;
    int Itime, Idepth, Ilat, Ilon;

    const int Ntime  = time.size();
    const int Ndepth = depth.size();
    const int Nlat   = latitude.size();
    const int Nlon   = longitude.size();

    const int OMP_chunksize = get_omp_chunksize(Nlat,Nlon);

    if (constants::CARTESIAN) {
        u_lon = u_x;
        u_lat = u_y;
        u_r   = u_z;
    } else {
        #pragma omp parallel default(none) \
        private(Itime, Idepth, Ilat, Ilon, index) \
        shared(u_x, u_y, u_z, u_r, u_lon, u_lat, longitude, latitude, mask)
        {
            #pragma omp for collapse(1) schedule(guided, OMP_chunksize)
            for (index = 0; index < u_lon.size(); ++index) {

                if ( mask.at(index) ) { // Skip land areas
                    Index1to4( index, Itime, Idepth, Ilat, Ilon,
                                      Ntime, Ndepth, Nlat, Nlon );

                    vel_Cart_to_Spher_at_point(     
                            u_r.at(index), u_lon.at(index), u_lat.at(index),
                            u_x.at(index), u_y.at(  index), u_z.at(  index),
                            longitude.at(Ilon), latitude.at(Ilat) );
                }
            }
        }
    }

}
