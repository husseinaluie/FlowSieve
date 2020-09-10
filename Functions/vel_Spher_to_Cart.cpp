#include <math.h>
#include "../functions.hpp"
#include "../constants.hpp"

void vel_Spher_to_Cart(
            std::vector<double> & u_x,
            std::vector<double> & u_y,
            std::vector<double> & u_z,
            const std::vector<double> & u_r,
            const std::vector<double> & u_lon,
            const std::vector<double> & u_lat,
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
        u_x = u_lon;
        u_y = u_lat;
        u_z = u_r;
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

                    vel_Spher_to_Cart_at_point(     
                            u_x.at(index), u_y.at(  index), u_z.at(  index),
                            u_r.at(index), u_lon.at(index), u_lat.at(index),
                            longitude.at(Ilon), latitude.at(Ilat));
                }
            }
        }
    }

}
