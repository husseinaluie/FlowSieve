#include <math.h>
#include <algorithm>
#include <vector>
#include <omp.h>
#include "../../constants.hpp"
#include "../../functions.hpp"
#include "../../particles.hpp"

void particles_initial_positions(
        std::vector<double> & starting_lat, 
        std::vector<double> & starting_lon, 
        const int Npts,
        const std::vector<double> & latitude, 
        const std::vector<double> & longitude, 
        const std::vector<double> & mask
        ) {

    const double lon_rng =        longitude.back() - longitude.front(),
                 lat_rng = 0.9 * (latitude.back()  - latitude.front());

    const double lon_mid = 0.5 * ( longitude.back() + longitude.front() ),
                 lat_mid = 0.5 * ( latitude.back()  + latitude.front()  );

    int left, right, bottom, top,
        BL_ind, BR_ind, TL_ind, TR_ind;

    const int Nlat = latitude.size(),
              Nlon = longitude.size();

    double part_lon, part_lat;

    for ( int II = 0; II < Npts; ++II ) {

        part_lon = ( ((double) rand() / (RAND_MAX)) - 0.5) * lon_rng + lon_mid;
        part_lat = ( ((double) rand() / (RAND_MAX)) - 0.5) * lat_rng + lat_mid;

        // Check if this particle is on land
        particles_get_edges(left, right, bottom, top, 
                part_lat, part_lon, latitude, longitude);

        BL_ind = Index(0, 0, bottom, left,
                       1, 1, Nlat,   Nlon);
        BR_ind = Index(0, 0, bottom, right,
                       1, 1, Nlat,   Nlon);
        TL_ind = Index(0, 0, top,    left,
                       1, 1, Nlat,   Nlon);
        TR_ind = Index(0, 0, top,    right,
                       1, 1, Nlat,   Nlon);

        // So long as we are on land, keep picking a new point
        while (   (   mask.at(BL_ind) 
                    * mask.at(TL_ind) 
                    * mask.at(BR_ind) 
                    * mask.at(TR_ind) == 0 
                  ) or (
                      std::isnan(bottom)
                  ) or ( 
                      std::isnan(top) 
                  )) {

            part_lon = ( ((double) rand() / (RAND_MAX)) - 0.5) * lon_rng + lon_mid;
            part_lat = ( ((double) rand() / (RAND_MAX)) - 0.5) * lat_rng + lat_mid;

            // Check if this particle is on land
            particles_get_edges(left, right, bottom, top, 
                    part_lat, part_lon, latitude, longitude);

            BL_ind = Index(0, 0, bottom, left,
                           1, 1, Nlat,   Nlon);
            BR_ind = Index(0, 0, bottom, right,
                           1, 1, Nlat,   Nlon);
            TL_ind = Index(0, 0, top,    left,
                           1, 1, Nlat,   Nlon);
            TR_ind = Index(0, 0, top,    right,
                           1, 1, Nlat,   Nlon);
        }

        // Store the position
        starting_lat.at(II) = part_lat;
        starting_lon.at(II) = part_lon;

    }
}
