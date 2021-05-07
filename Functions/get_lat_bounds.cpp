#include <stdio.h>
#include <math.h>    
#include <algorithm>
#include <vector>
#include "../constants.hpp"
#include "../functions.hpp"

/*!
 * \brief Get latitude bounds for coarse-graining.
 *
 *  See the methods documentation (\ref methods1) for more details.
 *
 *  @param[in,out]  LAT_lb,LAT_ub       logical index bounds for the filtering domain
 *  @param[in]      latitude            1D latitude vector
 *  @param[in]      Ilat                logical index of current position
 *  @param[in]      scale               filtering scale (metres)
 *
 */
void get_lat_bounds(
        int & LAT_lb,
        int & LAT_ub,
        const std::vector<double> & latitude,
        const int Ilat,
        const double scale) {

    const double KernPad = constants::KernPad;
    const double ref_lat = latitude.at(Ilat);
    const int    Nlat    = (int) latitude.size();
    
    if (KernPad < 0) {
        // KernPad < 0 means we use the entire domain, so don't bother
        //   calculating an integration region.
        LAT_lb = 0;
        LAT_ub = Nlat;
    } else {
        if (constants::UNIFORM_LAT_GRID) {
            // On spherical coordinates, the distance calculator is
            //      very expensive. So when the grid is uniform, we
            //      avoid extra calls to distance with the extra
            //      information that we have.

            const double dlat = fabs(latitude.at( 1) - latitude.at( 0));
            double dlat_m;
            int dlat_N;

            // The spacing (in metres and points) betwee latitude gridpoints
            //   The factor of 2 is diameter->radius 
            if (constants::CARTESIAN) { dlat_m = dlat; } 
            else                      { dlat_m = dlat * constants::R_earth; }

            dlat_N = ceil( ( KernPad * scale / dlat_m ) / 2.);
            dlat_N = std::min(Nlat, dlat_N);

            // Handle periodicity, if appropriate
            if (constants::PERIODIC_Y) {
                LAT_lb = Ilat - dlat_N;
                LAT_ub = Ilat + dlat_N;
                if (LAT_ub - LAT_lb > Nlat) { LAT_ub = LAT_lb + Nlat; }
            } else {
                LAT_lb = std::max(0,    Ilat - dlat_N);
                LAT_ub = std::min(Nlat, Ilat + dlat_N);
            }
        } else {
            // When the grid isn't uniform, we need to be a little more
            //      careful to find the optimal integration domain in
            //      logical space. To do this, we'll take advantage of
            //      built-in binary search routines

            // Binary search for the first time comparison returns false
            //      on the interval [0, Ilat]
            //   Embedded lambda operator gives comparison
            int tmp_lb = 
                std::lower_bound(
                        latitude.begin(), 
                        latitude.begin()+Ilat, 
                        scale * KernPad,
                        [=](double lat, double ref) {
                            double d = distance(0, ref_lat, 0, lat);
                            return d > ref;
                        }
                        )
                - latitude.begin();
            if (tmp_lb > 0) { tmp_lb--; }

            // Binary search for the first time comparison returns false
            //      on the interval [Ilat, Nlat]
            //   Embedded lambda operator gives comparison
            int tmp_ub = 
                std::lower_bound(
                        latitude.begin()+Ilat, 
                        latitude.end(), 
                        scale * KernPad,
                        [=](double lat, double ref) {
                            double d = distance(0, ref_lat, 0, lat);
                            return d < ref;
                        }
                    )
                - latitude.begin();

            LAT_lb = tmp_lb;
            LAT_ub = tmp_ub;
        }
    }
}
