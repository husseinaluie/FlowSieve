#include <stdio.h>
#include <math.h>    
#include <algorithm>
#include <vector>
#include <cassert>
#include "../constants.hpp"
#include "../functions.hpp"

/*!
 * \brief Get longitude integratation bounds for a specific latitude
 *
 *  See the methods documentation (\ref methods1) for more details.
 *
 *  @param[in,out]  LON_lb,LON_ub       logical index bounds for the filtering domain
 *  @param[in]      longitude           1D longitude vector
 *  @param[in]      Ilon                logical index of current position
 *  @param[in]      centre_lat          latitude for the centre of the filtering kernel
 *  @param[in]      curr_lat            current latitude
 *  @param[in]      scale               filtering scale (metres)
 *
 */
void get_lon_bounds(
        int & LON_lb,
        int & LON_ub,
        const std::vector<double> & longitude,
        const int Ilon,
        const double centre_lat,
        const double curr_lat,
        const double scale) {

    const double dlon    = longitude.at( 1) - longitude.at( 0);
    const double KernPad = constants::KernPad;
    const int    Nlon    = (int) longitude.size();

    // assumes uniform lon grid
    static_assert( constants::UNIFORM_LON_GRID, "Currently required uniform lon grid." );

    int dlon_N;
    
    if (KernPad < 0) {
        // KernPad < 0 means we use the entire domain, so don't bother
        //   calculating an integration region.
        LON_lb = 0;
        LON_ub = Nlon;
    } else {
        // This is how big of a 'circle' we actually want (i.e. padded out to capture kernel tapering)
        const double padded_scale = KernPad * scale;

        // The spacing (in metres) between longitude gridpoints at current latitude
        double dlon_m;
        if (constants::CARTESIAN) { dlon_m = dlon; } 
        else                      { dlon_m = dlon * constants::R_earth * cos(curr_lat); }

        // Next determine how far we need to go (since we're already 
        //   at a finite distance in latitude, we don't need to go
        //   the full scale distance in longitude).
        // Essentially, use 'circular' integration regions, not square
        //    this will reduce the number of cells
        //    required, which should improve performance.
        double local_scale;
        if (constants::CARTESIAN) { 
            // Simply use Cartesian pythagorean
            const double delta_lat = fabs( centre_lat - curr_lat ); 
            const double square_diff = pow(padded_scale, 2.) - pow(delta_lat, 2.);
            local_scale = sqrt( fmax( square_diff, 0. ) );

            dlon_N = int( std::min( (double) Nlon, ceil( ( local_scale / dlon_m ) / 2.) ) ); 
        }
        else { 

            // lambda0, phi0 are coordinates for centre of kernel "centre"
            // phi1 is latitude of current integration area "curr"

            
            // This is the distance between the centre point (centre_lat, centre_lon)
            // and the point at this latitude on the other side (curr_lat, centre_lon + pi)
            // dist( ( lambda0, phi0 ), ( lambda0 + pi, phi1 ) )
            const double dist_on_other_side = distance( longitude.at(Ilon), centre_lat,
                                                        longitude.at(Ilon) + M_PI, curr_lat );

            // dist( ( lambda0, phi0 ), ( lambda0, phi1 ) )
            const double dist_to_lat1 = distance( longitude.at(Ilon), centre_lat,
                                                  longitude.at(Ilon), curr_lat );

            if ( dist_on_other_side <= padded_scale / 2. ) {
                // Use all points on line of latitude phi1
                dlon_N = Nlon;
            } else if ( dist_to_lat1 >= padded_scale / 2 ) {
                // Do not use any points on line of latitude phi1
                dlon_N = 1;
            } else {
                // Else, there must be some del_lon in (0,pi) such that
                // dist( (lambda0, phi0), (lambda0+del_lon, phi1) ) = ell/2
                // solve for that del_lon using the equation for ditance

                const double cos_term = (padded_scale >= 2 * M_PI * constants::R_earth) ? 
                                         -1 : 
                                         cos( (padded_scale / 2.) / constants::R_earth );
                const double numer = cos_term - sin(centre_lat) * sin(curr_lat);
                const double denom = cos(centre_lat) * cos(curr_lat);
                const double del_lon = acos( numer / denom );

                dlon_N = (int) ceil(del_lon / dlon);
            }
        }

        if (constants::PERIODIC_X) {
            LON_lb = Ilon - dlon_N;
            LON_ub = Ilon + dlon_N;
            if (LON_ub - LON_lb >= Nlon) { 
                LON_lb = 0;
                LON_ub = Nlon;
            }
        } else {
            LON_lb = std::max(0,    Ilon - dlon_N);
            LON_ub = std::min(Nlon, Ilon + dlon_N);
        }
    }
}
