#include <stdio.h>
#include <math.h>    
#include <algorithm>
#include <vector>
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
        }
        else { 
            // Spherical law of cosines, simplified because we know
            //   that lines of longitude are perpendicular to
            //   lines of latitude
            const double delta_lat = constants::R_earth * fabs( centre_lat - curr_lat ); 
            if ( delta_lat < padded_scale ) { 
                const double cosrat =   cos( padded_scale     / constants::R_earth ) 
                                      / cos( delta_lat        / constants::R_earth );
                local_scale = acos( cosrat ) * constants::R_earth;
            } else {
                local_scale = 0.;
            }
        }


        // Now find the appropriate integration region
        //   The factor of 2 is diameter->radius 
        int dlon_N = std::min( Nlon, (int) ceil( ( local_scale / dlon_m ) / 2.) );

        if (constants::PERIODIC_X) {
            LON_lb = Ilon - dlon_N;
            LON_ub = Ilon + dlon_N;
            if (LON_ub - LON_lb > Nlon) { LON_ub = LON_lb + Nlon; }
        } else {
            LON_lb = std::max(0,    Ilon - dlon_N);
            LON_ub = std::min(Nlon, Ilon + dlon_N);
        }
    }
}
