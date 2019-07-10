#include <stdio.h>
#include <math.h>    
#include <algorithm>
#include <vector>
#include "../constants.hpp"
#include "../functions.hpp"

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

    double dlon_m, delta_lat, local_scale, cosrat;
    int dlon_N;
    
    if (KernPad < 0) {
        // KernPad < 0 means we use the entire domain, so don't bother
        //   calculating an integration region.
        LON_lb = 0;
        LON_ub = Nlon;
    } else {
        // On spherical coordinates, the distance calculator is
        //      very expensive. So when the grid is uniform, we
        //      avoid extra calls to distance with the extra
        //      information that we have.

        // The spacing (in metres and points) betwee longitude gridpoints
        //   The factor of 2 is diameter->radius 
        if (constants::CARTESIAN) { dlon_m = dlon; } 
        else                      { dlon_m = dlon * constants::R_earth * cos(curr_lat); }

        // Next determine how far we need to go (since we're already 
        //   at a finite distance in latitude, we don't need to go
        //   the full scale distance in longitude).
        // Essentially, use 'circular' integration regions, not square
        //    this will reduce the number of cells
        //    required, which should improve performance.
        // The abs in local_scale is to handle the 'comfort zone'
        //    where delta_lat > scale (from the KernPad factor in dlat_N)
        if (constants::CARTESIAN) { 
            delta_lat = centre_lat - curr_lat; 
            // Simply use Cartesian pythagorean
            local_scale = sqrt( fabs( scale*scale - delta_lat*delta_lat ));
        }
        else { 
            delta_lat = constants::R_earth * ( centre_lat - curr_lat ); 
            // Spherical law of cosines, simplified because we know
            //   that lines of longitude are perpendicular to
            //   lines of latitude
            if ( delta_lat < scale ) { 
                cosrat =   cos( scale     / constants::R_earth ) 
                         / cos( delta_lat / constants::R_earth );
            } else {
                cosrat =   cos( delta_lat / constants::R_earth ) 
                         / cos( scale     / constants::R_earth );
            }
            local_scale = acos( cosrat ) * constants::R_earth;
        }


        // Now find the appropriate integration region
        //   The factor of 2 is diameter->radius 
        dlon_N = ceil( ( KernPad * local_scale / dlon_m ) / 2.);
        dlon_N = std::min(Nlon, dlon_N);

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
