#include <stdio.h>
#include <math.h>    
#include "../constants.hpp"

/*!
 * \brief Compute the distance (in metres) between two points in the domain.
 *
 * In spherical coordinates,computes the distance between two points on
 *   a spherical shell along a great circle.
 *   (see https://en.wikipedia.org/wiki/Great-circle_distance)
 * It should avoid floating point issues for the grid scales that
 *   we're considering. 
 *
 * The last two arguments (Llon and Llat) give the physical 
 *   length of the two dimensions. This is used in the case
 *   of periodic Cartesian grids. They are otherwise unused.
 *
 * @param[in]   lon1,lat1   coordinates for the first position
 * @param[in]   lon2,lat2   coordinates for the second position
 * @param[in]   Llon,Llat   physical length of the dimensions
 *
 * @returns returns the distance (in metres) between two points.
 *
 */
double distance(
        const double lon1,
        const double lat1,
        const double lon2,
        const double lat2,
        const double Llon,
        const double Llat
        ) {

    // If they're the same point, just return zero. Otherwise, might encounter
    // some floating point issues in the inverse trigs (specifically the arccos call,
    // since rounding might push the argument to be slightly larger than 1.)
    if ( ( lon1 == lon2 ) and ( lat1 == lat2 ) ) { return 0.; }

    double distance;

    if (constants::CARTESIAN) {
        // If we're on a Cartesian grid, then just compute the straight-forward
        //   Cartesian distance. Account for periodicity if necessary
        double del_x, del_y;
        if (constants::PERIODIC_X) { del_x = fmin( fabs(lon1 - lon2), Llon - fabs(lon1 - lon2)); }
        else { del_x = lon1 - lon2; }

        if (constants::PERIODIC_Y) { del_y = fmin( fabs(lat1 - lat2), Llat - fabs(lat1 - lat2)); }
        else { del_y = lat1 - lat2; }

        distance = sqrt( pow(del_x, 2) + pow(del_y, 2) );

    } else {
        // If not on a Cartesian grid, then we're on a spherical grid.
        //   Compute distances along great circles.

        const long double   Delta_lon       = lon2 - lon1,
                            cos_lat2        = cos( lat2 ),
                            sin_lat2        = sin( lat2 ),
                            cos_lat1        = cos( lat1 ),
                            sin_lat1        = sin( lat1 ),
                            cos_Delta_lon   = cos( Delta_lon );
        if ( isinf( lon1 ) ) {
            fprintf( stdout, " %'g, %'g, %'g, %'g \n", lon1, lon2, lat1, lat2 );
        }
        // Since cos is even and cos(2pi - x) = cos(x), we don't need to worry about periodicity in computing Delta_lon for cos(Delta_lon)
        //      for sin(Delta_lon), sin is odd and sin(2pi - x) = -sin(x), but since the result is being squared, it's not a concern either

        long double Delta_sigma;

        if (not(constants::USE_HIGH_PRECISION_DISTANCE)) {
            // This is cheaper, and so long as our distances are at least a couple metres, the floating-point issues shouldn't arise.
            const double acos_argument = sin_lat1 * sin_lat2 + cos_lat1 * cos_lat2 * cos_Delta_lon;
            
            // Handle some rounding cases when distance is nearly maximal
            if ( ( acos_argument < -1 ) and ( fabs(acos_argument + 1) < 1e-10 ) ) { 
                Delta_sigma = M_PI;
            } else {
                Delta_sigma = acos( acos_argument );
            }
        } else {
            // If desired, use the more expensive distance calculator
            long double numer, denom;
            numer =   pow(                                    cos_lat2 * sin(Delta_lon), 2 ) 
                    + pow( cos_lat1 * sin_lat2  -  sin_lat1 * cos_lat2 * cos_Delta_lon , 2 );
            numer =  sqrt(numer);

            denom =        sin_lat1 * sin_lat2  +  cos_lat1 * cos_lat2 * cos_Delta_lon ;

            Delta_sigma = atan2(numer, denom);
        }
        distance = constants::R_earth * Delta_sigma;

    }

    return distance;
}
