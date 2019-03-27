#include <stdio.h>
#include <math.h>    
#include "../constants.hpp"

double distance(
        const double lon1,  /**< [in] Longitude of first position */
        const double lat1,  /**< [in] Latitude of first position */
        const double lon2,  /**< [in] Longitude of second position */
        const double lat2,  /**< [in] Latitude of second position */
        const double Llon,  /**< [in] Physical length of the longitude (first) dimension */
        const double Llat   /**< [in] Physical length of the latitude (second) dimension */
        ) {

    double distance;

    #if CARTESIAN
    // If we're on a Cartesian grid, then just compute the straight-forward
    //   Cartesian distance. Account for periodicity if necessary

    #if PERIODIC_X
    const double del_x = fmin( fabs(lon1 - lon2), Llon - fabs(lon1 - lon2));
    #else
    const double del_x = lon1 - lon2;
    #endif

    #if PERIODIC_Y
    const double del_y = fmin( fabs(lat1 - lat2), Llat - fabs(lat1 - lat2));
    #else
    const double del_y = lat1 - lat2;
    #endif

    distance = sqrt( pow(del_x, 2) + pow(del_y, 2) );

    #elif not(CARTESIAN)
    // If not on a Cartesian grid, then we're on a spherical grid.
    //   Compute distances along great circles.
    double numer, denom, Delta_sigma;

    const double Delta_lon = lon2 - lon1;
    const double cos_lat2 = cos(lat2);
    const double sin_lat2 = sin(lat2);
    const double cos_lat1 = cos(lat1);
    const double sin_lat1 = sin(lat1);
    const double cos_Delta_lon = cos(Delta_lon);

    numer =   pow(                                    cos_lat2 * sin(Delta_lon), 2 ) 
            + pow( cos_lat1 * sin_lat2  -  sin_lat1 * cos_lat2 * cos_Delta_lon , 2 );
    numer =  sqrt(numer);
    
    denom =        sin_lat1 * sin_lat2  +  cos_lat1 * cos_lat2 * cos_Delta_lon ;

    Delta_sigma = atan2(numer, denom);
    
    distance = constants::R_earth * Delta_sigma;

    #endif

    #if DEBUG >= 6
    fprintf(stdout, "dist(lon1 = %.5g, lat1 = %.5g, lon2 = %.5g, lat2 = %.5g) = %.5g\n",
            lon1, lat1, lon2, lat2, distance);
    #endif

    return distance;

}
