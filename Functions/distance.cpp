#include <stdio.h>
#include <math.h>    
#include "../constants.hpp"

#ifndef DEBUG
    #define DEBUG 0
#endif

double distance(
        const double lon1, /**< [in] Longitude of first position */
        const double lat1, /**< [in] Latitude of first position */
        const double lon2, /**< [in] Longitude of second position */
        const double lat2  /**< [in] Latitude of second position */
        ) {

    double numer, denom, Delta_sigma, distance;

    const double Delta_lon = lon2 - lon1;
    const double cos_lat2 = cos(lat2);
    const double sin_lat2 = sin(lat2);
    const double cos_lat1 = cos(lat1);
    const double sin_lat1 = sin(lat1);
    const double cos_Delta_lon = cos(Delta_lon);

    /*
    numer =   pow(                                              cos(lat2) * sin(Delta_lon), 2 ) 
                   + pow( cos(lat1) * sin(lat2)  -  sin(lat1) * cos(lat2) * cos(Delta_lon), 2 );
    numer =  sqrt(numer);
                
    denom =        sin(lat1) * sin(lat2)  +  cos(lat1) * cos(lat2) * cos(Delta_lon);
    */
    
    numer =   pow(                                    cos_lat2 * sin(Delta_lon), 2 ) 
            + pow( cos_lat1 * sin_lat2  -  sin_lat1 * cos_lat2 * cos_Delta_lon , 2 );
    numer =  sqrt(numer);
    
    denom =        sin_lat1 * sin_lat2  +  cos_lat1 * cos_lat2 * cos_Delta_lon ;

    Delta_sigma = atan2(numer, denom);
    
    distance = constants::R_earth * Delta_sigma;

    #if DEBUG >= 6
    fprintf(stdout, "dist(lon1 = %.5g, lat1 = %.5g, lon2 = %.5g, lat2 = %.5g) = %.5g\n",
            lon1, lat1, lon2, lat2, distance);
    #endif

    return distance;

}
