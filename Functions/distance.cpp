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

    double Delta_lon, numer, denom, Delta_sigma, distance;

    Delta_lon = lon2 - lon1;
    
    numer =   pow(                                       cos(lat2) * sin(Delta_lon), 2 ) 
            + pow( cos(lat1) * sin(lat2)  -  sin(lat1) * cos(lat2) * cos(Delta_lon), 2 );
    numer =  sqrt(numer);
    
    denom =        sin(lat1) * sin(lat2)  +  cos(lat1) * cos(lat2) * cos(Delta_lon);

    Delta_sigma = atan2(numer, denom);
    
    distance = constants::R_earth * Delta_sigma;

    #if DEBUG >= 5
    fprintf(stdout, "dist(lon1 = %.5g, lat1 = %.5g, lon2 = %.5g, lat2 = %.5g) = %.5g",
            lon1, lat1, lon2, lat2, distance);
    #endif

    return distance;

}
