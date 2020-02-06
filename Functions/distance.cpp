#include <stdio.h>
#include <math.h>    
#include "../constants.hpp"

double distance(
        const double lon1,
        const double lat1,
        const double lon2,
        const double lat2,
        const double Llon,
        const double Llat
        ) {

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

    }

    #if DEBUG >= 6
    fprintf(stdout, "dist(lon1 = %.5g, lat1 = %.5g, lon2 = %.5g, lat2 = %.5g) = %.5g\n",
            lon1, lat1, lon2, lat2, distance);
    #endif

    return distance;

}
