/*
 *
 * This function computes the distance between two points on
 *   a spherical shell along a great circle.
 *   (see https://en.wikipedia.org/wiki/Great-circle_distance)
 * It should avoid floating point issues for the grid scales that
 *   we're considering. Worst case, increase to quads for better
 *   accuracy.
 *
 */


#include <math.h>    
#include "../constants.hpp"

double dist(const double lon1, const double lat1, const double lon2, const double lat2) {

    double Delta_lon, numer, denom, Delta_sigma, distance;

    Delta_lon = lon2 - lon1;
    
    numer =   pow(                                       cos(lat2) * sin(Delta_lon), 2 ) 
            + pow( cos(lat1) * sin(lat2)  -  sin(lat1) * cos(lat2) * cos(Delta_lon), 2 );
    numer =  sqrt(numer);
    
    denom =        sin(lat1) * sin(lat2)  +  cos(lat1) * cos(lat2) * cos(Delta_lon);

    Delta_sigma = atan2(numer, denom);
    
    distance = constants::R_earth * Delta_sigma;

    return distance;
}
