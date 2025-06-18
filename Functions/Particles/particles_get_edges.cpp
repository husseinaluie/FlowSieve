#include <math.h>
#include <algorithm>
#include <vector>
#include <omp.h>
#include "../../constants.hpp"
#include "../../particles.hpp"

void particles_get_edges(
        int & left,
        int & right,
        int & bottom,
        int & top,
        const double & ref_lat,
        const double & ref_lon,
        const std::vector<double> & lat,
        const std::vector<double> & lon
        ){

    const double lon_min = lon.front(),
                 lon_max = lon.back(),
                 lat_min = lat.front(),
                 lat_max = lat.back();
    
    if ( ( ref_lon < lon_min ) or ( ref_lon > lon_max ) ) {
        left  = lon.size() - 1;
        right = 0;
    } else {
        right = std::upper_bound( lon.begin(), lon.end(), ref_lon ) - lon.begin();
        if (right == 0) {
            left = lon.size() - 1;
        } else {
            left = right - 1;
        }
    }
        
    if ( ref_lat < lat_min ) {
        // If we've gone outside of the lat bounds, then just use a nearest-neighbours
        bottom =  0;
        top    =  0;
    } else if ( ref_lat > lat_max ) {
        // If we've gone outside of the lat bounds, then just use a nearest-neighbours
        bottom = lat.size() - 1;
        top    = lat.size() - 1;
    } else {
        top = std::upper_bound( lat.begin(), lat.end(), ref_lat ) - lat.begin();
        if (top == 0) {
            bottom = 0;
        } else {
            bottom = top - 1;
        }
    }

}
