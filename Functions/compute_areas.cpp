#include <math.h>

#include "../functions.hpp"
#include "../constants.hpp"

#ifndef DEBUG
    #define DEBUG 0
#endif

void compute_areas(
        double * areas,           /**< [in] array in which areas will be stored */
        const double * longitude, /**< [in] array containing longitude dimension (1D) */
        const double * latitude,  /**< [in] array containing latitude dimension (1D) */
        const int Nlon,           /**< [in] size of the longitude dimension */
        const int Nlat            /**< [in] size of the latitude dimension */
        ) {

    // For the moment, assume a uniform grid
    double dlat = latitude[ 1] - latitude[ 0];
    double dlon = longitude[1] - longitude[0];
    double LAT;

    double R_earth = constants::R_earth;

    // Compute the area of each cell
    for (int ii = 0; ii < Nlat; ii++) {
        LAT = latitude[ii];
        for (int jj = 0; jj < Nlon; jj++) {
            areas[ii*Nlon + jj] = pow( R_earth, 2 ) * cos(LAT) * dlat * dlon;
        }
    }

    #if DEBUG >= 2
    fprintf(stdout, "  finished computing areas.\n\n");
    #endif
}
