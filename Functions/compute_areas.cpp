#include <math.h>
#include <vector>

#include "../functions.hpp"
#include "../constants.hpp"

#ifndef DEBUG
    #define DEBUG 0
#endif

void compute_areas(
        std::vector<double> &areas,           /**< [in] array in which areas will be stored */
        const std::vector<double> &longitude, /**< [in] array containing longitude dimension (1D) */
        const std::vector<double> &latitude  /**< [in] array containing latitude dimension (1D) */
        ) {

    // For the moment, assume a uniform grid
    double dlat = latitude.at( 1) - latitude.at( 0);
    double dlon = longitude.at(1) - longitude.at(0);

    // Get the array sizes
    const int Nlon = longitude.size();
    const int Nlat = latitude.size();

    double coeff = pow( constants::R_earth, 2) * dlat * dlon;

    // Compute the area of each cell
    for (int ii = 0; ii < Nlat; ii++) {
        for (int jj = 0; jj < Nlon; jj++) {
            areas.at(ii*Nlon + jj) = coeff * cos(latitude.at(ii));
        }
    }

    #if DEBUG >= 2
    fprintf(stdout, "  finished computing areas.\n\n");
    #endif
}
