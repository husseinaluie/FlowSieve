#include <math.h>
#include <vector>
#include <omp.h>

#include "../functions.hpp"
#include "../constants.hpp"

void compute_areas(
        std::vector<double> &areas,             /**< [in] array in which areas will be stored */
        const std::vector<double> &longitude,   /**< [in] array containing longitude dimension (1D) */
        const std::vector<double> &latitude     /**< [in] array containing latitude dimension (1D) */
        ) {

    // For the moment, assume a uniform grid
    const double dlat = latitude.at( 1) - latitude.at( 0);
    const double dlon = longitude.at(1) - longitude.at(0);

    // Get the array sizes
    const int Nlon = longitude.size();
    const int Nlat = latitude.size();

    double coeff;
    if (constants::CARTESIAN) { coeff = dlat * dlon; }
    else { coeff = pow( constants::R_earth, 2) * dlat * dlon; }

    // Compute the area of each cell
    int ii, jj;
    #pragma omp parallel default(none) private(ii, jj) shared(areas, coeff, latitude)
    {
        #pragma omp for collapse(2) schedule(static)
        for (ii = 0; ii < Nlat; ii++) {
            for (jj = 0; jj < Nlon; jj++) {
                if (constants::CARTESIAN) {
                    areas.at(ii*Nlon + jj) = coeff;
                } else {
                    areas.at(ii*Nlon + jj) = coeff * cos(latitude.at(ii));
                }
            }
        }
    }
}
