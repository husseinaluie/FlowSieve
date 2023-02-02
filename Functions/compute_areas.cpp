#include <math.h>
#include <vector>
#include <omp.h>

#include "../functions.hpp"
#include "../constants.hpp"

/*! 
 * \brief Compute the area of each computational cell
 *
 * Works for both spherical and Cartesian coordinated (based on constants.hpp)
 * as well as uniform and non-uniform grids (also based on constants.hpp).
 *
 * @param[in,out]   areas       array in which areas will be stored
 * @param[in]       longitude   longitude vector (1D)
 * @param[in]       latitude    latitude vector (1D)
 *
 */
void compute_areas(
        std::vector<double> &areas,
        const std::vector<double> &longitude,
        const std::vector<double> &latitude
        ) {

    double dlat, dlon, coeff;
    int ii, jj;

    // Get the array sizes
    const int Nlon = longitude.size();
    const int Nlat = latitude.size();

    if (constants::UNIFORM_LAT_GRID) {
        // For the moment, assume a uniform grid
        dlat = latitude.at( 1) - latitude.at( 0);
        dlon = longitude.at(1) - longitude.at(0);

        if (constants::CARTESIAN) { coeff = dlat * dlon; }
        else { coeff = pow( constants::R_earth, 2) * dlat * dlon; }

        // Compute the area of each cell
        #pragma omp parallel default(none) private(ii, jj) shared(areas, coeff, latitude) \
        firstprivate( Nlat, Nlon )
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
    } else {

        // Compute the area of each cell
        #pragma omp parallel default(none) private(ii, jj, dlat, dlon)\
        shared(areas, coeff, longitude, latitude) \
        firstprivate( Nlat, Nlon )
        {
            dlon = longitude.at(1) - longitude.at(0);
            #pragma omp for collapse(2) schedule(guided)
            for (ii = 0; ii < Nlat; ii++) {
                for (jj = 0; jj < Nlon; jj++) {

                    // Get the lat grid spacing
                    if (ii == 0) {
                        dlat = ( latitude.at(1) - latitude.at(0) ) / 2.;
                    } else if (ii == Nlat - 1) {
                        dlat = ( latitude.at(Nlat-1) - latitude.at(Nlat-2) ) / 2.;
                    } else {
                        dlat =   ( latitude.at(ii) - latitude.at(ii-1) ) / 2.
                               + ( latitude.at(ii) - latitude.at(ii-1) ) / 2.;
                    }

                    // Get coefficient
                    if (constants::CARTESIAN) { coeff = dlat * dlon; }
                    else { coeff = pow( constants::R_earth, 2) * dlat * dlon; }

                    // Compute cell area
                    if (constants::CARTESIAN) {
                        areas.at(ii*Nlon + jj) = coeff;
                    } else {
                        areas.at(ii*Nlon + jj) = coeff * cos(latitude.at(ii));
                    }
                }
            }
        }
    }

}
