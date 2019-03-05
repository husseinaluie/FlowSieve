#include <vector>
#include "../../differentiation_tools.hpp"
#include "../../constants.hpp"
#include "../../functions.hpp"

double latitude_derivative_at_point(
        const std::vector<double> & field,      /**< [in] field to differentiate */
        const std::vector<double> & latitude,   /**< [in] (1D) latitude array */
        const int Itime,                        /**< [in] time index at which to differentiate */
        const int Idepth,                       /**< [in] depth index at which to differentiate */
        const int Ilat,                         /**< [in] latitude index at which to differentiate */
        const int Ilon,                         /**< [in] longitude at which to differentiate */
        const int Ntime,                        /**< [in] size of time dimension */
        const int Ndepth,                       /**< [in] size of depth dimension */
        const int Nlat,                         /**< [in] size of latitude dimension */
        const int Nlon,                         /**< [in] size of longitude dimension */
        const std::vector<double> & mask        /**< [in] (2D) array to distinguish land/water cells */
        ) {

    int UP, DOWN, index;
    double deriv_val = 0.;

    // Differentiation vectors
    std::vector<double> ddlat(5);

    // Assuming uniform grid
    const double dlat = latitude.at(1) - latitude.at(0);

    // Build outwards to try and build the stencil, but stop when
    //   we either hit a land cell or have gone far enough.
    for (DOWN = Ilat; DOWN > 0; DOWN--) {
        if ( (Ilat - DOWN) > constants::DiffOrd + 2 ) { break; }
        index = Index(0, 0, DOWN, Ilon, Ntime, Ndepth, Nlat, Nlon);
        if (mask.at(index) == 0) { break; }
    }
    for (UP = Ilat; UP < Nlat; UP++) {
        if ( (UP - Ilat) > constants::DiffOrd + 2 ) { break; }
        index = Index(0, 0, UP, Ilon, Ntime, Ndepth, Nlat, Nlon);
        if (mask.at(index) == 0) { break; }
    }

    // We've possibly made too large of a stencil, so now collapse it back down
    while (UP - DOWN > constants::DiffOrd) {
        if (UP - Ilat > Ilat - DOWN) { UP--; }
        else { DOWN++; }
    }

    if (UP - DOWN >= constants::DiffOrd) {
        // If we have enough cells for differentiation, do it
        differentiation_vector(ddlat, dlat, Ilat-DOWN);
        for (int IND = DOWN; IND < UP; IND++) {
            index = Index(Itime, Idepth, IND, Ilon, Ntime, Ndepth, Nlat, Nlon);
            deriv_val +=  field.at(index) * ddlat.at(IND-DOWN);
        }
    } else {
        // Otherwise, return 0
        deriv_val = 0.;
    }

    return deriv_val;
}
