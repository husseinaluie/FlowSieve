#include "../functions.hpp"
#include "../constants.hpp"

double latitude_derivative_at_point(
        const double * field,    /**< [in] field to differentiate */
        const double * latitude, /**< [in] (1D) latitude array */
        const int Itime,         /**< [in] time index at which to differentiate */
        const int Idepth,        /**< [in] depth index at which to differentiate */
        const int Ilat,          /**< [in] latitude index at which to differentiate */
        const int Ilon,          /**< [in] longitude at which to differentiate */
        const int Ntime,         /**< [in] size of time dimension */
        const int Ndepth,        /**< [in] size of depth dimension */
        const int Nlat,          /**< [in] size of latitude dimension */
        const int Nlon,          /**< [in] size of longitude dimension */
        const double * mask      /**< [in] (2D) array to distinguish land/water cells */
        ) {

    int UP, DOWN, index;
    double deriv_val = 0.;

    // Differentiation vectors
    double *ddlat;
    ddlat = new double[5];

    // Assuming uniform grid
    double dlat = latitude[1] - latitude[0];

    // Build outwards to try and build the stencil, but stop when
    //   we either hit a land cell or have gone far enough.
    for (DOWN = Ilat; DOWN > 0; DOWN--) {
        if ( (Ilat - DOWN) > constants::DiffOrd + 2 ) { break; }
        index = Index(0, 0, DOWN, Ilon, Ntime, Ndepth, Nlat, Nlon);
        if (mask[index] == 0) { break; }
    }
    for (UP = Ilat; UP < Nlat; UP++) {
        if ( (UP - Ilat) > constants::DiffOrd + 2 ) { break; }
        index = Index(0, 0, UP, Ilon, Ntime, Ndepth, Nlat, Nlon);
        if (mask[index] == 0) { break; }
    }

    // We've possibly made too large of a stencil, so now collapse it back down
    while (UP - DOWN > constants::DiffOrd + 1) {
        if (UP - Ilat > Ilat - DOWN) { UP--; }
        else { DOWN++; }
    }

    if (UP - DOWN > constants::DiffOrd) {
        // If we have enough cells for differentiation, do it
        differentiation_vector(ddlat, dlat, Ilat-DOWN);
        for (int IND = DOWN; IND < UP; IND++) {
            index = Index(Itime, Idepth, IND, Ilon, Ntime, Ndepth, Nlat, Nlon);
            deriv_val +=  field[index] * ddlat[IND-DOWN];
        }
    } else {
        // Otherwise, return 0
        deriv_val = 0.;
    }

    delete[] ddlat;

    return deriv_val;
}
