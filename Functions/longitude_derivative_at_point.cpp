#include "../functions.hpp"
#include "../constants.hpp"

double longitude_derivative_at_point(
        const double * field,     /**< [in] field to differentiate */
        const double * longitude, /**< [in] (1D) longitude array */
        const int Itime,          /**< [in] time index at which to differentiate */
        const int Idepth,         /**< [in] depth index at which to differentiate */
        const int Ilat,           /**< [in] latitude index at which to differentiate */
        const int Ilon,           /**< [in] longitude at which to differentiate */
        const int Ntime,          /**< [in] size of time dimension */
        const int Ndepth,         /**< [in] size of depth dimension */
        const int Nlat,           /**< [in] size of latitude dimension */
        const int Nlon,           /**< [in] size of longitude dimension */
        const double * mask       /**< [in] (2D) array to distinguish land/water cells */
        ) {

    int LEFT, RIGHT, index;
    double deriv_val = 0.;

    // Differentiation vector
    double *ddlon;
    ddlon = new double[5];

    // Assuming a uniform grid
    double dlon = longitude[1] - longitude[0];

    // Build outwards to try and build the stencil, but stop when
    //   we either hit a land cell or have gone far enough.
    for (LEFT = Ilon; LEFT > 0; LEFT--) {
        if ( (Ilon - LEFT) > constants::DiffOrd + 2 ) { break; }
        index = Index(0, 0, Ilat, LEFT, Ntime, Ndepth, Nlat, Nlon);
        if (mask[index] == 0) { break; }
    }
    for (RIGHT = Ilon; RIGHT < Nlon; RIGHT++) {
        if ( (RIGHT - Ilon) > constants::DiffOrd + 2 ) { break; }
        index = Index(0, 0, Ilat, RIGHT, Ntime, Ndepth, Nlat, Nlon);
        if (mask[index] == 0) { break; }
    }

    // We've possibly made too large of a stencil, so now collapse it back down
    while (RIGHT - LEFT > constants::DiffOrd + 1) {
        if (RIGHT - Ilon > Ilon - LEFT) { RIGHT--; }
        else { LEFT++; }
    }

    if (RIGHT - LEFT > constants::DiffOrd) {
        // If we have enough cells for differentiation, do it
        differentiation_vector(ddlon, dlon, Ilon-LEFT);
        for (int IND = LEFT; IND < RIGHT; IND++) {
            index = Index(Itime, Idepth, Ilat, IND, Ntime, Ndepth, Nlat, Nlon);
            deriv_val += ddlon[IND-LEFT] * field[index];
        }
    } else {
        // Otherwise, return 0
        deriv_val = 0.;
    }

    delete[] ddlon;

    return deriv_val;
}
