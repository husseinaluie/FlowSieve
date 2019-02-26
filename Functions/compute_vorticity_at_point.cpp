#include <algorithm>
#include <math.h>
#include "../functions.hpp"
#include "../constants.hpp"

#ifndef DEBUG
    #define DEBUG 0
#endif

void compute_vorticity_at_point(
        double & vort_r_tmp,      /**< [in] where to store vort_r */
        double & vort_lon_tmp,    /**< [in] where to store vort_lon */
        double & vort_lat_tmp,    /**< [in] where to store vort_lat */
        const double * u_r,       /**< [in] full u_r for calculation */
        const double * u_lon,     /**< [in] full u_lon for calculation */
        const double * u_lat,     /**< [in] full u_lat for calculation */
        const int Ntime,          /**< [in] Length of time dimension */
        const int Ndepth,         /**< [in] Length of depth dimension */
        const int Nlat,           /**< [in] Length of latitude dimension */
        const int Nlon,           /**< [in] Length of longitude dimension */
        const int Itime,          /**< [in] Current position in time dimension */
        const int Idepth,         /**< [in] Current position in depth dimension */
        const int Ilat,           /**< [in] Current position in latitude dimension */
        const int Ilon,           /**< [in] Current position in longitude dimension */
        const double * longitude, /**< [in] Longitude dimension (1D) */
        const double * latitude,  /**< [in] Latitude dimension (1D) */
        const double * mask       /**< [in] Mask array (2D) to distinguish land from water*/
        ) {


    vort_r_tmp   = 0.;
    vort_lon_tmp = 0.;
    vort_lat_tmp = 0.;

    // Differentiation vectors
    double *ddlon, *ddlat;
    ddlon = new double[5];
    ddlat = new double[5];

    // Assuming uniform grid
    double dlon = longitude[1] - longitude[0];
    double dlat = latitude[1]  - latitude[0];

    //
    //// Longitudinal derivative component
    //
    int index;

    // Find E-W spacing for differentiation
    //   (stop when we hit the first land cell
    //     or when we have enough )
    int LEFT, RIGHT;

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

    // If we have enough cells for differentiation, do it
    if (RIGHT - LEFT > constants::DiffOrd) {
        differentiation_vector(ddlon, dlon, Ilon-LEFT);
        for (int IND = LEFT; IND < RIGHT; IND++) {
            index = Index(Itime, Idepth, Ilat, IND, Ntime, Ndepth, Nlat, Nlon);
            vort_r_tmp += constants::R_earth * ddlon[IND-LEFT] * u_lat[index];
        }
    }

    //
    //// Latitudinal derivative component
    //

    // Find N-S spacing for differentiation
    int UP, DOWN;

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

    // If we have enough cells for differentiation, do it
    if (UP - DOWN > constants::DiffOrd) {
        differentiation_vector(ddlat, dlat, Ilat-DOWN);
        for (int IND = DOWN; IND < UP; IND++) {
            index = Index(Itime, Idepth, IND, Ilon, Ntime, Ndepth, Nlat, Nlon);
            vort_r_tmp +=  -constants::R_earth * u_lon[index] * ddlat[IND-DOWN] * cos(latitude[Ilat]);
        }
        index = Index(Itime, Idepth, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon);
        vort_r_tmp += constants::R_earth * u_lon[index] * sin(latitude[Ilat]);
    }

    // Scale
    vort_r_tmp *= 1. / ( pow(constants::R_earth, 2) * cos(latitude[Ilat]) );

    // Free the differentiation arrays
    delete[] ddlon;
    delete[] ddlat;

}

