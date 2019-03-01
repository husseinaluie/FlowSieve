#include "../functions.hpp"
#include "../constants.hpp"
#include <math.h>

double y_derivative_at_point(
        const double * field,     /**< [in] field to differentiate */
        const double * latitude,  /**< [in] (1D) latitude array */
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

    // Currently assuming ddr = 0
    // ddy = ( cos(lon) / (r cos(lat)) ) * ddlon   - ( sin(lon) * sin(lat) / r ) * ddlat
    double dfield_dlon, dfield_dlat;
    double lon = longitude[Ilon];
    double lat = latitude[Ilat];
    double r = constants::R_earth;

    dfield_dlon = longitude_derivative_at_point(
                        field, longitude,
                        Itime, Idepth, Ilat, Ilon,
                        Ntime, Ndepth, Nlat, Nlon,
                        mask);

    dfield_dlat = latitude_derivative_at_point(
                        field, latitude,
                        Itime, Idepth, Ilat, Ilon,
                        Ntime, Ndepth, Nlat, Nlon,
                        mask);

    double deriv_val =   ( cos(lon) / (r * cos(lat) ) ) * dfield_dlon
                       - ( sin(lon) * sin(lat) / r )    * dfield_dlat;

    return deriv_val;
}
