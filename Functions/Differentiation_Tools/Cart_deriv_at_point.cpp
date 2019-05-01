#include <vector>
#include "../../differentiation_tools.hpp"
#include "../../constants.hpp"
#include "../../functions.hpp"
#include <math.h>
#include <string>

double Cart_derivative_at_point(
        const std::vector<double> & field,      /**< [in] field to differentiate */
        const std::vector<double> & latitude,   /**< [in] (1D) latitude array */
        const std::vector<double> & longitude,  /**< [in] (1D) longitude array */
        const std::string & dim,                /**< [in] "x", "y", or "z" (dimension along which to differentiate) */
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

    // Currently assuming ddr = 0
    // ddx = -( sin(lon) / (r cos(lat)) ) * ddlon   - ( cos(lon) * sin(lat) / r ) * ddlat
    // ddy =  ( cos(lon) / (r cos(lat)) ) * ddlon   - ( sin(lon) * sin(lat) / r ) * ddlat
    // ddz =                                          (            cos(lat) / r ) * ddlat
    double dfield_dlon, dfield_dlat;

    dfield_dlon = spher_derivative_at_point(
                        field, longitude, "lon",
                        Itime, Idepth, Ilat, Ilon,
                        Ntime, Ndepth, Nlat, Nlon,
                        mask);

    dfield_dlat = spher_derivative_at_point(
                        field, latitude, "lat",
                        Itime, Idepth, Ilat, Ilon,
                        Ntime, Ndepth, Nlat, Nlon,
                        mask);

    double deriv_val = 0.;

    if (not(constants::CARTESIAN)) {
        const double lon = longitude.at(Ilon);
        const double lat = latitude.at(Ilat);
        const double r = constants::R_earth;
        const double sin_lon = sin(lon);
        const double sin_lat = sin(lat);
        const double cos_lon = cos(lon);
        const double cos_lat = cos(lat);
        if (dim == "x") {
            deriv_val = - ( sin_lon / (r * cos_lat ) ) * dfield_dlon
                        - ( cos_lon * sin_lat / r )    * dfield_dlat;
        } else if (dim == "y") {
            deriv_val =   ( cos_lon / (r * cos_lat ) ) * dfield_dlon
                        - ( sin_lon * sin_lat / r )    * dfield_dlat;
        } else if (dim == "z") {
            deriv_val = ( cos_lat / r ) * dfield_dlat;
        }
    } else {
        // If we're on an underlying Cartesian grid, then
        //    the 'spherical' derivatives are actually the
        //    Cartesian derivatives, so we don't need to 
        //    do anything.
        if      (dim == "x") { deriv_val = dfield_dlon; } 
        else if (dim == "y") { deriv_val = dfield_dlat; } 
        else if (dim == "z") { deriv_val = 0.; }
    }

    return deriv_val;
}
