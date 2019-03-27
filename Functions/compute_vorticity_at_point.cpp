#include <algorithm>
#include <vector>
#include <math.h>
#include "../functions.hpp"
#include "../differentiation_tools.hpp"
#include "../constants.hpp"

void compute_vorticity_at_point(
        double & vort_r_tmp,                    /**< [in] where to store vort_r */
        double & vort_lon_tmp,                  /**< [in] where to store vort_lon */
        double & vort_lat_tmp,                  /**< [in] where to store vort_lat */
        const std::vector<double> & u_r,        /**< [in] full u_r for calculation */
        const std::vector<double> & u_lon,      /**< [in] full u_lon for calculation */
        const std::vector<double> & u_lat,      /**< [in] full u_lat for calculation */
        const int Ntime,                        /**< [in] Length of time dimension */
        const int Ndepth,                       /**< [in] Length of depth dimension */
        const int Nlat,                         /**< [in] Length of latitude dimension */
        const int Nlon,                         /**< [in] Length of longitude dimension */
        const int Itime,                        /**< [in] Current position in time dimension */
        const int Idepth,                       /**< [in] Current position in depth dimension */
        const int Ilat,                         /**< [in] Current position in latitude dimension */
        const int Ilon,                         /**< [in] Current position in longitude dimension */
        const std::vector<double> & longitude,  /**< [in] Longitude dimension (1D) */
        const std::vector<double> & latitude,   /**< [in] Latitude dimension (1D) */
        const std::vector<double> & mask        /**< [in] Mask array (2D) to distinguish land from water*/
        ) {

    // For the moment, only compute vort_r
    vort_r_tmp   = 0.;
    vort_lon_tmp = 0.;
    vort_lat_tmp = 0.;


    #if CARTESIAN
    vort_r_tmp += Cart_derivative_at_point(
            u_lat, latitude, longitude, "x",
            Itime, Idepth, Ilat, Ilon,
            Ntime, Ndepth, Nlat, Nlon,
            mask);
    vort_r_tmp -= Cart_derivative_at_point(
            u_lon, latitude, longitude, "y",
            Itime, Idepth, Ilat, Ilon,
            Ntime, Ndepth, Nlat, Nlon,
            mask);
    #elif not(CARTESIAN)
    // Longitudinal derivative component
    vort_r_tmp += constants::R_earth
                    * spher_derivative_at_point(
                        u_lat, longitude, "lon",
                        Itime, Idepth, Ilat, Ilon,
                        Ntime, Ndepth, Nlat, Nlon,
                        mask);

    // Latitudinal derivative component
    //  - ddlat ( u_lon * cos(lat) ) = u_lon * sin(lat) - ddlat( u_lon ) * cos(lat)
    int index = Index(Itime, Idepth, Ilat, Ilon,
                      Ntime, Ndepth, Nlat, Nlon);

    vort_r_tmp += constants::R_earth
                    * (   sin(latitude.at(Ilat)) * u_lon.at(index)
                        - cos(latitude.at(Ilat)) * spher_derivative_at_point(
                                                    u_lon, latitude, "lat",
                                                    Itime, Idepth, Ilat, Ilon,
                                                    Ntime, Ndepth, Nlat, Nlon,
                                                    mask)
                      );

    // Scale
    vort_r_tmp *= 1. / ( pow(constants::R_earth, 2) * cos(latitude.at(Ilat)) );
    #endif
}

