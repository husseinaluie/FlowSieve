#include <algorithm>
#include <vector>
#include <math.h>
#include "../functions.hpp"
#include "../constants.hpp"

#ifndef DEBUG
    #define DEBUG 0
#endif

void compute_largescale_strain(
        double & S_xx,                          /**< [in] where to store xx component */
        double & S_xy,                          /**< [in] where to store xy component */
        double & S_xz,                          /**< [in] where to store xz component */
        double & S_yy,                          /**< [in] where to store yy component */
        double & S_yz,                          /**< [in] where to store yz component */
        double & S_zz,                          /**< [in] where to store zz component */
        const std::vector<double> & u_x,        /**< [in] full (4D) u_x for calculation */
        const std::vector<double> & u_y,        /**< [in] full (4D) u_y for calculation */
        const std::vector<double> & u_z,        /**< [in] full (4D) u_z for calculation */
        const int Itime,                        /**< [in] Current position in time dimension */
        const int Idepth,                       /**< [in] Current position in depth dimension */
        const int Ilat,                         /**< [in] Current position in latitude dimension */
        const int Ilon,                         /**< [in] Current position in longitude dimension */
        const int Ntime,                        /**< [in] Length of time dimension */
        const int Ndepth,                       /**< [in] Length of depth dimension */
        const int Nlat,                         /**< [in] Length of latitude dimension */
        const int Nlon,                         /**< [in] Length of longitude dimension */
        const std::vector<double> & longitude,  /**< [in] Longitude dimension (1D) */
        const std::vector<double> & latitude,   /**< [in] Latitude dimension (1D) */
        const std::vector<double> & mask        /**< [in] Mask array (2D) to distinguish land from water*/
        ) {

    //// For the sake of convenience, compute all of the derivatives first
    
    // x derivatives
    double ux_x, uy_x, uz_x;
    ux_x = x_derivative_at_point(u_x, latitude, longitude, 
            Itime, Idepth, Ilat, Ilon,
            Ntime, Ndepth, Nlat, Nlon,
            mask);
    uy_x = x_derivative_at_point(u_y, latitude, longitude, 
            Itime, Idepth, Ilat, Ilon,
            Ntime, Ndepth, Nlat, Nlon,
            mask);
    uz_x = x_derivative_at_point(u_z, latitude, longitude, 
            Itime, Idepth, Ilat, Ilon,
            Ntime, Ndepth, Nlat, Nlon,
            mask);

    // y derivatives
    double ux_y, uy_y, uz_y;
    ux_y = y_derivative_at_point(u_x, latitude, longitude, 
            Itime, Idepth, Ilat, Ilon,
            Ntime, Ndepth, Nlat, Nlon,
            mask);
    uy_y = y_derivative_at_point(u_y, latitude, longitude, 
            Itime, Idepth, Ilat, Ilon,
            Ntime, Ndepth, Nlat, Nlon,
            mask);
    uz_y = y_derivative_at_point(u_z, latitude, longitude, 
            Itime, Idepth, Ilat, Ilon,
            Ntime, Ndepth, Nlat, Nlon,
            mask);

    // z derivatives
    double ux_z, uy_z, uz_z;
    ux_z = z_derivative_at_point(u_x, latitude, longitude, 
            Itime, Idepth, Ilat, Ilon,
            Ntime, Ndepth, Nlat, Nlon,
            mask);
    uy_z = z_derivative_at_point(u_y, latitude, longitude, 
            Itime, Idepth, Ilat, Ilon,
            Ntime, Ndepth, Nlat, Nlon,
            mask);
    uz_z = z_derivative_at_point(u_z, latitude, longitude, 
            Itime, Idepth, Ilat, Ilon,
            Ntime, Ndepth, Nlat, Nlon,
            mask);


    // Now actually compute the strain terms
    S_xx = ux_x;
    S_yy = uy_y;
    S_zz = uz_z;
    S_xy = 0.5 * (ux_y + uy_x);
    S_xz = 0.5 * (ux_z + uz_x);
    S_yz = 0.5 * (uy_z + uz_y);

}

