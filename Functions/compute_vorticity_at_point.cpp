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
        const std::vector<bool> & mask          /**< [in] Mask array (2D) to distinguish land from water*/
        ) {

    // For the moment, only compute vort_r
    vort_r_tmp   = 0.;
    vort_lon_tmp = 0.;
    vort_lat_tmp = 0.;

    std::vector<const std::vector<double>*> deriv_fields {&u_lon, &u_lat, &u_r};

    if (constants::CARTESIAN) {

        double ux_y, ux_z, uy_x, uy_z, uz_x, uz_y;
        std::vector<double*> x_deriv_vals {NULL,  &uy_x, &uz_x};
        std::vector<double*> y_deriv_vals {&ux_y, NULL,  &uz_y};
        std::vector<double*> z_deriv_vals {&ux_z, &uy_z, NULL };

        Cart_derivatives_at_point(
           x_deriv_vals, y_deriv_vals,
           z_deriv_vals, deriv_fields,
           latitude, longitude,
           Itime, Idepth, Ilat, Ilon,
           Ntime, Ndepth, Nlat, Nlon,
           mask);

        vort_lon_tmp = uz_y - uy_z;
        vort_lat_tmp = ux_z - uz_x;
        vort_r_tmp   = uy_x - ux_y;
    } else {

        int index = Index(Itime, Idepth, Ilat, Ilon,
                          Ntime, Ndepth, Nlat, Nlon);
        double ur_lon, ur_lat, ulon_r, ulon_lat, ulat_r, ulat_lon;
        double lat = latitude.at(Ilat);

        std::vector<double*> lon_deriv_vals {NULL,      &ulat_lon, &ur_lon};
        std::vector<double*> lat_deriv_vals {&ulon_lat, NULL,      &ur_lat};
        std::vector<double*> r_deriv_vals   {&ulon_r,   &ulat_r,   NULL   };

        spher_derivative_at_point(
                lat_deriv_vals, deriv_fields,
                latitude, "lat",
                Itime, Idepth, Ilat, Ilon,
                Ntime, Ndepth, Nlat, Nlon,
                mask);

        spher_derivative_at_point(
                lon_deriv_vals, deriv_fields,
                longitude, "lon",
                Itime, Idepth, Ilat, Ilon,
                Ntime, Ndepth, Nlat, Nlon,
                mask);

        // Longitudinal derivative component
        vort_r_tmp += constants::R_earth * ulat_lon;

        // Latitudinal derivative component
        //  - ddlat ( u_lon * cos(lat) ) = u_lon * sin(lat) - ddlat( u_lon ) * cos(lat)
        vort_r_tmp += constants::R_earth
            * (sin(lat) * u_lon.at(index)  -  cos(lat) * ulon_lat);

        // Scale
        vort_r_tmp *= 1. / ( pow(constants::R_earth, 2) * cos(lat) );
    }
}
