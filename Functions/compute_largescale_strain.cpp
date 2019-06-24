#include <algorithm>
#include <vector>
#include <math.h>
#include "../functions.hpp"
#include "../differentiation_tools.hpp"
#include "../constants.hpp"

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

    // For the sake of convenience, compute all of the derivatives first
    std::vector<double*> x_deriv_vals, y_deriv_vals, z_deriv_vals;
    std::vector<const std::vector<double>*> deriv_fields;

    double ux_x, uy_x, uz_x;
    double ux_y, uy_y, uz_y;
    double ux_z, uy_z, uz_z;

    deriv_fields.push_back(&u_x);
    deriv_fields.push_back(&u_y);
    deriv_fields.push_back(&u_z);

    x_deriv_vals.push_back(&ux_x);
    x_deriv_vals.push_back(&uy_x);
    x_deriv_vals.push_back(&uz_x);

    y_deriv_vals.push_back(&ux_y);
    y_deriv_vals.push_back(&uy_y);
    y_deriv_vals.push_back(&uz_y);

    z_deriv_vals.push_back(&ux_z);
    z_deriv_vals.push_back(&uy_z);
    z_deriv_vals.push_back(&uz_z);
    
    Cart_derivatives_at_point(
            x_deriv_vals, y_deriv_vals,
            z_deriv_vals, deriv_fields,
            latitude, longitude,
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
