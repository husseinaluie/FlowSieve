#include <math.h>
#include "../functions.hpp"
#include "../constants.hpp"

/*!
 * \brief Convert single Cartesian velocity to spherical velocity
 *
 * Convert Cartesian velocities to spherical
 *   velocities.
 *   (u_x, u_y, u_z) -> (u_r, u_lon, u_lat)
 *
 * Note: we are using linear velocities (m/s),
 *   not angular (rad/s), so 
 *   \f{eqnarray*}{
 *   u_\lambda = r\cos(\phi)\cdot\hat{u}_\lambda \\
 *   u_\phi    = r\cdot\hat{u}_\phi
 *   \f}
 *
 * Note: we are still using a Spherical
 *   coordinate system, we are only converting
 *   the velocity fields.
 *
 * @param[in,out]   u_r,u_lon,u_lat     Computed Spherical velocities
 * @param[in]       u_x,u_y,u_z         Cartesian velocities to be converted
 * @param[in]       lon,lat             coordinates of the location of conversion
 *
 */
void vel_Cart_to_Spher_at_point(
            double & u_r,
            double & u_lon,
            double & u_lat,
            const double u_x,
            const double u_y,
            const double u_z,
            const double lon,
            const double lat
        ) {

    if (constants::CARTESIAN) {
        u_lon = u_x;
        u_lat = u_y;
        u_r   = u_z;
    } else {
        const double cos_lon = cos(lon);
        const double cos_lat = cos(lat);
        const double sin_lon = sin(lon);
        const double sin_lat = sin(lat);

        u_r   =   u_x * cos_lon * cos_lat
                + u_y * sin_lon * cos_lat
                + u_z           * sin_lat;

        u_lon = - u_x * sin_lon
                + u_y * cos_lon;

        u_lat = - u_x * cos_lon * sin_lat
                - u_y * sin_lon * sin_lat
                + u_z           * cos_lat;
    }
}
