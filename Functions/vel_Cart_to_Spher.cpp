/*
 *
 * Convert Cartesian velocities to spherical
 *   velocities.
 *
 * Note: we are using linear velocities (m/s),
 *   not angular (rad/s), so 
 *   u_\lambda = r\cos(\phi)\cdot\hat{u}_\lambda
 *   u_\phi    = r\cdot\hat{u}_\phi
 *
 * Note: we are still using a Spherical
 *   coordinate system, we are only converting
 *   the velocity fields.
 *
 */

#include <math.h>
#include "../functions.hpp"

// If the DEBUG flag hasn't been set,
//   then use default value of 0 
#ifndef DEBUG
    #define DEBUG 0
#endif

void vel_Cart_to_Spher(
            double & u_r, double & u_lon, double & u_lat,
            const double u_x, const double u_y, const double u_z,
            const double lon, const double lat
        ) {

    #if DEBUG >= 4
    fprintf(stdout, "Velocity conversion: Cartesian to Spherical");
    #endif

    u_r   =   u_x * cos(lon) * cos(lat)
            + u_y * sin(lon) * cos(lat)
            + u_z            * sin(lat);
    
    u_lon = - u_x * sin(lon)
            + u_y * cos(lon);
    
    u_lat = - u_x * cos(lon) * sin(lat)
            - u_y * sin(lon) * sin(lat)
            + u_z            * cos(lat);

}
