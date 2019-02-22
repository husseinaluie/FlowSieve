/*
 *
 * Convert Spherical velocities to Cartesian
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

void vel_Spher_to_Cart(
            double & u_x, double & u_y, double & u_z,
            const double u_r, const double u_lon, const double u_lat,
            const double lon, const double lat
        ) {

    #if DEBUG >= 4
    fprintf(stdout, "Velocity conversion: Spherical to Cartesian");
    #endif

    u_x =   u_r   * cos(lon) * cos(lat)
          - u_lon * sin(lon)
          - u_lat * cos(lon) * sin(lat);
        
    u_y =   u_r   * sin(lon) * cos(lat)
          + u_lon * cos(lon)
          - u_lat * sin(lon) * sin(lat);
                     
    u_z =   u_r   * sin(lat)
          + u_lat * cos(lat);

}
