#include <math.h>
#include "../functions.hpp"

void vel_Spher_to_Cart(
            double & u_x, double & u_y, double & u_z,
            const double u_r, const double u_lon, const double u_lat,
            const double lon, const double lat
        ) {

    u_x =     u_r   * cos(lon) * cos(lat)
            - u_lon * sin(lon)
            - u_lat * cos(lon) * sin(lat);
        
    u_y =     u_r   * sin(lon) * cos(lat)
            + u_lon * cos(lon)
            - u_lat * sin(lon) * sin(lat);
                     
    u_z =   u_r   * sin(lat)
          + u_lat * cos(lat);

}
