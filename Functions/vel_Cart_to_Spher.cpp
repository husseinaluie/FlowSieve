#include <math.h>
#include "../functions.hpp"

void vel_Cart_to_Spher(
            double & u_r, double & u_lon, double & u_lat,
            const double u_x, const double u_y, const double u_z,
            const double lon, const double lat
        ) {

    u_r   =   u_x * cos(lon) * cos(lat)
            + u_y * sin(lon) * cos(lat)
            + u_z            * sin(lat);
    
    u_lon = - u_x * sin(lon)
            + u_y * cos(lon);
    
    u_lat =  - u_x * cos(lon) * sin(lat)
             - u_y * sin(lon) * sin(lat)
             + u_z            * cos(lat);

}
