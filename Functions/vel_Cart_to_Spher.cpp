#include <math.h>
#include "../functions.hpp"
#include "../constants.hpp"

void vel_Cart_to_Spher(
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
