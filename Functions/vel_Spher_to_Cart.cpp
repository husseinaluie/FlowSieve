#include <math.h>
#include "../functions.hpp"
#include "../constants.hpp"

void vel_Spher_to_Cart(
            double & u_x,
            double & u_y,
            double & u_z,
            const double u_r,
            const double u_lon,
            const double u_lat,
            const double lon,
            const double lat
        ) {

    if (constants::CARTESIAN) {
        u_x = u_lon;
        u_y = u_lat;
        u_z = u_r;
    } else {
        const double cos_lon = cos(lon);
        const double cos_lat = cos(lat);
        const double sin_lon = sin(lon);
        const double sin_lat = sin(lat);
        u_x =   u_r * cos_lon * cos_lat
            - u_lon * sin_lon
            - u_lat * cos_lon * sin_lat;

        u_y =   u_r * sin_lon * cos_lat
            + u_lon * cos_lon
            - u_lat * sin_lon * sin_lat;

        u_z =   u_r           * sin_lat
            + u_lat           * cos_lat;
    }
}
