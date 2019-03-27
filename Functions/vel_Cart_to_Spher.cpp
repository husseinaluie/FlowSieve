#include <math.h>
#include "../functions.hpp"
#include "../constants.hpp"

void vel_Cart_to_Spher(
            double & u_r,     /**< [in] u_r to be returned */
            double & u_lon,   /**< [in] u_lon to be returned */
            double & u_lat,   /**< [in] u_lat to be returned */
            const double u_x, /**< [in] u_x to convert */
            const double u_y, /**< [in] u_y to convert */
            const double u_z, /**< [in] u_z to convert */
            const double lon, /**< [in] longitude of location for conversion */
            const double lat  /**< [in] latitude of location for conversion */
        ) {

    #if CARTESIAN
    u_lon = u_x;
    u_lat = u_y;
    u_r   = u_z;
    #elif not(CARTESIAN)
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
    #endif

}
