#include <math.h>
#include "../functions.hpp"

// If the DEBUG flag hasn't been set,
//   then use default value of 0 
#ifndef DEBUG
    #define DEBUG 0
#endif

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

    u_r   =   u_x * cos(lon) * cos(lat)
            + u_y * sin(lon) * cos(lat)
            + u_z            * sin(lat);
    
    u_lon = - u_x * sin(lon)
            + u_y * cos(lon);
    
    u_lat = - u_x * cos(lon) * sin(lat)
            - u_y * sin(lon) * sin(lat)
            + u_z            * cos(lat);

}
