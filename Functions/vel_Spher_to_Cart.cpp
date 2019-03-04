#include <math.h>
#include "../functions.hpp"

#ifndef DEBUG
    #define DEBUG 0
#endif

void vel_Spher_to_Cart(
            double & u_x,       /**< [in] u_x to be returned */
            double & u_y,       /**< [in] u_y to be returned */
            double & u_z,       /**< [in] u_z to be returned */
            const double u_r,   /**< [in] u_r to convert */
            const double u_lon, /**< [in] u_lon to convert */
            const double u_lat, /**< [in] u_lat to convert */
            const double lon,   /**< [in] longitude of location for conversion */
            const double lat    /**< [in] latitude of location for conversion */
        ) {

    u_x =   u_r   * cos(lon) * cos(lat)
          - u_lon * sin(lon)
          - u_lat * cos(lon) * sin(lat);
        
    u_y =   u_r   * sin(lon) * cos(lat)
          + u_lon * cos(lon)
          - u_lat * sin(lon) * sin(lat);
                     
    u_z =   u_r   * sin(lat)
          + u_lat * cos(lat);

}
