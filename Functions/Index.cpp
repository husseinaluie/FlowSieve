#include "../functions.hpp"
#include "../constants.hpp"

int Index( 
        const int Itime,  /**< [in] Time index */
        const int Idepth, /**< [in] Depth index */
        const int Ilat,   /**< [in] Latitude index */
        const int Ilon,   /**< [in] Longitude index */
        const int Ntime,  /**< [in] Length of time dimension */
        const int Ndepth, /**< [in] Length of depth dimension */
        const int Nlat,   /**< [in] Length of latitude dimension */
        const int Nlon    /**< [in] Length of longitude dimension */
        ){

    /*
    int index =   Itime  * ( Ndepth * Nlat * Nlon )
                + Idepth * (          Nlat * Nlon )
                + Ilat   * (                 Nlon )
                + Ilon;
    */

    int index =  ( ( Itime * Ndepth + Idepth) * Nlat + Ilat ) * Nlon + Ilon;

    #if DEBUG >= 5
    fprintf(stdout, "      Index conversion (%d, %d, %d, %d) -> %d\n ", 
            Itime, Idepth, Ilat, Ilon, index);
    #endif

    return index;

}
