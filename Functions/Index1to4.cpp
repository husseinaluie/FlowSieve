#include "../functions.hpp"
#include "../constants.hpp"

void Index1to4( 
        const size_t index,    /**< [in]  1-index */
        int & Itime,        /**< [out] Time index */
        int & Idepth,       /**< [out] Depth index */
        int & Ilat,         /**< [out] Latitude index */
        int & Ilon,         /**< [out] Longitude index */
        const int & Ntime,        /**< [out] Length of time dimension */
        const int & Ndepth,       /**< [out] Length of depth dimension */
        const int & Nlat,         /**< [out] Length of latitude dimension */
        const int & Nlon          /**< [out] Length of longitude dimension */
        ){

    /*
    int index =   Itime  * ( Ndepth * Nlat * Nlon )
                + Idepth * (          Nlat * Nlon )
                + Ilat   * (                 Nlon )
                + Ilon;
    */

    //int index =  ( ( Itime * Ndepth + Idepth) * Nlat + Ilat ) * Nlon + Ilon;

    Ilon   = index % Nlon;
    Ilat   =   (index - Ilon) / Nlon  %  Nlat;
    Idepth =  ((index - Ilon) / Nlon - Ilat) / Nlat  %  Ndepth;
    Itime  = (((index - Ilon) / Nlon - Ilat) / Nlat - Idepth) / Ndepth;

    #if DEBUG >= 5
    fprintf(stdout, "      Index conversion %d -> (%d, %d, %d, %d)\n ", 
            index, Itime, Idepth, Ilat, Ilon);
    #endif
}
