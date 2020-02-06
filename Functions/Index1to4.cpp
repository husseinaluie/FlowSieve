#include "../functions.hpp"
#include "../constants.hpp"

void Index1to4( 
        const size_t index,
        int & Itime,
        int & Idepth,
        int & Ilat,
        int & Ilon,
        const int & Ntime,
        const int & Ndepth,
        const int & Nlat,
        const int & Nlon
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
