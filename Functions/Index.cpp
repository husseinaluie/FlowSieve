#include "../functions.hpp"
#include "../constants.hpp"

int Index( 
        const int Itime,
        const int Idepth,
        const int Ilat,
        const int Ilon,
        const int Ntime,
        const int Ndepth,
        const int Nlat,
        const int Nlon
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
