#include "../functions.hpp"
#include "../constants.hpp"

size_t Index( 
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

    size_t index =  ( ( (size_t) Itime * (size_t) Ndepth + (size_t) Idepth) * (size_t) Nlat + (size_t) Ilat ) * (size_t) Nlon + (size_t) Ilon;

    #if DEBUG >= 5
    fprintf(stdout, "      Index conversion (%d, %d, %d, %d) -> %d\n ", 
            Itime, Idepth, Ilat, Ilon, index);
    #endif

    return index;

}
