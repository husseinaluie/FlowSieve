#include "../functions.hpp"
#include "../constants.hpp"

/*!
 * \brief Convenience tool to convert physical index (time, depth, lat, lon) to a logical index.
 *
 * Index is a function to convert a four-point (physical) index
 *   (Itime, Idepth, Ilat, Ilon) into a one-point (logical) index
 *   to access the double arrays.
 *
 * Assumes standard CF ordering: time-depth-lat-lon
 *
 * @param[in] Itime,Idepth,Ilat,Ilon    4-indices to be converted to a 1-index
 * @param[in] Ntime,Ndepth,Nlat,Nlon    dimension sizes
 *
 * @returns The effective 1-index that corresponds to the 4-index tuplet
 *
 */
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
