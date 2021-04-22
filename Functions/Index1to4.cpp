#include "../functions.hpp"
#include "../constants.hpp"

/*!
 * \brief Convenience tool to convert logical index to physical index (time, depth, lat, lon).
 *
 * Index is a function to convert a one-point (logical) index
 *   into a four-point (physical) index (Itime, Idepth, Ilat, Ilon)
 *   to access the double arrays.
 *
 * Assumes standard CF ordering: time-depth-lat-lon
 *
 * @param[in]       index                   logical index to be converted to a 4-index
 * @param[in,out]   Itime,Idepth,Ilat,Ilon  4-indices to be returned
 * @param[in]       Ntime,Ndepth,Nlat,Nlon  dimension sizes
 *
 */
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

    Ilon   = index % Nlon;
    Ilat   =   (index - Ilon) / Nlon  %  Nlat;
    Idepth =  ((index - Ilon) / Nlon - Ilat) / Nlat  %  Ndepth;
    Itime  = (((index - Ilon) / Nlon - Ilat) / Nlat - Idepth) / Ndepth;

}
