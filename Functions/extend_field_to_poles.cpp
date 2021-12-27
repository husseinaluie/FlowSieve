#include <math.h>
#include <algorithm>
#include <vector>
#include <cassert>
#include "../functions.hpp"
#include "../constants.hpp"


void extend_field_to_poles(
        std::vector<double> & array_to_extend,
        const dataset & source_data,
        const std::vector<double> & extended_latitude,
        const int Ilat_start
        ) {

    // Get grid sizes
    const size_t    Ntime   = source_data.Ntime,
                    Ndepth  = source_data.Ndepth,
                    Nlat    = source_data.Nlat,
                    Nlon    = source_data.Nlon;

    const size_t extended_Nlat = extended_latitude.size();

    const size_t size_to_extend = array_to_extend.size(),
                 extended_size  = Ntime * Ndepth * extended_Nlat * Nlon;

    std::vector<double> extended_array( extended_size, 0. );

    int Itime, Idepth, Ilat, Ilon;
    size_t extended_index;

    // Copy from the original array into a padded one
    for ( size_t index = 0; index < size_to_extend; index++ ) {
        Index1to4( index, Itime, Idepth, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon );

        extended_index = Index( Itime, Idepth, Ilat + Ilat_start, Ilon, Ntime, Ndepth, extended_Nlat, Nlon );

        extended_array.at( extended_index ) = array_to_extend.at( index );
    }

    // Now swap out the contents of the extended array
    //  into the array_to_extend
    array_to_extend.swap( extended_array );

}
