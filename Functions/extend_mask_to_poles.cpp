#include <math.h>
#include <algorithm>
#include <vector>
#include <cassert>
#include <omp.h>
#include "../functions.hpp"
#include "../constants.hpp"


void extend_mask_to_poles(
        std::vector<bool> & mask_to_extend,
        const dataset & source_data,
        const std::vector<double> & extended_latitude,
        const int Ilat_start,
        const bool extend_val
        ) {

    // Get grid sizes
    const size_t    Ntime   = source_data.Ntime,
                    Ndepth  = source_data.Ndepth,
                    Nlat    = source_data.Nlat,
                    Nlon    = source_data.Nlon;

    const size_t extended_Nlat = extended_latitude.size();

    const size_t size_to_extend = mask_to_extend.size(),
                 extended_size  = Ntime * Ndepth * extended_Nlat * Nlon;

    // Start extended mask if 'land' everywhere if using land,
    //  otherwise zero-velocity water
    std::vector<bool> extended_mask( extended_size, extend_val );

    int Itime, Idepth, Ilat, Ilon;
    size_t extended_index, index;

    /*
    #if DEBUG >= 2
    fprintf( stdout, "    Original mask size: %'zu, vs extended size: %'zu\n",
                mask_to_extend.size(), extended_mask.size() );
    fprintf( stdout, "    Orininal lat (%'zu) vs extended lat (%'zu) and shifted start (%d)\n",
                     Nlat, extended_latitude.size(), Ilat_start );
    #endif
    */

    // Copy from the original mask into a padded one
    #pragma omp parallel \
    default(none) \
    shared( extended_mask, mask_to_extend ) \
    private( index, Itime, Idepth, Ilat, Ilon, extended_index ) \
    firstprivate( size_to_extend, Nlon, Nlat, Ndepth, Ntime, extended_Nlat,\
                    Ilat_start )
    {
        #pragma omp for collapse(1) schedule(static)
        for ( index = 0; index < size_to_extend; index++ ) {
            Index1to4( index, Itime, Idepth, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon );

            extended_index = Index( Itime, Idepth, Ilat + Ilat_start, Ilon, Ntime, Ndepth, extended_Nlat, Nlon );

            extended_mask.at( extended_index ) = mask_to_extend.at( index );
        }
    }

    // Now swap out the contents of the extended mask 
    //  into the mask_to_extend
    mask_to_extend.swap( extended_mask );

}

