#include <math.h>
#include <vector>
#include <omp.h>
#include <stdlib.h>
#include <mpi.h>

#include "../constants.hpp"
#include "../functions.hpp"

/*!
 *
 * \brief If the grid includes the poles (lat = 90 degrees), then mask it out
 *
 * Modifies mask in-place
 *
 * @param[in]       latitude                    Latitude grid
 * @param[in,out]   mask                        Mask array to differentiate land/water
 * @param[in]       Ntime,Ndepth,Nlat,Nlon      Dimension sizes
 *
 */
void mask_out_pole(
        const std::vector<double> & latitude,
        std::vector<bool> & mask,
        const int Ntime,
        const int Ndepth,
        const int Nlat,
        const int Nlon
        ){

    int wRank=-1, wSize=-1;
    MPI_Comm_rank( MPI_COMM_WORLD, &wRank );
    MPI_Comm_size( MPI_COMM_WORLD, &wSize );

    int Itime, Idepth, Ilat, Ilon;
    size_t index;
    
    const double D2R = M_PI / 180.;
    const double pole_cut = (90. - 0.05) * D2R;

    std::vector<int> masked_latitudes;

    if (not(constants::CARTESIAN)) {
        #pragma omp parallel default(none) \
            private(Itime, Idepth, Ilat, Ilon, index) \
            shared(latitude, mask, masked_latitudes) \
            firstprivate( Nlon, Nlat, Ntime, Ndepth )
        { 
            #pragma omp for collapse(1) schedule(static)
            for (index = 0; index < mask.size(); ++index) {
                Index1to4(index, Itime, Idepth, Ilat, Ilon,
                                 Ntime, Ndepth, Nlat, Nlon);
                if ( fabs(latitude.at(Ilat)) >= pole_cut) {
                    if ( (Itime == 0) and (Idepth == 0) and (Ilon == 0) ) {
                        masked_latitudes.push_back(Ilat);
                    }
                    mask.at(index) = false;
                }
            }
        }
    }

    #if DEBUG >= 0
    for ( size_t ind = 0; ind < masked_latitudes.size(); ++ind ) {
        Ilat = masked_latitudes.at( ind );
        if (wRank == 0) { fprintf(stdout, "Masking out Ilat %d (too close to pole (%g))\n", Ilat, latitude.at(Ilat)/D2R); }
    }
    #endif
}
