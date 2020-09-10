#include <math.h>
#include <vector>
#include <omp.h>
#include <stdlib.h>

#include "../constants.hpp"
#include "../functions.hpp"

void mask_out_pole(
        const std::vector<double> & latitude,
        std::vector<bool> & mask,
        const int Ntime,
        const int Ndepth,
        const int Nlat,
        const int Nlon
        ){

    int Itime, Idepth, Ilat, Ilon;
    size_t index;
    
    // quarter of a degree from pole
    //const double pole_cut = M_PI * ( 1. - (0.25 / 180.) );
    const double pole_cut = M_PI * ( 0.5 - (0.5 / 180.) );

    #if not(CARTESIAN)
    #pragma omp parallel default(none) \
        private(Itime, Idepth, Ilat, Ilon, index) \
        shared(latitude, mask)
    { 
        // Mask out top quarter of a degree
        #pragma omp for collapse(1) schedule(static)
        for (index = 0; index < mask.size(); ++index) {
            Index1to4(index, Itime, Idepth, Ilat, Ilon,
                             Ntime, Ndepth, Nlat, Nlon);
            if ( fabs(latitude.at(Ilat)) >= pole_cut) {
                mask.at(index) = false;
            }
        }
    }
    #endif
}
