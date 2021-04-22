#include <math.h>
#include <vector>
#include <omp.h>
#include <stdlib.h>

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

    int Itime, Idepth, Ilat, Ilon;
    size_t index;
    
    const double D2R = M_PI / 180.;
    const double pole_cut = (90. - 0.1) * D2R;

    if (not(constants::CARTESIAN)) {
        #pragma omp parallel default(none) \
            private(Itime, Idepth, Ilat, Ilon, index) \
            shared(latitude, mask)
        { 
            #pragma omp for collapse(1) schedule(static)
            for (index = 0; index < mask.size(); ++index) {
                Index1to4(index, Itime, Idepth, Ilat, Ilon,
                                 Ntime, Ndepth, Nlat, Nlon);
                if ( fabs(latitude.at(Ilat)) >= pole_cut) {
                    mask.at(index) = false;
                }
            }
        }
    }
}
