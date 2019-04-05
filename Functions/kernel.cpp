#include "../functions.hpp"
#include "../constants.hpp"
#include <math.h>

double kernel(
        const double dist,  /**< [in] Distance as argument to the kernel */
        const double scale  /**< [in] Filtering scale */
        ) {
   
    double kern;

    #if KERNEL_OPT == 0
    if ( dist < (scale / 2) ) {
        kern =  1.;
    } else {
        kern =  0.;
    }
    #elif KERNEL_OPT == 1
    kern = exp( -pow( dist / (scale/2)  , 4) );
    #elif KERNEL_OPT == 2
    kern = exp( -pow( dist / (scale/2)  , 2) );
    #endif

    #if DEBUG >= 6
    fprintf(stdout, "Kernel(dist=%.4g, scale=%.4g) = %.4g\n",
            dist, scale, kern);
    #endif

    return kern;
}
