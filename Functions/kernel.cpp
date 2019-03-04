#include "../functions.hpp"

#ifndef DEBUG
    #define DEBUG 0
#endif

double kernel(
        const double dist,  /**< [in] Distance as argument to the kernel */
        const double scale  /**< [in] Filtering scale */
        ) {
   
    double kern;

    if ( dist < (scale / 2) ) {
        kern =  1.;
    } else {
        kern =  0.;
    }

    #if DEBUG >= 6
    fprintf(stdout, "Kernel(dist=%.4g, scale=%.4g) = %.4g\n",
            dist, scale, kern);
    #endif

    return kern;
}
