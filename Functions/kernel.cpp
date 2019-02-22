#include "../functions.hpp"

// If the DEBUG flag hasn't been set,
//   then use default value of 0 
#ifndef DEBUG
    #define DEBUG 0
#endif

double kernel(const double dist, const double scale) {
   
    double kern;

    if ( dist < (scale / 2) ) {
        kern =  1.;
    } else {
        kern =  0.;
    }

    #if DEBUG >= 5
    fprintf(stdout, "Kernel(dist=%.4g, scale=%.4g) = %.4g\n",
            dist, scale, kern);
    #endif

    return kern;
}
