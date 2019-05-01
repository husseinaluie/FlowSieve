#include "../functions.hpp"
#include "../constants.hpp"
#include <math.h>

double sinc(const double &x) {
    if (x == 0.) {
        return 1; // account for zero
    } else if (fabs(x) < 1e-10) {
        return 1 - x*x/6; // to avoid some numerical issues, use series expansion for small x
    } else {
        return sin(x) / x;
    }
}

double kernel(
        const double dist,  /**< [in] Distance as argument to the kernel */
        const double scale  /**< [in] Filtering scale */
        ) {
   
    double kern;

    switch (constants::KERNEL_OPT) {
        case 0: if ( dist < (scale / 2) ) { kern =  1.; }
                else { kern =  0.; }
                break;
        case 1: kern = exp( -pow( dist / (scale/2)  , 4) );
                break;
        case 2: kern = exp( -pow( dist / (scale/2)  , 2) );
                break;
        case 3: kern = sinc( dist / (scale/2) );
                break;
    }

    #if DEBUG >= 6
    fprintf(stdout, "Kernel(dist=%.4g, scale=%.4g) = %.4g\n",
            dist, scale, kern);
    #endif

    return kern;
}
