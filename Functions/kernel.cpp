#include "../functions.hpp"
#include "../constants.hpp"
#include <math.h>

double sinc(const double &x) {
    double s;

    // to avoid some numerical issues, use series expansion for small x
    if (fabs(x) < 1e-9) { s = 1. - x*x/6.; }
    else                { s = sin(x) / x; }

    return s;
}

double kernel(
        const double dist,  /**< [in] Distance as argument to the kernel */
        const double scale  /**< [in] Filtering scale */
        ) {
   
    double kern;
    const double D = dist / ( scale / 2. );

    switch (constants::KERNEL_OPT) {
        case 0: kern = D < 1 ? 1. : 0;
                break;
        case 1: kern = exp( -pow( D, 4) );
                break;
        case 2: kern = exp( -pow( D, 2) );
                break;
        case 3: kern = sinc( M_PI * D );
                break;
        case 4: kern = 0.5 * (1 - tanh( (D - 1) / (0.1) ));
                break;
    }

    #if DEBUG >= 6
    fprintf(stdout, "Kernel(dist=%.4g, scale=%.4g) = %.4g\n",
            dist, scale, kern);
    #endif

    return kern;
}
