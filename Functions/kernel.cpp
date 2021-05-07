#include "../functions.hpp"
#include "../constants.hpp"
#include <math.h>

/*!
 * \brief sinc function that uses Taylor approximation near zero for stability.
 *
 * This is only used in the case of a sinc kernel
 *
 */
double sinc(const double &x) {
    double s;

    // to avoid some numerical issues, use series expansion for small x
    if (fabs(x) < 1e-9) { s = 1. - x*x/6.; }
    else                { s = sin(x) / x; }

    return s;
}

/*!
 * \brief Primary kernel function coarse-graining procedure (G in publications)
 *
 * @param[in]   distance    distance for evaluating the kernel
 * @param[in]   scale       filter scale (in metres)
 * 
 * @returns The kernel value for a given distance and filter scale
 *
 */
double kernel(
        const double dist,
        const double scale
        ) {
   
    double kern;
    const double D = ( scale > 0 ) ? ( dist / ( scale / 2. ) ) : ( dist == 0 ) ? 1. : 0.;

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
