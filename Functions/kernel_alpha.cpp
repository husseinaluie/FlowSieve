#include "../functions.hpp"
#include "../constants.hpp"
#include <math.h>


double kernel_alpha(void) {

    // Get integration bounds
    const double LB = 0.;
    // If we have an 'unbounded' kernel, then use 100, and hope it's enough
    // Otherwise, use the pre-set kernel bounds
    const double UB = (constants::KernPad < 0) ? 100. : constants::KernPad;

    const int N_int_pts = (constants::KernPad < 0) ? 100000 : 10000;

    // Compute int_{LB}^{UB} ( kernel(r) * r^3 ) dr
    double alpha = 0., norm_fact = 0., r_loc = 0.;
    const double dr = (UB - LB) / N_int_pts;
    for (int II = 0; II < N_int_pts; ++II) {
        r_loc = II*dr;
        alpha     += M_PI * dr * ( kernel(r_loc, 1.) * pow(r_loc, 3) );
        norm_fact += M_PI * dr * ( kernel(r_loc, 1.) *     r_loc     );
    }

    return 0.5 * alpha / norm_fact;
}
