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

double dsincdx(const double &x) {
    double s;

    // to avoid some numerical issues, use series expansion for small x
    if (fabs(x) < 1e-9) { s = - x/6.; }
    else                { s = ( x*cos(x) -  sin(x)) / pow(x,2); }

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
        const double scale,
        const int deriv_order
        ) {
   
    double kern = 0.;
    const double D        = ( scale > 0 ) ? ( dist / ( scale / 2. ) )   : ( dist == 0 ) ? 1. : 0.,
                 dDdell   = ( scale > 0 ) ? ( -2 * dist / pow(scale,2) ) : 0.,
                 d2Ddell2 = ( scale > 0 ) ? (  4 * dist / pow(scale,3) ) : 0.;

    const double arg = (D - 1) / 0.1;

    switch (constants::KERNEL_OPT) {
        case constants::KernelType::TopHat: 
                kern = D < 1 ? 1. : 0;
                // tophat has ill-defined ell-derivatives
                break;
        case constants::KernelType::HyperGaussian: 
                if ( deriv_order == 0 ) {
                    kern = exp( -pow( D, 4) );
                } else if (deriv_order == 1 ) {
                    kern = exp( -pow( D, 4) ) * (-4*pow(D,3)) * dDdell;
                } else if (deriv_order == 2 ) {
                    kern = (4*D*D) * exp( -pow( D, 4) ) * ( 
                            (4*pow(D,4) - 3) * pow(dDdell,2) - D * d2Ddell2 ) ;
                }
                break;
        case constants::KernelType::Gaussian: 
                if ( deriv_order == 0 ) {
                    kern = exp( -pow( D, 2) );
                } else if (deriv_order == 1 ) {
                    kern = exp( -pow( D, 2) ) * (-2*D) * dDdell;
                } else if (deriv_order == 2 ) {
                    kern = exp( -pow( D, 2) ) * ( 
                            (4*D*D-2) * pow(dDdell,2) - 2 * D * d2Ddell2 ) ;
                }
                break;
        case constants::KernelType::JohnsonGaussian: 
                if ( deriv_order == 0 ) {
                    kern = exp( -pow( D, 2)/8 );
                } else if (deriv_order == 1 ) {
                    kern = exp( -pow( D, 2)/8 ) * (-2*D/8) * dDdell;
                } else if (deriv_order == 2 ) {
                    kern = 0.;
                }
                break;
        case constants::KernelType::Sinc: 
                // likely to be discontinued
                if ( deriv_order == 0 ) {
                    kern = sinc( M_PI * D );
                } else if (deriv_order == 1 ) {
                    kern = dsincdx( M_PI * D ) * M_PI * dDdell;
                } else if (deriv_order == 2 ) {
                }
                break;
        case constants::KernelType::SmoothHat: 
                if ( deriv_order == 0 ) {
                    kern = 0.5 * (1 - tanh( arg ));
                } else if (deriv_order == 1 ) {
                    kern = -5 * pow(1./cosh(arg), 2) * dDdell;
                } else if (deriv_order == 2 ) {
                    kern = 5 * pow(1./cosh(arg), 2) * ( 20 * tanh(arg) * pow(dDdell,2) - d2Ddell2 );
                }
                break;
        case constants::KernelType::HighOrder: 
                const double arg2 = (D - 1) / 0.5;
                const double c2 = 0.21534029041474162;
                if ( deriv_order == 0 ) {
                    kern = 0.5 * (1 - tanh( arg )) - c2 * exp( -pow(arg2,2.) );
                } else if (deriv_order == 1 ) {
                    kern = -5 * pow(1./cosh(arg), 2) * dDdell
                        + c2 * 4 * arg2 * exp( -pow(arg2,2.) ) * dDdell;
                } else if (deriv_order == 2 ) {
                    kern = 5 * pow(1./cosh(arg), 2) * ( 20 * tanh(arg) * pow(dDdell,2) - d2Ddell2 )
                        + c2 * 4 * exp( -pow(arg2,2.) ) * ( 2 * (1 - 2*pow(arg2,2)) * pow( dDdell, 2 ) 
                                                            + arg2 * d2Ddell2 
                                                          );
                }
                break;
    }

    #if DEBUG >= 6
    fprintf(stdout, "Kernel(dist=%.4g, scale=%.4g) = %.4g\n",
            dist, scale, kern);
    #endif

    return kern;
}
