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

double my_sech( const double &x ) {
    if ( fabs(x) > 100 ) { return 0.; }
    return 1. / cosh(x);
}

double dsincdx(const double &x) {
    double s;

    // to avoid some numerical issues, use series expansion for small x
    if (fabs(x) < 1e-9) { s = - x/6.; }
    else                { s = ( x*cos(x) -  sin(x)) / my_pow<2>(x); }

    return s;
}

double tanh_approx( const double x ) {

    double P0, P1, P2, P3, P4, P5, P6,
           Q0, Q1, Q2, Q3, Q4, Q5, Q6;

    if ( std::fabs(x) < 4 ) {
        P0 = 9.999999974526363467e-01;
        P1 = -1.599443322746878782e-01;
        P2 = 3.595127484720697653e-02;
        P3 = -7.424073964463595576e-03;
        P4 = -1.778547808392130693e-03;
        P5 = -4.007298748504722896e-05;
        P6 = -9.141379148092099172e-08;

        Q0 = 1.;
        Q1 = 1.733889779940582243e-01;
        Q2 = -3.958566605905809027e-02;
        Q3 = 1.023033368843799008e-02;
        Q4 = -5.602254840994548395e-03;
        Q5 = -3.366473730445655385e-04;
        Q6 = -2.824768502461549868e-06;
    } else if ( std::fabs(x) < 7 ) {
        P0 = 3.449407438866221830e-04;
        P1 = 1.826548470619919224e-03;
        P2 = 5.229269672913945850e-03;
        P3 = 1.158050293906958328e-03;
        P4 = 4.112533925346342767e-03;
        P5 = 1.635148123074126525e-04;
        P6 = 3.832501200165625704e-07;

        Q0 = 1.;
        Q1 = -8.289969959998851360e-04;
        Q2 = -3.387462347792928076e-03;
        Q3 = -4.684906547614595997e-03;
        Q4 = 6.837164721711731785e-03;
        Q5 = 1.124417323940959295e-03;
        Q6 = 1.235244647821191148e-05;
    } else if ( std::fabs(x) < 10 ) {
        P0 = 7.059619037217820540e-09;
        P1 = 1.645788529795822215e-07;
        P2 = 2.675108090246749099e-06;
        P3 = 2.007672402802097752e-05;
        P4 = -1.573553261269088378e-05;
        P5 = -3.672550758982300500e-07;
        P6 = -4.080115189633750529e-10;

        Q0 = 1.;
        Q1 = -3.523136566620521056e-08;
        Q2 = -7.046421480757668838e-07;
        Q3 = -8.499830267785125977e-06;
        Q4 = -2.450700574089997782e-05;
        Q5 = -3.450225273652721158e-06;
        Q6 = -1.937966035912260994e-08;
    } else {
        return ( x > 0 ) ? 1 : -1;
    }

    const double z = x*x;
    const double P = ( ( ( ( ( P6 * z + P5 ) * z + P4 ) * z + P3 ) * z + P2 ) * z + P1 ) * z + P0;
    const double Q = ( ( ( ( ( Q6 * z + Q5 ) * z + Q4 ) * z + Q3 ) * z + Q2 ) * z + Q1 ) * z + Q0;

    return x * P / Q;

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
    const double D        = ( scale > 0 ) ? ( dist / ( scale / 2. ) )   : ( dist == 0 ) ? 1. : 0.;
    const double dDdell   = ( ( scale > 0 ) and (deriv_order > 0) ) ? ( -2 * dist / my_pow<2>(scale) ) : 0.,
                 d2Ddell2 = ( ( scale > 0 ) and (deriv_order > 1) ) ? (  4 * dist / my_pow<3>(scale) ) : 0.;

    const double arg = (D - 1) / 0.1;

    switch (constants::KERNEL_OPT) {
        case constants::KernelType::TopHat: 
                kern = D < 1 ? 1. : 0;
                // tophat has ill-defined ell-derivatives
                break;
        case constants::KernelType::HyperGaussian: 
                if ( deriv_order == 0 ) {
                    kern = exp( -my_pow<4>(D) );
                } else if (deriv_order == 1 ) {
                    kern = exp( -my_pow<4>(D) ) * (-4*my_pow<3>(D)) * dDdell;
                } else if (deriv_order == 2 ) {
                    kern = (4*D*D) * exp( -my_pow<4>(D) ) * ( 
                            (4*my_pow<4>(D) - 3) * my_pow<2>(dDdell) - D * d2Ddell2 ) ;
                }
                break;
        case constants::KernelType::Gaussian: 
                if ( deriv_order == 0 ) {
                    kern = exp( -my_pow<2>(D) );
                } else if (deriv_order == 1 ) {
                    kern = exp( -my_pow<2>(D) ) * (-2*D) * dDdell;
                } else if (deriv_order == 2 ) {
                    kern = exp( -my_pow<2>(D) ) * ( 
                            (4*D*D-2) * my_pow<2>(dDdell) - 2 * D * d2Ddell2 ) ;
                }
                break;
        case constants::KernelType::JohnsonGaussian: 
                if ( deriv_order == 0 ) {
                    kern = exp( -my_pow<2>(D)/8 );
                } else if (deriv_order == 1 ) {
                    kern = exp( -my_pow<2>(D) / 8 ) * (-2*D/8) * dDdell;
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
                    /*
                    if ( arg >= 10 ) { kern = 0.; }
                    else if ( arg <= -10 ) { kern = 1.; }
                    else { kern = 0.5 * (1 - tanh( arg )); }
                    */
                    kern = 0.5 * ( 1 - tanh_approx( arg ) );
                } else if (deriv_order == 1 ) {
                    kern = -5 * my_pow<2>(my_sech(arg)) * dDdell;
                } else if (deriv_order == 2 ) {
                    kern = 5 * my_pow<2>(my_sech(arg)) * ( 20 * tanh_approx(arg) * my_pow<2>(dDdell) - d2Ddell2 );
                }
                break;
        case constants::KernelType::HighOrder: 
                const double arg2 = (D - 1) / 0.5;
                const double c2 = 0.21534029041474162;
                if ( deriv_order == 0 ) {
                    kern = 0.5 * (1 - tanh( arg )) - c2 * exp( -my_pow<2>(arg2) );
                } else if (deriv_order == 1 ) {
                    kern = -5 * my_pow<2>(1./cosh(arg)) * dDdell
                        + c2 * 4 * arg2 * exp( -my_pow<2>(arg2) ) * dDdell;
                } else if (deriv_order == 2 ) {
                    kern = 5 * my_pow<2>(1./cosh(arg)) * ( 20 * tanh(arg) * my_pow<2>(dDdell) - d2Ddell2 )
                        + c2 * 4 * exp( -my_pow<2>(arg2) ) * ( 2 * (1 - 2*my_pow<2>(arg2)) * my_pow<2>(dDdell) 
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
