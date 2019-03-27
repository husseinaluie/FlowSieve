#include <math.h>
#include "../functions.hpp"

double equation_of_state(
        const double T, /**< [in] temperature (degrees C) */
        const double S, /**< [in] salinity (PSU) */
        const double p  /**< [in] pressure (bars) */
        ) {


    // UNESCO 1981
    //   Valid for  0 < T < 40  (Celcius)
    //              0 < S < 42  (PSU)
    //   https://link.springer.com/content/pdf/bbm%3A978-3-319-18908-6%2F1.pdf

    // First, define a whole whack of constants

    const double a0 =  999.842594;
    const double a1 =    6.793953e-2;
    const double a2 =   -9.095290e-3;
    const double a3 =    1.001685e-4;
    const double a4 =   -1.120083e-6;
    const double a5 =    6.536332e-9;

    const double b0 =  8.2449e-1;
    const double b1 = -4.0899e-3;
    const double b2 =  7.6438e-5;
    const double b3 = -8.2467e-7;
    const double b4 =  5.3875e-9;

    const double c0 = -5.7246e-3;
    const double c1 =  1.0227e-4;
    const double c2 = -1.6546e-6;

    const double d0 =  4.8314e-4;

    const double e0 =  19652.210000;
    const double e1 =    148.420600;
    const double e2 =     -2.327105;
    const double e3 =      1.360477e-2;
    const double e4 =     -5.155288e-5;

    const double f0 =  54.674600;
    const double f1 =  -0.603459;
    const double f2 =   1.099870e-2;
    const double f3 =  -6.167000e-5;

    const double g0 =  7.9440e-2;
    const double g1 =  1.6483e-2;
    const double g2 = -5.3009e-4;

    const double h0 =  3.23990;
    const double h1 =  1.43714e-3;
    const double h2 =  1.16092e-4;
    const double h3 = -5.77905e-7;

    const double i0 =  2.28380e-3;
    const double i1 = -1.09810e-5;
    const double i2 = -1.60780e-6;
    
    const double j0 =  1.91075e-4;

    const double k0 =  8.50935e-5;
    const double k1 = -6.12293e-6;
    const double k2 =  5.27870e-8;

    const double m0 = -9.9348e-7;
    const double m1 =  2.0816e-8;
    const double m2 =  9.1697e-10;


    // Compute SMOW (Standard Mean Ocean Water)
    const double SMOW =   a0 
                        + a1 * T 
                        + a2 * pow(T, 2)
                        + a3 * pow(T, 3)
                        + a4 * pow(T, 4)
                        + a5 * pow(T, 5);


    // Compute density at normal atmospheric pressure (p=0)
    double B1 =   b0
                + b1 * T
                + b2 * pow(T,2)
                + b3 * pow(T,3)
                + b4 * pow(T,4);

    double C1 =   c0
                + c1 * T
                + c2 * pow(T,2);


    const double rho_ST0 = 
        SMOW + B1 * S + C1 * pow(S,1.5) + d0 * pow(S,2);


    // Determine compressible module at p = 0
    const double Kw =   e0
                      + e1 * T
                      + e2 * pow(T,2)
                      + e3 * pow(T,3)
                      + e4 * pow(T,4);

    const double F1 =   f0
                      + f1 * T
                      + f2 * pow(T,2)
                      + f3 * pow(T,3);

    const double G1 =   g0
                      + g1 * T
                      + g2 * pow(T,2);

    double K_ST0 = Kw + F1 * S + G1 * pow(S,1.5);


    // Determine compressible module of the sea water 
    const double Aw =   h0
                      + h1 * T
                      + h2 * pow(T,2)
                      + h3 * pow(T,3);

    const double A1 =   Aw
                      + (i0 + i1 * T + i2 * pow(T,2)) * S
                      + j0 * pow(S, 1.5);

    const double Bw =   k0
                      + k1 * T
                      + k2 * pow(T,2);

    const double B2 =   Bw
                      + (m0 + m1 * T + m2 * pow(T,2)) * S;

    const double K_STrho = K_ST0 + A1 * p + B2 * pow(p,2);

    // Compute the density
    double rho;
    rho = rho_ST0 / (1 - p / K_STrho );

    return rho;
}


