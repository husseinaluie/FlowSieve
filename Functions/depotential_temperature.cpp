#include "../functions.hpp"

double depotential_temperature(
        const double p,     /**< [in] pressure at point  */
        const double theta  /**< [in] potential temperature at point */
        ) {

    // Following equation 1.154a from Vallis (p.36)

    // 0C = 273.15K

    const double T_0         = 283;      // K
    const double alpha_0     = 9.738e-4; // (m^3 / kg) 
    const double beta_T      = 1.67e-4;  // K^-1
    const double c_p0        = 3986.;    // J / (kg K)
    const double gamma_star  = 1.1e-8;   // Pa^-1
    const double beta_T_star = 1.00e-5;  // K^-2

    double theta_prime, c1, c2, T_prime, T;

    theta_prime = (theta + 273.15) - T_0;

    c1 = T_0 * alpha_0 * beta_T      / c_p0;
    c2 = T_0 * alpha_0 * beta_T_star / c_p0;
    T_prime =   c1 * p * (1 + 0.5 * gamma_star * p + c2 * p) 
              + theta_prime * (1 + c2 * p);

    T = T_prime + T_0 - 273.15;

    return T;
}
