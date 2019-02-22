#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP 1

/*
 * This header file provides a namespace to give a consistent source of constants.
 *
 * Usage:
 * #include "functions.hpp"
 * ...
 * double R_earth = constants::R_earth;
 *
 */

namespace constants
{
    // Mean radius of the Earth
    const double R_earth = 6371e3;

    // Mean fluid density
    const double rho0 = 1e3;
}

#endif
