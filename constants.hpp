#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP 1

// Default debug value
#ifndef DEBUG
    #define DEBUG 0
#endif

/*!
 * \file
 * \brief Provide namespace for global constants (physical and computational).
 *
 * This header file provides a namespace to give a consistent source of constants.
 *
 * Usage:
 * @code
 * #include "constants.hpp"
 * ...
 * double R_earth = constants::R_earth;
 * @endcode
 */

/*!
 *  \addtogroup constants
 *  @{
 * \brief Provide namespace for global constants (physical and computational).
 */

/*! 
 * \namespace consants 
 * \brief Provide namespace for global constants (physical and computational).
 */
namespace constants
{

    /*!
     * \param R_earth
     * \brief Mean radius of the Earth
     * @ingroup constants
     */
    const double R_earth = 6371e3;

    /*!
     * \param rho0
     * \brief Mean fluid density
     * @ingroup constants
     */
    const double rho0 = 1025;

    /*!
     * \param g
     * \brief (constant) acceleration due to gravity.
     * @ingroup constants
     */
    const double g = 9.81;

    /*!
     * \param DiffOrd
     * \brief Differentiation order for finite differencing (currently must be 2, 4, or 6).
     * @ingroup constants
     */
    const int DiffOrd = 6;

    /*!
     * \param fill_value
     * \brief Fill value used to indicate land values in output files.
     * @ingroup constants
     */
    const double fill_value = -32767;

    /*!
     * \param fill_value_s
     * \brief Fill value used to indicate land values in output files (signed short)
     * @ingroup constants
     */
    const signed short fill_value_s = -32767;

    /*!
     * \param CARTESIAN
     * \brief Boolean indicating if the coordinate system is Cartesian. (If false, then spherical)
     * @ingroup constants
     */
    const bool CARTESIAN = false;

    /*!
     * \param PERIODIC_X
     * \brief Boolean indicating if the coordinate system is periodic in x / longitude.
     * @ingroup constants
     */
    const bool PERIODIC_X = true;

    /*!
     * \param PERIODIC_Y
     * \brief Boolean indicating if the coordinate system is periodic in y / latitute.
     * @ingroup constants
     */
    const bool PERIODIC_Y = false;

    /*!
     * \param UNIFORM_LAT_GRID
     * \brief Boolean indicating if the latitude grid is uniform
     * @ingroup constants
     */
    const bool UNIFORM_LAT_GRID = true;

    /*!
     * \param COMP_VORT
     * \brief Boolean indicating if vorticity should be computed.
     * @ingroup constants
     */
    const bool COMP_VORT = false;

    /*!
     * \param COMP_TRANSFERS
     * \brief Boolean indicating if non-linear transfers (Pi) should be computed.
     * @ingroup constants
     */
    const bool COMP_TRANSFERS = true;

    /*!
     * \param COMP_BC_TRANSFERS
     * \brief Boolean indicating if baroclinic transfers (Lambda^m) should be computed.
     * @ingroup constants
     */
    const bool COMP_BC_TRANSFERS = false;

    /*!
     * \param MINIMAL_OUTPUT
     * \brief Boolean indicating if user wants a minimal output
     *
     * Removes: fine velocities, coarse/fine KE, divergences, transport
     *
     * @ingroup constants
     */
    const bool MINIMAL_OUTPUT = true;

    /*!
     * \param CAST_TO_INT
     * \brief Boolean indicating if user wants to cast to int output
     * (reduces output size by factor 2, but also reduces precision)
     *
     * @ingroup constants
     */
    const bool CAST_TO_INT = true;

    /*!
     * \param DO_TIMING
     * \brief Boolean indicating if we want to output internal timings
     * @ingroup constants
     */
    const bool DO_TIMING = false;

    /*!
     * \param KERNEL_OPT
     * \brief Integer flag indicating the choice of kernel
     *
     *  0 = tophat
     *  1 = Hyper gaus (exp(-x^4))
     *  2 = Gaus  (exp(-x^2))
     *  3 = sinc (sharp-spectral)
     * @ingroup constants
     */
    const int KERNEL_OPT = 1;

    /*!
     * \param KernPad
     * \brief Scale factor for kernel search radius 
     *
     * Filter integral applied in circle of radius (filt_scale/2) * KernPad
     *
     * @ingroup constants
     */
    const double KernPad = 2.5;
    /*
    switch (KERNEL_OPT) {
        case 0: const double KernPad =  1.1;
        case 1: const double KernPad =  2.5;  // exp(-2.5^4) ~1e-17
        case 2: const double KernPad =  5.;   // exp(-5^2)   ~1e-11
        case 3: const double KernPad = -1.;
    }
    */

}
/*! @} End of Doxygen Groups*/

#endif
