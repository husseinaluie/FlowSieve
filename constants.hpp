#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP 1

#ifndef DEBUG
    #define DEBUG 0
#endif

#ifndef CARTESIAN
    #define CARTESIAN false
#endif

#ifndef PERIODIC_X
    #define PERIODIC_X true
#endif

#ifndef PERIODIC_Y
    #define PERIODIC_Y false
#endif

#ifndef COMP_VORT
    #define COMP_VORT true
#endif

#ifndef COMP_TRANSFERS
    #define COMP_TRANSFERS true
#endif

#ifndef COMP_BC_TRANSFERS
    #define COMP_BC_TRANSFERS false
#endif

/*!
 * \file
 * \brief Provide namespace for global constants (physical and computational).
 *
 * This header file provides a namespace to give a consistent source of constants.
 *
 * Usage:
 * @code
 * #include "functions.hpp"
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
     * \brief (constant) acceleration due to gravity
     * @ingroup constants
     */
    const double g = 9.81;

    /*!
     * \param DiffOrd
     * \brief Differentiation order for finite differencing (currently must be 4)
     * @ingroup constants
     */
    const int DiffOrd = 6;

}
/*! @} End of Doxygen Groups*/

#endif
