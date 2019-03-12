#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP 1

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
    const double rho0 = 1e3;

    /*!
     * \param DiffOrd
     * \brief Differentiation order for finite differencing (currently must be 4)
     * @ingroup constants
     */
    const int DiffOrd = 6;

}
/*! @} End of Doxygen Groups*/

#endif
