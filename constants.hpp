#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP 1

#include <map>

// Default debug value
#ifndef DEBUG
    #define DEBUG 1
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
     * \param USE_HIGH_PRECISION_DISTANCE
     * \brief If spatial resolution is very high on sphere (less than ~50 metres or so), 
     * and high presision is needed in distance calculations, turn this on.
     *
     * The default spherical distance caluclation, spherical law of cosines, has floating point
     * errors for small distances (a couple metres or so). However, so long as your grid is very high resolution,
     * this shouldn't be an issue. This is particularly true if: you are also using a continuous kernel and the filtering scales
     * themselves also aren't very short (couple of metres or so).
     *
     * The test routine (Tests/distance_formulas.cpp) uses both methods on a specified grid and outputs the result.
     * If you think you might need the high precision scheme, test it out there first.
     *
     * @ingroup constants
     */
    const bool USE_HIGH_PRECISION_DISTANCE = false;

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
    const int DiffOrd = 4;

    /*!
     * \param fill_value
     * \brief Fill value used to indicate land values in output files.
     * @ingroup constants
     */
    const double fill_value = -1e8;

    /*!
     * \param fill_value_s
     * \brief Fill value used to indicate land values in output files (signed short)
     * @ingroup constants
     */
    const signed short fill_value_s = -32767;

    /*!
     * \param DEFORM_AROUND_LAND
     * \brief Boolean indicating whether or not the kernel should deform around land.
     *
     * If true, land is given zero weight (i.e. not factored into local average).
     * If false, land is included in local average (treated as zero velocity
     * @ingroup constants
     */
    const bool DEFORM_AROUND_LAND = false;

    /*!
     * \param FILTER_OVER_LAND
     * \brief Boolean to indicate whether or not land values should be filled in with coarse-grained results.
     *
     * If true, land values are filled in with the results of coarse-graining. 
     * Masks are removed when variables are loaded and filled in with constants
     * If false, land values are maintained as land, and masked out
     *
     * Note that this does NOT affect the shape of filtering kernels.
     * @ingroup constants
     */
    const bool FILTER_OVER_LAND = true;

    /*!
     * \param CARTESIAN
     * \brief Boolean indicating if the coordinate system is Cartesian. (If false, then spherical)
     * @ingroup constants
     */
    const bool CARTESIAN = false;

    /*!
     * \param PERIODIC_X
     * \brief Boolean indicating if the coordinate system is periodic in x / longitude.
     *
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
     * \param UNIFORM_LON_GRID
     * \brief Boolean indicating if the longitude grid is uniform
     * @ingroup constants
     */
    const bool UNIFORM_LON_GRID = true;

    /*!
     * \param UNIFORM_LAT_GRID
     * \brief Boolean indicating if the latitude grid is uniform
     * @ingroup constants
     */
    const bool UNIFORM_LAT_GRID = true;

    /*!
     * \param FULL_LON_SPAN
     * \brief Boolean indicating if the provided longitude grid spans the full periodic domain.
     *
     * i.e. is the first longitude point beside the last longitude point?
     *
     * This flag enables computational optimizations by allowing the kernel to simply be translated
     * through lon, instead of having to recompute for each longitude index.
     *
     * @ingroup constants
     */
    const bool FULL_LON_SPAN = true;

    /*!
     * \param COMP_VORT
     * \brief Boolean indicating if vorticity should be computed.
     * @ingroup constants
     */
    const bool COMP_VORT = true;

    /*!
     * \param COMP_TRANSFERS
     * \brief Boolean indicating if non-linear transfers (Pi) should be computed.
     * For coarse_grain.x, this is required to get fine (sub-filter) KE
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
     * \param DO_OKUBOWEISS_ANALYSIS
     * \brief Boolean to indicate if post-processing should also bin by Okubo-Weiss
     * @ingroup constants
     */
    const bool DO_OKUBOWEISS_ANALYSIS = false;

    /*!
     * \param MINIMAL_OUTPUT
     * \brief Boolean indicating if user wants a minimal output
     *
     * Removes: fine velocities, coarse/fine KE, divergences
     *
     * @ingroup constants
     */
    const bool MINIMAL_OUTPUT = true;

    /*! 
     * \param NO_FULL_OUTPUTS
     * \brief Indicates that no full fields should be produced.
     *   
     * Only results from APPLY_POSTPROCESS will be produced.
     *
     * If NO_FULL_OUTPUTS and not(APPLY_POSTPROCESS), the no outputs
     * at all will be produced, so the code will simply halt immediately.
     *
     * @ingroup constants
     */
    const bool NO_FULL_OUTPUTS = true;

    /*!
     * \param CAST_TO_SINGLE
     * \brief Boolean indicating if user wants to cast to float (single) output
     * (reduces output size by factor 2, but also reduces precision)
     *
     * @ingroup constants
     */
    const bool CAST_TO_SINGLE = true;

    /*!
     * \param CAST_TO_INT
     * \brief Boolean indicating if user wants to cast to int output
     * (further reduces output size by factor 2, but also reduces precision)
     *
     * @ingroup constants
     */
    const bool CAST_TO_INT = false;

    /*!
     * \param DO_TIMING
     * \brief Boolean indicating if we want to output internal timings
     * @ingroup constants
     */
    const bool DO_TIMING = true;

    /*!
     * \param APPLY_POSTPROCESS
     * \brief Boolean indicating whether or not the postprocess routines should be applied
     * @ingroup constants
     */
    const bool APPLY_POSTPROCESS = true;

    /*!
     * \param POSTPROCESS_DO_TIME_MEANS
     * \brief Boolean indicating whether or not the postprocess routines should include time means (i.e. spatial maps)
     * @ingroup constants
     */
    const bool POSTPROCESS_DO_TIME_MEANS = true;

    /*!
     * \param KERNEL_OPT
     * \brief Integer flag indicating the choice of kernel
     *
     *  0 = tophat
     *  1 = Hyper gaus  ( exp( -x^4 )                  )
     *  2 = Gaus        ( exp( -x^2 )                  )
     *  3 = sinc        ( sinc( pi * x )               )
     *  4 = tanh        ( 1 - tanh( (x - 1) / (0.1) )  )
     * @ingroup constants
     */
    const int KERNEL_OPT = 4;

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
        case 4: const double KernPad =  2.5;
    }
    */


    /*!
     * \param variable_descriptions
     * \brief A dictionary of variable descriptions to provide details in netcdf outputs
     *
     * @ingroup constants
     */
    const std::map< std::string, std::string > variable_descriptions = {
        { "coarse_u_r",         "Coarse-grained vertical/radial velocity." },
        { "coarse_u_lon",       "Coarse-grained zonal velocity." },
        { "coarse_u_lat",       "Coarse-grained meridional velocity." },
        { "coarse_KE",          "Kinetic energy of coarse-grained velocity" },
        { "fine_KE",            "Small-scale kinetic energy ( filter(KE(u)) - KE(filter(u)) )" },
        { "enstrophy",          "Enstrophy of coarse-grained velocity" },
        { "Pi",                 "Non-linear energy transfer from large-scales to small-scales" },
        { "Z",                  "Non-linear enstrophy transfer from large-scales to small-scales" },
        { "OkuboWeiss",         "Okubo-Weiss parameter ( positive -> strain dominated, negative -> vortex dominated )" },
        { "div_Jtransport",     "divergence of energy transport term" },
        { "coarse_vort_r",      "Radial (z) vorticity of coarse-grained velocity." },
        { "coarse_vel_div",     "Divergence of the coarse-grained velocity" }
    };

    /*!
     * \param variable_units
     * \brief A dictionary of variable units to provide details in netcdf outputs
     *
     * @ingroup constants
     */
    const std::map< std::string, std::string > variable_units = {
        { "coarse_u_r",         "m / s" },
        { "coarse_u_lon",       "m / s" },
        { "coarse_u_lat",       "m / s" },
        { "coarse_KE",          "J / (m^3)" },
        { "fine_KE",            "J / (m^3)" },
        { "enstrophy",          "J / (m^5)" },
        { "Pi",                 "Watt / (m^3)" },
        { "Z",                  "Watt / (m^5)" },
        { "OkuboWeiss",         "1 / (s^2)" },
        { "div_Jtransport",     "Watt / (m^3)" },
        { "coarse_vort_r",      "1 / s" },
        { "coarse_vel_div",     "1 / s" }
    };

    const std::string spatial_average_description       = "The lat/lon average computed over each defined region (see region dimension).";
    const std::string time_average_description          = "Time average over the entire provided dataset.";
    const std::string OkuboWeiss_average_description    = "Variable binned by Okubo-Weiss parameter (i.e. histogram). "
                                                          "Values in OkuboWeiss dimension indicate the lower bound of each bin.";

}
/*! @} End of Doxygen Groups*/

#endif
