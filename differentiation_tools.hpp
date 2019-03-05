#ifndef DIFFERENTIATION_TOOLS_HPP
#define DIFFERENTIATION_TOOLS_HPP 1

#include <stdio.h>
#include <stdlib.h>
#include <vector>

/*!
 * \file
 * \brief Collection of all differentiation-related functions.
 */

/*!
 * \brief Assign appropriate differentiation vector
 *
 * This function produces the differentiation vector necessary 
 *   for computing the value of a derivative at a single point.
 *
 * Currently assumes fourth-order finite differencing on a five-point uniform stencil.
 *
 * The grid need not be centred (to account for coastlines).
 */
void differentiation_vector(std::vector<double> & diff_array, const double delta, const int index);

/*!
 * \brief Computes latitudinal derivative at a specific point.
 */
double latitude_derivative_at_point(
        const std::vector<double> & field, const std::vector<double> & latitude,
        const int Itime, const int Idepth, const int Ilat, const int Ilon,
        const int Ntime, const int Ndepth, const int Nlat, const int Nlon,
        const std::vector<double> & mask);

/*!
 * \brief Computes longitudinal derivative at a specific point.
 */
double longitude_derivative_at_point(
        const std::vector<double> & field, const std::vector<double> & longitude,
        const int Itime, const int Idepth, const int Ilat, const int Ilon,
        const int Ntime, const int Ndepth, const int Nlat, const int Nlon,
        const std::vector<double> & mask);

/*!
 * \brief Computes Cartesian x derivative at a specific point.
 *
 * Computation is done via chain rule on spherical coordinate derivatives.
 *
 * Calls latitude_derivative_at_point() and longitude_derivative_at_point()
 */
double x_derivative_at_point(
        const std::vector<double> & field, 
        const std::vector<double> & latitude, 
        const std::vector<double> & longitude,
        const int Itime, const int Idepth, const int Ilat, const int Ilon,
        const int Ntime, const int Ndepth, const int Nlat, const int Nlon,
        const std::vector<double> & mask);

/*!
 * \brief Computes Cartesian y derivative at a specific point.
 *
 * Computation is done via chain rule on spherical coordinate derivatives.
 *
 * Calls latitude_derivative_at_point() and longitude_derivative_at_point()
 */
double y_derivative_at_point(
        const std::vector<double> & field, 
        const std::vector<double> & latitude, 
        const std::vector<double> & longitude,
        const int Itime, const int Idepth, const int Ilat, const int Ilon,
        const int Ntime, const int Ndepth, const int Nlat, const int Nlon,
        const std::vector<double> & mask);

/*!
 * \brief Computes Cartesian z derivative at a specific point.
 *
 * Computation is done via chain rule on spherical coordinate derivatives.
 *
 * Calls latitude_derivative_at_point() and longitude_derivative_at_point()
 */
double z_derivative_at_point(
        const std::vector<double> & field, 
        const std::vector<double> & latitude, 
        const std::vector<double> & longitude,
        const int Itime, const int Idepth, const int Ilat, const int Ilon,
        const int Ntime, const int Ndepth, const int Nlat, const int Nlon,
        const std::vector<double> & mask);

#endif
