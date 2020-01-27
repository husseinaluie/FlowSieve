#ifndef DIFFERENTIATION_TOOLS_HPP
#define DIFFERENTIATION_TOOLS_HPP 1

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <string>
#include "constants.hpp"

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
void differentiation_vector(
        std::vector<double> & diff_array, 
        const double delta, 
        const int index,
        const int order_of_deriv = 1,
        const int diff_ord = constants::DiffOrd);


void non_uniform_diff_vector(
        std::vector<double> & diff_array,
        const std::vector<double> & grid,
        const int Iref,
        const int LB,
        const int UB,
        const int diff_ord);

/*!
 * \brief Computes latitudinal derivative at a specific point.
 *
 * \param vector of pointers to return values
 * \param vector of pointers to fields to differentiate
 * \param (1D) latitude  or longitude array
 * \param "lat" or "lon": dimensional along which to differentate
 * \param time index at which to differentiate
 * \param depth index at which to differentiate
 * \param latitude index at which to differentiate
 * \param longitude at which to differentiate
 * \param size of time dimension
 * \param size of depth dimension
 * \param size of latitude dimension
 * \param size of longitude dimension
 * \param (2D) array to distinguish land/water cells 
 */
void spher_derivative_at_point(
        const std::vector<double*> & deriv_vals,
        const std::vector<const std::vector<double>*> & fields,
        const std::vector<double> & grid,
        const std::string & dim,
        const int Itime, const int Idepth, const int Ilat, const int Ilon,
        const int Ntime, const int Ndepth, const int Nlat, const int Nlon,
        const std::vector<double> & mask,
        const int order_of_deriv = 1,
        const int diff_ord = constants::DiffOrd);


/*!
 * \brief Computes all Cartesian derivative at a specific point.
 *
 * Computation is done via chain rule on spherical coordinate derivatives.
 * Since both lat- and lon- derivatives are taken, the cost increase to 
 * compute all three Cartesian derivatives in comparison to computing a 
 * single Cartesian derivative is negligible.
 *
 * Calls latitude_derivative_at_point() and longitude_derivative_at_point()
 *
 * \param vector of pointers to return x-deriv values
 * \param vector of pointers to return y-deriv values
 * \param vector of pointers to return z-deriv values
 * \param vector of pointers to fields to differentiate
 * \param (1D) latitude array
 * \param (1D) longitude array
 * \param time index at which to differentiate
 * \param depth index at which to differentiate
 * \param latitude index at which to diffentiate
 * \param longitude at which to differentiate
 * \param size of time dimension
 * \param size of depth dimension
 * \param size of latitude dimension
 * \param size of longitude dimension
 * \param (2D) array to distinguish land/water cells 
 */
void Cart_derivatives_at_point(
        const std::vector<double*> & x_deriv_vals,
        const std::vector<double*> & y_deriv_vals,
        const std::vector<double*> & z_deriv_vals,
        const std::vector<const std::vector<double>*> & fields,
        const std::vector<double> & latitude, 
        const std::vector<double> & longitude,
        const int Itime, const int Idepth, const int Ilat, const int Ilon,
        const int Ntime, const int Ndepth, const int Nlat, const int Nlon,
        const std::vector<double> & mask,
        const int order_of_deriv = 1,
        const int diff_ord = constants::DiffOrd);

// Return the differentiation vector
void get_diff_vector(
        std::vector<double> & diff_vector,
        int & LB_ret,
        const std::vector<double> & grid,
        const std::string & dim,
        const int Itime,
        const int Idepth,
        const int Ilat,
        const int Ilon,
        const int Ntime,
        const int Ndepth,
        const int Nlat,
        const int Nlon,
        const std::vector<double> & mask,
        const int order_of_deriv,
        const int diff_ord
        );

#endif
