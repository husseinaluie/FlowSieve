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
 * This function implicitly assumes a uniform grid.
 *
 * For first derivatives, can use 2nd, 4th, and 6th order convergence.
 *
 * For second derivatives, can use 2nd and 4th order convergence
 *
 * The grid need not be centred (to account for coastlines).
 *
 * @param[in,out]   diff_array      vector into which to store the differentiation coefficients
 * @param[in]       delta           the (constant) grid spacing
 * @param[in]       index           index for where in the stencil you want the derivative
 * @param[in]       order_of_deriv  order of the derivative (default is first derivative)
 * @param[in]       diff_ord        convergence order (default is specified in constants.hpp)
 *
 */
void differentiation_vector(
        std::vector<double> & diff_array, 
        const double delta, 
        const int index,
        const int order_of_deriv = 1,
        const int diff_ord = constants::DiffOrd);



/*!
 * \brief Assign appropriate differentiation vector
 *
 * This function produces the differentiation vector necessary 
 *   for computing the value of a derivative at a single point.
 *
 * The grid need not the uniform.
 *
 * For first derivatives, 2nd order convergence can be used.
 *
 * Second derivatives not yet implemented.
 *
 * The grid need not be centred (to account for coastlines).
 *
 * @param[in,out]   diff_array  vector into which to store the differentiation coefficients
 * @param[in]       grid        vector giving the grid on which the derivative is taken
 * @param[in]       Iref        integer giving the index for the point at which you want the derivative
 * @param[in]       LB          integer giving the lower bound for the differentiation stencil
 * @param[in]       UB          integer giving the upper bound for the differentiation stencil
 * @param[in]       diff_ord    convergence order (default is specified in constants.hpp)
 *
 */
void non_uniform_diff_vector(
        std::vector<double> & diff_array,
        const std::vector<double> & grid,
        const int Iref,
        const int LB,
        const int UB,
        const int diff_ord = constants::DiffOrd);


/*!
 * \brief Computes latitudinal derivative at a specific point.
 *
 * @param[in,out]   deriv_vals              vector of pointers to return values
 * @param[in]       fields                  vector of pointers to fields to differentiate
 * @param[in]       grid                    (1D) latitude  or longitude array
 * @param[in]       dim                     "lat" or "lon": dimensional along which to differentate
 * @param[in]       Itime,Idepth,Ilat,Ilon  Indices for the target point in space
 * @param[in]       Ntime,Ndepth,Nlat,Nlon  Sizes of the dimensions
 * @param[in]       mask                    array to distinguish land/water cells 
 * @param[in]       order_of_deriv          order of the derivative (1st deriv, 2nd deriv, etc)
 * @param[in]       diff_ord                convergence order (default is specified in constants.hpp)
 *
 */
void spher_derivative_at_point(
        const std::vector<double*> & deriv_vals,
        const std::vector<const std::vector<double>*> & fields,
        const std::vector<double> & grid,
        const std::string & dim,
        const int Itime, const int Idepth, const int Ilat, const int Ilon,
        const int Ntime, const int Ndepth, const int Nlat, const int Nlon,
        const std::vector<bool> & mask,
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
 * @param[in,out]   x_deriv_vals            vector of pointers to return x-deriv values
 * @param[in,out]   y_deriv_vals            vector of pointers to return y-deriv values
 * @param[in,out]   z_deriv_vals            vector of pointers to return z-deriv values
 * @param[in]       fields                  vector of pointers to fields to differentiate
 * @param[in]       latitude                (1D) latitude array
 * @param[in]       longitude               (1D) longitude array
 * @param[in]       Itime,Idepth,Ilat,Ilon  Indices for the target point in space
 * @param[in]       Ntime,Ndepth,Nlat,Nlon  Sizes of the dimensions
 * @param[in]       mask                    array to distinguish land/water cells 
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
        const std::vector<bool> & mask,
        const int order_of_deriv = 1,
        const int diff_ord = constants::DiffOrd);



/*!
 * \brief Function to return coefficients for differentiation
 * at a target point in space
 *
 * This is used for the Toroidal projection, which requires
 * building a Laplacian matrix.
 *
 * @param[in,out]   diff_vector             vector into which to store the coefficients
 * @param[in,out]   LB_ret                  index of the lower bound of the stencil (returned)
 * @param[in]       grid                    grid on which to differentiate
 * @param[in]       dim                     name of the differentiation dimension ("lon", "lat")
 * @param[in]       Itime,Idepth,Ilat,Ilon  Indices for the target point in space
 * @param[in]       Ntime,Ndepth,Nlat,Nlon  Sizes of the dimensions
 * @param[in]       mask                    mask to distinguish land/water cells
 * @param[in]       order_of_deriv          order of the derivative
 * @param[in]       diff_ord                convergence order of the derivative (default from constants.hpp)
 *
 */
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
        const std::vector<bool> & mask,
        const int order_of_deriv,
        const int diff_ord = constants::DiffOrd
        );

#endif
