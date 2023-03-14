#ifndef PREPROCESS_HPP
#define PREPROCESS_HPP 1

#include <stdio.h>
#include <stdlib.h>
#include "ALGLIB/linalg.h"
#include <mpi.h>
#include <vector>

/*!
 * \file
 * \brief Various pre-processing routines related to coarse-graining.
 * - Apply ALGLIB interpolation routines to fill in masked areas
 * - Apply ALGLIB least-squares solvers to apply toroidal projections
 */


/*!
 *  \addtogroup InterpolationRoutines
 *  @{
 * \brief Functions directly pertaining to toroidal projection.
 */

/*!
 * @ingroup InterpolationRoutines
 */
void interpolate_over_land(
        std::vector<double> &interp_field,
        const std::vector<double> &field,
        const std::vector<double> &time,
        const std::vector<double> &depth,
        const std::vector<double> &latitude,
        const std::vector<double> &longitude,
        const std::vector<bool> &mask);

/*!
 * @ingroup InterpolationRoutines
 */
void interpolate_over_land_from_coast(
        std::vector<double> &interp_field,
        const std::vector<double> &field,
        const int                 nlayers,
        const std::vector<double> &time,
        const std::vector<double> &depth,
        const std::vector<double> &latitude,
        const std::vector<double> &longitude,
        const std::vector<bool> &mask,
        const std::vector<int>    &myCounts,
        const MPI_Comm comm = MPI_COMM_WORLD
        );

/*!
 * @ingroup InterpolationRoutines
 */
void get_coast(
        std::vector<double> &lon_coast,
        std::vector<double> &lat_coast,
        std::vector<double> &field_coast,
        const std::vector<double> &lon_full,
        const std::vector<double> &lat_full,
        const std::vector<double> &field_full,
        const std::vector<bool> &mask,
        const int Itime,
        const int Idepth,
        const int Ntime,
        const int Ndepth,
        const int Nlat,
        const int Nlon);

void depth_integrate(
        std::vector<double> & depth_integral,
        const std::vector<double> & field_to_integrate,
        const dataset & source_data,
        const MPI_Comm comm = MPI_COMM_WORLD
        );

void Apply_Helmholtz_Projection(
        const std::string output_fname,
        dataset & source_data,
        const std::vector<double> & seed_tor,
        const std::vector<double> & seed_pot,
        const bool single_seed,
        const double rel_tol,
        const int max_iters,
        const bool weight_err,
        const bool use_mask,
        const double Tikhov_Laplace,
        const MPI_Comm comm = MPI_COMM_WORLD
        );

void Apply_LLC_Helmholtz_Projection(
        const std::string output_fname,
        dataset & source_data,
        const std::vector<double> & seed_tor,
        const std::vector<double> & seed_pot,
        const bool single_seed,
        const double rel_tol,
        const int max_iters,
        const bool weight_err,
        const bool use_mask,
        const double Tikhov_Laplace,
        const MPI_Comm comm = MPI_COMM_WORLD
        );

void Apply_Helmholtz_Projection_uiuj(
        const std::string output_fname,
        dataset & source_data,
        const std::vector<double> & seed_v_r,
        const std::vector<double> & seed_v_lon,
        const std::vector<double> & seed_v_lat,
        const bool single_seed,
        const double Tikhov_Lambda,
        const double Tikhov_Laplace,
        const double rel_tol,
        const int max_iters,
        const bool weight_err,
        const bool use_mask,
        const MPI_Comm comm = MPI_COMM_WORLD
        );

void Apply_Helmholtz_Projection_SymTensor(
        const std::string output_fname,
        dataset & source_data,
        const std::vector<double> & seed_v_r,
        const std::vector<double> & seed_v_lon,
        const std::vector<double> & seed_v_lat,
        const bool single_seed,
        const double rel_tol,
        const int max_iters,
        const bool weight_err,
        const bool use_mask,
        const MPI_Comm comm = MPI_COMM_WORLD
        );

/*!
 * \brief Computes the (toroidal) velocity corresponding to field F.
 * @ingroup ToroidalProjection
 *
 * Specifically, \f$ \hat{e}_r \times \frac{1}{r}\nabla_H F \f$
 *
 * @param[in,out]   vel_lon,vel_lat         Where to store the toroidal velocities
 * @param[in]       F                       Field from which to compute the velocities
 * @param[in]       longitude,latitude      Grid vectors (1D)
 * @param[in]       Ntime,Ndepth,Nlat,Nlon  Dimension sizes
 * @param[in]       mask                    Array to distinguish land/water
 *
 */
void toroidal_vel_from_F(  
        std::vector<double> & vel_lon,
        std::vector<double> & vel_lat,
        const std::vector<double> & F,
        const dataset & source_data,
        const std::vector<bool> & mask
    );

void potential_vel_from_F(  
        std::vector<double> & vel_lon,
        std::vector<double> & vel_lat,
        const std::vector<double> & F,
        const dataset & source_data,
        const std::vector<bool> & mask
    );

void uiuj_from_Helmholtz(  
        std::vector<double> & ulon_ulon,
        std::vector<double> & ulon_ulat,
        std::vector<double> & ulat_ulat,
        const std::vector<double> & v_r,
        const std::vector<double> & Phi_v,
        const std::vector<double> & Psi_v,
        const dataset & source_data
    );


/*!
 * \brief Computes the curl term (RHS) of the projection operation.
 * @ingroup ToroidalProjection
 *
 * Computes \f$ \nabla_H \times \vec{u} \cdot \hat{e}_r \f$, which
 * corresponds to the RHS of the least-squares problem.
 *
 * Uses the differentation order specified in contants.hpp
 *
 * *seed*: I couldn't figure out how to provide a seed directly to the solver. Instead,
 * seeds are done in the following way. Call the seed \f$ x_0 \f$ and write \f$ x = x' + x_0 \f$.
 * If \f$ Ax=b\f$ then \f$Ax' = b - Ax_0 \f$.
 * The seed is applied in exactly this way to modify the RHS of the problem. 
 * The seed is then added back in afterwards.
 *
 * @param[in,out]   Lap                     Where to store the (sparse) differentiation matrix
 * @param[in]       longitude,latitude      Grid vectors (1D)
 * @param[in]       Itime,Idepth            Indicates current time/depth iteration  
 * @param[in]       Ntime,Ndepth,Nlat,Nlon  Dimension sizes
 * @param[in]       mask                    Array to distinguish land/water
 * @param[in]       seed                    (optional) seed for the solver
 *
 */
void toroidal_curl_u_dot_er(
        std::vector<double> & out_arr,
        const std::vector<double> & u_lon,
        const std::vector<double> & u_lat,
        const dataset & source_data,
        const std::vector<bool> & mask,
        const std::vector<double> * seed = NULL
        );


void toroidal_sparse_Lap(
        alglib::sparsematrix & Lap,
        const dataset & source_data,
        const int Itime,
        const int Idepth,
        const std::vector<bool> & mask,
        const bool area_weight = false,
        const size_t row_skip = 0,
        const size_t column_skip = 0
        );

void sparse_vel_from_PsiPhi(
        alglib::sparsematrix & LHS_matr,
        const dataset & source_data,
        const int Itime,
        const int Idepth,
        const std::vector<bool> & mask,
        const bool area_weight
        );


/*!
 * \brief This is just a helper to compute Lap(F). It's provided as an output for diagnostic purposes.
 * @ingroup ToroidalProjection
 *
 * @param[in,out]   out_arr                 Where to store the laplacian
 * @param[in]       F                       Field to differentiate
 * @param[in]       longitude,latitude      Grid vectors (1D)
 * @param[in]       Ntime,Ndepth,Nlat,Nlon  Dimension sizes
 * @param[in]       mask                    Array to distinguish land/water
 *
 */
void toroidal_Lap_F(
        std::vector<double> & out_arr,
        const std::vector<double> & F,
        const std::vector<double> & longitude,
        const std::vector<double> & latitude,
        const int Ntime,
        const int Ndepth,
        const int Nlat,
        const int Nlon,
        const std::vector<bool> & mask
        );


/*!
 * \brief This is just a helper to compute div(vel). It's provided as an output for diagnostic purposes.
 * @ingroup ToroidalProjection
 *
 * @param[in,out]   div                     Where to store the divergence
 * @param[in]       vel_lon,vel_lat         Velocity fields
 * @param[in]       longitude,latitude      Grid vectors (1D)
 * @param[in]       Ntime,Ndepth,Nlat,Nlon  Dimension sizes
 * @param[in]       mask                    Array to distinguish land/water
 *
 */
void toroidal_vel_div(  
        std::vector<double> & div,
        const std::vector<double> & vel_lon,
        const std::vector<double> & vel_lat,
        const dataset & source_data,
        const std::vector<bool> & mask
    );

void Extract_Beta_Geos_Vel(
        std::vector<double> & u_beta,
        std::vector<double> & v_beta,
        const std::vector<double> & ssh,
        const std::vector<bool> & mask,
        dataset & source_data,
        const double rel_tol,
        const int max_iters,
        const MPI_Comm comm = MPI_COMM_WORLD
        );

#endif
