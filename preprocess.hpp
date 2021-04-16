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
        const std::vector<double> &time,
        const std::vector<double> &depth,
        const std::vector<double> &latitude,
        const std::vector<double> &longitude,
        const std::vector<bool> &mask,
        const std::vector<int>    &myCounts
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

/*!
 *  \addtogroup ToroidalProjection
 *  @{
 * \brief Functions directly pertaining to toroidal projection.
 */

/*!
 * \brief Applies toroidal projection to the specified velocity field.
 * @ingroup ToroidalProjection
 *
 * Specifically, solves for (scalar field) F that minimizes error in 
 * \f$ \nabla_H^2 F = \nabla_H \times \vec{u} \cdot \hat{e}_r \f$
 *
 * This uses sparse differentiation matrices (of the order specified in constants.hpp) and the ALGLIB least-squares solver.
 *
 * *single_seed*: If true, then only one seed is provided and should be used as the seed for all different times. If false, then the provided seed incorporates a different seed for each time.
 *
 * @param[in]       output_fname                    Name for the output file
 * @param[in,out]   u_lon,u_lat                     velocity field to be projected
 * @param[in]       time,depth,latitude,longitude   dimension vectors (1D)
 * @param[in]       mask                            array to distinguish land/water
 * @param[in]       myCounts                        Local (to MPI process) dimension sizes
 * @param[in]       myStarts                        Vector indicating where the local (to MPI process) region fits in the whole
 * @param[in]       seed                            Seed for the least-squares solver
 * @param[in]       single_seed                     Indicates if a single seed is used - see notes
 * @param[in]       comm                            MPI communicator (default MPI_COMM_WORLD)
 *
 */
void Apply_Toroidal_Projection(
        const std::string output_fname,
        dataset & source_data,
        const std::vector<double> & seed,
        const bool single_seed,
        const double rel_tol,
        const int max_iters,
        const bool weight_err,
        const bool use_mask,
        const MPI_Comm comm = MPI_COMM_WORLD
        );

void Apply_Potential_Projection(
        const std::string output_fname,
        dataset & source_data,
        const std::vector<double> & seed,
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
        const std::vector<double> & longitude,
        const std::vector<double> & latitude,
        const int Ntime,
        const int Ndepth,
        const int Nlat,
        const int Nlon,
        const std::vector<bool> & mask
    );

void potential_vel_from_F(  
        std::vector<double> & vel_lon,
        std::vector<double> & vel_lat,
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
 * Then \f$ Ax=b \implies Ax' = b - Ax_0 \f$.
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
        const std::vector<double> & longitude,
        const std::vector<double> & latitude,
        const int Itime,
        const int Idepth,
        const int Ntime,
        const int Ndepth,
        const int Nlat,
        const int Nlon,
        const std::vector<bool> & mask,
        const std::vector<double> * seed = NULL
        );


/*!
 * \brief Builds the (sparse!) Laplacian differentiation operator.
 * @ingroup ToroidalProjection
 *
 * Uses the differentation order specified in contants.hpp
 *
 * Currently only handles spherical coordinates.
 *
 * @param[in,out]   Lap                     Where to store the (sparse) differentiation matrix
 * @param[in]       longitude,latitude      Grid vectors (1D)
 * @param[in]       Itime,Idepth            Current time-depth iteration
 * @param[in]       Ntime,Ndepth,Nlat,Nlon  Dimension sizes
 * @param[in]       mask                    Array to distinguish land/water
 * @param[in]       areas                   Array giving the size of each cell
 * @param[in]       area_weight             Bool indicating if the Laplacian should be weighted by cell-size (i.e. weight error by cell size). Default is false.
 *
 */
void toroidal_sparse_Lap(
        alglib::sparsematrix & Lap,
        const std::vector<double> & latitude,
        const std::vector<double> & longitude,
        const int Itime,
        const int Idepth,
        const int Ntime,
        const int Ndepth,
        const int Nlat,
        const int Nlon,
        const std::vector<bool> & mask,
        const std::vector<double> & areas,
        const bool area_weight = false
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
        const std::vector<double> & longitude,
        const std::vector<double> & latitude,
        const int Ntime,
        const int Ndepth,
        const int Nlat,
        const int Nlon,
        const std::vector<bool> & mask
    );

#endif
