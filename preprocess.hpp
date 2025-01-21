#ifndef PREPROCESS_HPP
#define PREPROCESS_HPP 1

#include <stdio.h>
#include <stdlib.h>
#include "ALGLIB/linalg.h"
#include <mpi.h>
#include <vector>

#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/OrderingMethods>

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
        const std::vector<short int> &mask);

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
        const std::vector<short int> &mask,
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
        const std::vector<short int> &mask,
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

void map_grid_to_grid(
        const dataset & source_data,
        dataset & target_data,
        std::vector<std::string> vars_to_map,
        const MPI_Comm comm = MPI_COMM_WORLD
        );

/*!
 * \brief Class to store relevant variables for Helmholtz projections
 *
 */
class HelmholtzDataClass {

    public:

        //
        //// Variables
        //

        std::vector<short int> all_land_neighbours;
        size_t num_coastal, num_all_land;
        std::vector<size_t> pt_maps_to;
        size_t Ncol, Nrow, Npts_mapped;
        std::vector<size_t> num_mapped_before_col, num_mapped_before_row, num_land_before;
        std::vector<size_t> island_reps;
        std::map< size_t, std::vector<size_t> > coastal_boundaries;

        // Eigen-related
        Eigen::SparseMatrix<double> LHS;
        Eigen::LeastSquaresConjugateGradient< Eigen::SparseMatrix<double> > solver;
        Eigen::SparseQR< Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > direct_solver;
        Eigen::VectorXd RHS, x0, soln;
        std::vector<double> RHS_vector;

        //
        bool weight_err = true, use_vel = true, use_vort_div = true, collapse_land = false;
        double Tikhov = 1, tolerance = 0;
        unsigned int iteration_max, iterations_per_cycle;
        double stagnation_tolerance = 1e-30;


        // Convergence tracking
        std::vector<double> vel_2_errors, vort_2_errors, div_2_errors,
                            vel_inf_errors, vort_inf_errors, div_inf_errors,
                            vel_2_norms, vort_2_norms, div_2_norms,
                            vel_inf_norms, vort_inf_norms, div_inf_norms;


        //
        //// Functions
        //

        // Constructor
        HelmholtzDataClass();

        // Clear
        void clear();

        // Identify which points have land-only neighbours
        void IdentifyLandlockedPoints( const dataset & data, 
                                       const bool second_order_adjacency = false );

        // Create land-collapsing map [i.e. map continguous land to single point]
        void CreateLandCollapsingMap( const dataset & data );

        // Build the LHS and RHS part of the problem
        void Build_LHS( const dataset & data );
        void Build_RHS( const dataset & data );
        void VerifyRowNorms();

        void Build_Scalar_LHS_and_RHS( const dataset & data );

        // Set the seed for the solver
        void Set_Seed( const dataset & data );

        // Set the solver
        void InitializeSolver();

        // Apply the solver with guess x0
        void Solve();

        // Apply direct solver to sparse system (no guess)
        void DirectSolve();

        // Extract Psi/Phi from the solver grid onto the physical grid
        void Extract_PsiPhi( dataset & data );

        // Extract the projected vels
        void Extract_ProjectedVars( dataset & data );

        // Extract Land-filled Scalar
        void Extract_LandfilledScalar( dataset & data );

        // Extract the residuals
        void Extract_Residuals( dataset & data );

        // Compute projection errors
        void ComputeProjectionErrors( const dataset & data );

        // Determine convergence
        bool IsConverged() const;
        bool StagnationTestScalar() const;
        bool StagnationTestVector() const;

};

void Helmholtz_Solver(
        const std::string output_fname,
        dataset & source_data,
        const double rel_tol,
        const unsigned int max_iters,
        const unsigned int iters_per_batch,
        const bool weight_err,
        const bool use_mask,
        const bool use_vel,
        const bool use_vort_div,
        const bool collapse_land, 
        const int num_refinements,
        const MPI_Comm comm = MPI_COMM_WORLD
        );

void Helmholtz_Solver_wSPH(
        const std::string output_fname,
        dataset & source_data,
        const double rel_tol,
        const unsigned int max_iters,
        const unsigned int iters_per_batch,
        const bool weight_err,
        const bool use_mask,
        const bool use_vel,
        const bool use_vort_div,
        const bool collapse_land, 
        const int num_refinements,
        const MPI_Comm comm = MPI_COMM_WORLD
        );

void Helmholtz_Solver_Diffusion(
        const std::string output_fname,
        dataset & source_data,
        const double rel_tol,
        const unsigned int max_iters,
        const unsigned int iters_per_batch,
        const bool weight_err,
        const bool use_mask,
        const bool use_vel,
        const bool use_vort_div,
        const bool collapse_land, 
        const int num_refinements,
        const double CFL,
        const double hyper_visc,
        const MPI_Comm comm = MPI_COMM_WORLD
        );

void Scalar_Solver(
        const std::string output_fname,
        dataset & source_data,
        const double rel_tol,
        const unsigned int max_iters,
        const unsigned int iters_per_batch,
        const bool weight_err,
        const bool use_mask,
        const int num_refinements,
        const MPI_Comm comm = MPI_COMM_WORLD
        );

void initialize_coarsened_grid(
        dataset & coarsened_grid,
        const dataset & reference_grid,
        const int Nlat_coarse,
        const int Nlon_coarse
        );

void BuildPolyhedralGrid(
        dataset & polyhedral_grid,
        const size_t target_num_points 
        );

void SphericalHarmonicSolver(
        std::vector<double> & reconstructed_field,
        const dataset & source_data,
        const std::vector<double> & Laplaclian_of_field
        );



void Apply_Helmholtz_Projection_Eigen(
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
        const double filter_scale = -1,
        const MPI_Comm comm = MPI_COMM_WORLD
        );

void Apply_LLC_Helmholtz_Projection_Eigen(
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

void Apply_LLC_Helmholtz_Projection_Eigen_vels(
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

void Apply_LLC_Helmholtz_Projection_Eigen_vels_DeltaLand(
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

void Apply_LLC_Helmholtz_Projection_Eigen_PsiPhi_DeltaLand(
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

void Apply_LLC_Helmholtz_Projection_ALGLIB_PsiPhi_DeltaLand(
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

void Apply_LLC_Helmholtz_Projection_Eigen_both_DeltaLand(
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

void Apply_LLC_Helmholtz_Projection_Eigen_both(
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

void Apply_LLC_Helmholtz_Projection_AMGCL(
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
        const std::vector<short int> & mask
    );

void potential_vel_from_F(  
        std::vector<double> & vel_lon,
        std::vector<double> & vel_lat,
        const std::vector<double> & F,
        const dataset & source_data,
        const std::vector<short int> & mask
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
        const std::vector<short int> & mask,
        const std::vector<double> * seed = NULL
        );

void scalar_laplacian(  
        std::vector<double> & Lap_scalar,
        const std::vector<double> & scalar,
        const dataset & source_data,
        const std::vector<short int> & mask
    );


void toroidal_sparse_Lap(
        alglib::sparsematrix & Lap,
        const dataset & source_data,
        const int Itime,
        const int Idepth,
        const std::vector<short int> & mask,
        const bool area_weight = false,
        const size_t row_skip = 0,
        const size_t column_skip = 0
        );

void sparse_vel_from_PsiPhi(
        alglib::sparsematrix & LHS_matr,
        const dataset & source_data,
        const int Itime,
        const int Idepth,
        const std::vector<short int> & mask,
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
        const std::vector<short int> & mask
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
        const std::vector<short int> & mask
    );

void Extract_Beta_Geos_Vel(
        std::vector<double> & u_beta,
        std::vector<double> & v_beta,
        const std::vector<double> & ssh,
        const std::vector<short int> & mask,
        dataset & source_data,
        const double rel_tol,
        const int max_iters,
        const MPI_Comm comm = MPI_COMM_WORLD
        );

#endif
