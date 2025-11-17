#include "../constants.hpp"
#include "../functions.hpp"
#include "../netcdf_io.hpp"
#include "../preprocess.hpp"
#include "../differentiation_tools.hpp"
#include <algorithm>
#include <utility> // std::swap
#include <vector>
#include <deque>
#include <omp.h>
#include <math.h>
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>


void Scalar_Solver(
        const std::string output_fname,
        dataset & source_data,
        const double rel_tol,
        const unsigned int max_iters,
        const unsigned int iters_per_batch,
        const bool weight_err,
        const bool use_mask,
        const int num_refinements,
        const MPI_Comm comm
        ) {

    int wRank, wSize;
    MPI_Comm_rank( comm, &wRank );
    MPI_Comm_size( comm, &wSize );

    const size_t Npts = source_data.mask.size();
    size_t index, Ipt;

    // Get number of land points
    size_t num_land = 0;
    #pragma omp parallel default(none) \
    shared( source_data ) \
    private( index ) \
    firstprivate( Npts ) \
    reduction( +:num_land )
    {
        num_land = 0;
        #pragma omp for collapse(1) schedule(static)
        for (index = 0; index < Npts; index++) {
            if (not(source_data.mask[index])) {
                num_land += 1;
            }
        }
    }

    #if DEBUG >= 2
    fprintf(stdout, "Initializing class objects for data storage.\n");
    #endif
    HelmholtzDataClass Helmholtz_data;
    dataset solution_grid;
    dataset *coarsened_grid = new dataset();

    // Beginning the v-cycle [or, more acurrately, just a / cycle?]
    // Start at the coarsest resolution, solve, refine, repeat
    int Nlat_coarse = 0, Nlon_coarse = 0;
    size_t Npts_coarse = 0;
    for ( int refine_level = num_refinements; refine_level >= 0; refine_level-- ) {

        #if DEBUG >= 2
        fprintf(stdout, "\n\nBeginning refinement level %d.\n", refine_level);
        #endif

        // If refinement is turned on, do it
        if ( refine_level > 0 ) {
            Nlat_coarse = floor(sqrt( 0.5 * (Npts / pow(4., refine_level)) ));
            Nlon_coarse = 2 * Nlat_coarse;
            Npts_coarse = (size_t) Nlat_coarse * Nlon_coarse;

            #if DEBUG >= 2
            fprintf(stdout, "Refinement level has Nlat x Nlon of %d x %d.\n", Nlat_coarse, Nlon_coarse);
            #endif


            // Downsample the fields onto the coarse grid
            std::vector<std::string> vars_to_map = { "scalar" };

            initialize_coarsened_grid( *coarsened_grid, source_data, Nlat_coarse, Nlon_coarse );
            for ( size_t Ivar = 0; Ivar < vars_to_map.size(); Ivar++ ) {
                coarsened_grid->variables.insert( std::pair< std::string, std::vector<double> >(
                            vars_to_map[Ivar], std::vector<double>(Npts_coarse, 0.) ) );
            }

            #if DEBUG >= 2
            fprintf(stdout, "Mapping onto the coarse grid.\n");
            #endif
            map_grid_to_grid( source_data, *coarsened_grid, vars_to_map );

        } else {
            // Otherwise we're on the original grid, so just use the source
            coarsened_grid = &source_data;
            Npts_coarse = Npts;
            if (constants::GRID_TYPE == constants::GridType::MeshGrid) {
                Nlat_coarse = source_data.Nlat;
                Nlon_coarse = source_data.Nlon;
            }
        }

        coarsened_grid->variables.insert( std::pair< std::string, std::vector<double> >(
                    "land_filled_scalar", std::vector<double>(Npts_coarse, 0.) ) );

        // Now, up-sample the previous solutions onto the new grid
        if ( refine_level < num_refinements ) {
            #if DEBUG >= 2
            fprintf(stdout, "Up-sampling from the solution grid.\n");
            #endif
            std::vector<std::string> vars_to_map = { "land_filled_scalar" };
            map_grid_to_grid( solution_grid, *coarsened_grid, vars_to_map );
        }

        #if DEBUG >= 2
        fprintf(stdout, "(Re)setting up the solution grid.\n");
        #endif
        // And reset the solution grid for the new resolution
        if (constants::GRID_TYPE == constants::GridType::MeshGrid) {
            initialize_coarsened_grid( solution_grid, source_data, Nlat_coarse, Nlon_coarse );
        } else {
            initialize_coarsened_grid( solution_grid, source_data, 
                    refine_level == 0 ? 1 : Nlat_coarse, 
                    refine_level == 0 ? Npts_coarse : Nlon_coarse );
        }
        solution_grid.variables.insert( std::pair< std::string, std::vector<double> >(
                    "land_filled_scalar", std::vector<double>(Npts_coarse, 0.) ) );

        if ( refine_level == 0 ) {
            // If we're on the original grid, then enforce matching grids
            solution_grid.mask = source_data.mask;
            solution_grid.longitude = source_data.longitude;
            solution_grid.latitude = source_data.latitude;
        }

        // Reset the Helmholtz data for this resolution
        Helmholtz_data.clear();
        Helmholtz_data.weight_err = weight_err;
        Helmholtz_data.tolerance = rel_tol;
        Helmholtz_data.iterations_per_cycle = iters_per_batch;
        Helmholtz_data.iteration_max = max_iters;

        // Identify which points have land-only neighbours
        #if DEBUG >= 2
        fprintf(stdout, "Building map of land-locked points.\n");
        #endif
        Helmholtz_data.IdentifyLandlockedPoints( *coarsened_grid );

        // Create land-collapsing map [i.e. map contiguous land to single point]
        #if DEBUG >= 2
        fprintf(stdout, "Identifying contiguous land masses and creating map.\n");
        #endif
        Helmholtz_data.CreateLandCollapsingMap( *coarsened_grid );
        // Add one if last point is land
        //Helmholtz_data.Nrow = Helmholtz_data.num_land_before.back() + 1;
        Helmholtz_data.Nrow = Helmholtz_data.num_land_before.back() + 0;
        Helmholtz_data.Ncol = Helmholtz_data.Nrow;

        // Build the LHS and RHS part of the problem
        #if DEBUG >= 2
        fprintf(stdout, "Building least-squares LHS and RHS\n");
        #endif
        Helmholtz_data.Build_Scalar_LHS_and_RHS( *coarsened_grid );
        if ( refine_level < num_refinements ) {
            Helmholtz_data.x0.resize(   Helmholtz_data.Ncol );
            for ( size_t ii = 0; ii < Npts_coarse; ii++ ) {
                if ( coarsened_grid->mask[ii] ) { continue; } // skip water
                Ipt = Helmholtz_data.num_land_before[ii];
                Helmholtz_data.x0[Ipt] = coarsened_grid->variables.at("land_filled_scalar").at(ii);
            }
        } else {
            Helmholtz_data.x0.resize(   Helmholtz_data.Ncol );
            Helmholtz_data.soln.resize( Helmholtz_data.Ncol );
            for ( size_t ii = 0; ii < Helmholtz_data.Ncol; ii++ ) {
                Helmholtz_data.x0[ii] = 0.;
                Helmholtz_data.soln[ii] = 0.;
            }
        }

        // Set the solver
        Helmholtz_data.InitializeSolver();
        Helmholtz_data.solver.setMaxIterations(iters_per_batch);

        bool keep_solving = true;
        unsigned int total_iters = 0, iter_cycle = 0;
        while ( keep_solving ) {

            // Apply solver for base_iters number of iterations, and then test for convergence
            // using b - A*x as the rhs
            Helmholtz_data.Solve();

            // Increment counter
            iter_cycle++;

            // Compute how much the solution has changed
            bool is_stagnated = Helmholtz_data.StagnationTestScalar();
            if ( is_stagnated ) { 
                keep_solving = false; 
                fprintf( stdout,"  Halting at solver cycle %d. Solver has stagnated.\n", iter_cycle );
            }

            // Increment the seed
            //      note that x0 stores the accumulated solution across all solver iterations
            Helmholtz_data.x0 = Helmholtz_data.soln + Helmholtz_data.x0;

            // If too many iterations, stop
            total_iters += iters_per_batch;
            if ( total_iters >= max_iters ) { 
                keep_solving = false; 
                fprintf( stdout,"  Halting at solver cycle %d. Iteration limit reached.\n", iter_cycle );
            }

            if (keep_solving) {
                fprintf( stdout,"  Solver cycle %d complete.\n", iter_cycle );
            }
        }

        // Now that solving is done, output the solution
        //  we'll do this at each resolution for tracking convergence etc.
        solution_grid.variables.insert( std::pair< std::string, std::vector<double> >(
                    "scalar", std::vector<double>(Npts_coarse, 0.) ) );
        solution_grid.variables.at("scalar") = coarsened_grid->variables.at("scalar");
        Helmholtz_data.Extract_LandfilledScalar( solution_grid );
        std::ostringstream name_stream;
        name_stream << "projection_R" << (refine_level) << ".nc";
        const std::string output_filename = name_stream.str();
        write_HelmholtzScalar_output( output_filename, Helmholtz_data, 
                *coarsened_grid, solution_grid );

    }
}
