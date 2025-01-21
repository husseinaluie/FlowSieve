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
        const MPI_Comm comm
        ) {

    int wRank, wSize;
    MPI_Comm_rank( comm, &wRank );
    MPI_Comm_size( comm, &wSize );

    // If we've passed the DO_TIMING flag, then create some timing vars
    Timing_Records timing_records;
    double clock_on;

    const size_t Npts = source_data.mask.size();

    // Fill in the land areas with zero velocity
    size_t index;
    #pragma omp parallel default(none) \
    shared( source_data ) \
    private( index ) \
    firstprivate( Npts )
    {
        #pragma omp for collapse(1) schedule(static)
        for (index = 0; index < Npts; index++) {
            if (not(source_data.mask[index])) {
                source_data.variables["u_lon"][index] = 0.;
                source_data.variables["u_lat"][index] = 0.;
            }
        }
    }

    #if DEBUG >= 2
    fprintf(stdout, "Initializing class objects for data storage.\n");
    #endif
    HelmholtzDataClass Helmholtz_data;
    Helmholtz_data.collapse_land = collapse_land;
    dataset solution_grid;
    dataset *coarsened_grid = new dataset();

    // Beginning the v-cycle [or, more acurrately, just a / cycle?]
    // Start at the coarsest resolution, solve, refine, repeat
    int Nlat_coarse = 0, Nlon_coarse = 0;
    size_t Npts_coarse = 0;
    for ( int refine_level = 0; refine_level >= 0; refine_level-- ) {

        // Otherwise we're on the original grid, so just use the source
        coarsened_grid = &source_data;
        Npts_coarse = Npts;

        coarsened_grid->variables.insert( std::pair< std::string, std::vector<double> >(
                    "vort", std::vector<double>(Npts_coarse, 0.) ) );
        coarsened_grid->variables.insert( std::pair< std::string, std::vector<double> >(
                    "div", std::vector<double>(Npts_coarse, 0.) ) );

        // Get vorticity and divergence
        #if DEBUG >= 2
        fprintf(stdout, "Computing vorticity on coarse grid.\n");
        #endif
        if (constants::DO_TIMING) { clock_on = MPI_Wtime(); }
        toroidal_curl_u_dot_er( 
                coarsened_grid->variables.at("vort"),
                coarsened_grid->variables.at("u_lon"),
                coarsened_grid->variables.at("u_lat"),
                *coarsened_grid,
                coarsened_grid->mask
                );

        #if DEBUG >= 2
        fprintf(stdout, "Computing divergence on coarse grid.\n");
        #endif
        toroidal_vel_div( 
                coarsened_grid->variables.at("div"),
                coarsened_grid->variables.at("u_lon"),
                coarsened_grid->variables.at("u_lat"),
                *coarsened_grid,
                coarsened_grid->mask
                );

        // Set vort and div to zero on land
        size_t Ipt;
        if ( use_vort_div ) {
            #pragma omp parallel default(none) \
            shared( coarsened_grid ) \
            private( Ipt ) \
            firstprivate( Npts_coarse )
            {
                #pragma omp for collapse(1) schedule(static)
                for ( Ipt = 0; Ipt < Npts_coarse; ++Ipt) {
                    if ( not( coarsened_grid->mask[Ipt] ) ) {
                        coarsened_grid->variables.at("vort")[Ipt] = 0.;
                        coarsened_grid->variables.at( "div")[Ipt] = 0.;
                    }
                }
            }
        }

        // Compute norms of velocities, vort, and div
        //   we'll use these later to normalize the components of
        //   the RHS
        // NOTE: This is a linear-algebra norm, not a physical norm
        //       [i.e. not weighted by space]
        double uo_norm = 0, vo_norm = 0, vort_norm = 0, div_norm = 0;
        if ( use_vort_div and use_vel ) {
            #pragma omp parallel default(none) \
            shared( coarsened_grid ) \
            private( Ipt ) \
            firstprivate( Npts_coarse ) \
            reduction( +:uo_norm,vo_norm,vort_norm,div_norm )
            {
                #pragma omp for collapse(1) schedule(static)
                for ( Ipt = 0; Ipt < Npts_coarse; ++Ipt) {
                    uo_norm   += pow( coarsened_grid->variables.at("u_lon")[Ipt], 2.);
                    vo_norm   += pow( coarsened_grid->variables.at("u_lat")[Ipt], 2.);
                    vort_norm += pow( coarsened_grid->variables.at( "vort")[Ipt], 2.);
                    div_norm  += pow( coarsened_grid->variables.at(  "div")[Ipt], 2.);
                }
            }
            uo_norm   = sqrt( uo_norm   / Npts_coarse );
            vo_norm   = sqrt( vo_norm   / Npts_coarse );
            vort_norm = sqrt( vort_norm / Npts_coarse );
            div_norm  = sqrt( div_norm  / Npts_coarse );
        }


        // Reset the Helmholtz data for this resolution
        Helmholtz_data.clear();
        Helmholtz_data.collapse_land = collapse_land;
        Helmholtz_data.use_vort_div = use_vort_div;
        Helmholtz_data.use_vel = use_vel;
        if ( use_vort_div and use_vel ) {
            Helmholtz_data.Tikhov = ( uo_norm + vo_norm ) / ( vort_norm + div_norm ); // scale by avg. of norms
        } else {
            Helmholtz_data.Tikhov = 1.;
        }
        Helmholtz_data.weight_err = weight_err;
        Helmholtz_data.tolerance = rel_tol;
        Helmholtz_data.iterations_per_cycle = iters_per_batch;
        Helmholtz_data.iteration_max = max_iters;

        if ( collapse_land ) {
            // Identify which points have land-only neighbours
            #if DEBUG >= 2
            fprintf(stdout, "Building map of land-locked points.\n");
            #endif
            if (constants::DO_TIMING) { clock_on = MPI_Wtime(); }
            Helmholtz_data.IdentifyLandlockedPoints( *coarsened_grid, use_vort_div );
            if (constants::DO_TIMING) { timing_records.add_to_record(MPI_Wtime() - clock_on, "IdentifyLandlockedPoints"); }

            // Create land-collapsing map [i.e. map contiguous land to single point]
            #if DEBUG >= 2
            fprintf(stdout, "Identifying contiguous land masses and creating map.\n");
            #endif
            if (constants::DO_TIMING) { clock_on = MPI_Wtime(); }
            Helmholtz_data.CreateLandCollapsingMap( *coarsened_grid );
            if (constants::DO_TIMING) { timing_records.add_to_record(MPI_Wtime() - clock_on, "CreateLandCollapsingMap"); }
        } else {
            // If we've turned off land collapsing, then just make some dummy variables
            Helmholtz_data.pt_maps_to.resize(Npts_coarse, 0);
            for ( size_t II = 0; II < Npts_coarse; II++ ) { Helmholtz_data.pt_maps_to[II] = II; }
            Helmholtz_data.all_land_neighbours.resize(Npts_coarse, false);
            Helmholtz_data.num_mapped_before_col.resize( Npts_coarse, 1 ); 
            Helmholtz_data.num_mapped_before_col[0] = 0;
            Helmholtz_data.num_mapped_before_row.resize( Npts_coarse, 0 );
            Helmholtz_data.Nrow = Npts_coarse;
            Helmholtz_data.Ncol = Npts_coarse-1;
        }

        // Build the LHS and RHS part of the problem
        #if DEBUG >= 2
        fprintf(stdout, "Building least-squares LHS and RHS\n");
        #endif

        if (constants::DO_TIMING) { clock_on = MPI_Wtime(); }
        Helmholtz_data.Build_LHS( *coarsened_grid );
        Helmholtz_data.Build_RHS( *coarsened_grid );
        Helmholtz_data.VerifyRowNorms();
        if (constants::DO_TIMING) { timing_records.add_to_record(MPI_Wtime() - clock_on, "Build_LHS_RHS"); }

        // Use spherical harmonics to set the seed
        if (constants::DO_TIMING) { clock_on = MPI_Wtime(); }
        coarsened_grid->variables.insert( std::pair< std::string, std::vector<double> >(
                    "Psi", std::vector<double>(Npts_coarse, 0.) ) );
        SphericalHarmonicSolver(  coarsened_grid->variables.at("Psi"),
                                 *coarsened_grid,
                                  coarsened_grid->variables.at("vort") );

        coarsened_grid->variables.insert( std::pair< std::string, std::vector<double> >(
                    "Phi", std::vector<double>(Npts_coarse, 0.) ) );
        SphericalHarmonicSolver(  coarsened_grid->variables.at("Phi"),
                                 *coarsened_grid,
                                  coarsened_grid->variables.at("div") );
        Helmholtz_data.Set_Seed( *coarsened_grid );
        if (constants::DO_TIMING) { timing_records.add_to_record(MPI_Wtime() - clock_on, "SphericalHarmonicSolver"); }

        //
        //// Output the seed / starting guess
        //

        // Now that solving is done, output the solution
        std::vector<std::string> vars_to_add = { 
            "projected_u_lon", "projected_u_lat", 
            "proj_uo_tor", "proj_uo_pot", "proj_vo_tor", "proj_vo_pot" ,
            "proj_vort", "proj_div"
        };
        if ( use_vort_div ) {
            vars_to_add.push_back( "projected_vort" );
            vars_to_add.push_back( "projected_div" );
        }
        for ( size_t Ivar = 0; Ivar < vars_to_add.size(); Ivar++ ) {
            coarsened_grid->variables.insert( std::pair< std::string, std::vector<double> >(
                        vars_to_add[Ivar], std::vector<double>(Npts_coarse, 0.) ) );
        }

        // Compute velocity components
        toroidal_vel_from_F( coarsened_grid->variables.at("proj_uo_tor"),
                             coarsened_grid->variables.at("proj_vo_tor"),
                             coarsened_grid->variables.at("Psi"),
                             *coarsened_grid,
                             coarsened_grid->mask );

        potential_vel_from_F( coarsened_grid->variables.at("proj_uo_pot"),
                              coarsened_grid->variables.at("proj_vo_pot"),
                              coarsened_grid->variables.at("Phi"),
                              *coarsened_grid,
                              coarsened_grid->mask );
        
        // Compute vorticity and divergence
        toroidal_curl_u_dot_er( 
                coarsened_grid->variables.at("proj_vort"),
                coarsened_grid->variables.at("proj_uo_tor"),
                coarsened_grid->variables.at("proj_vo_tor"),
               *coarsened_grid,
                coarsened_grid->mask
                );

        toroidal_vel_div( 
                coarsened_grid->variables.at("proj_div"),
                coarsened_grid->variables.at("proj_uo_pot"),
                coarsened_grid->variables.at("proj_vo_pot"),
               *coarsened_grid,
                coarsened_grid->mask
                );

        //  we'll do this at each resolution for tracking convergence etc.
        Helmholtz_data.Extract_ProjectedVars( *coarsened_grid );
        std::ostringstream seed_stream;
        seed_stream << "seed_R" << (refine_level) << ".nc";
        const std::string seed_filename = seed_stream.str();
        write_Helmholtz_output( seed_filename, Helmholtz_data, *coarsened_grid, *coarsened_grid );

        // Get the starting errors, for reference
        Helmholtz_data.ComputeProjectionErrors( *coarsened_grid );
        if (use_vort_div and use_vel) {
            fprintf( stdout,"\nSolver initializing with relative 2-norm errors of vels, vort, and div of %.2e, %.2e, and %.2e\n", 
                    Helmholtz_data.vel_2_errors.back()  / Helmholtz_data.vel_2_norms.back(),
                    Helmholtz_data.vort_2_errors.back() / Helmholtz_data.vort_2_norms.back(), 
                    Helmholtz_data.div_2_errors.back()  / Helmholtz_data.div_2_norms.back() 
                   );
        } else if (use_vort_div)  {
            fprintf( stdout,"\nSolver initializing with relative 2-norm errors of vort and div of %.2e, and %.2e\n", 
                    Helmholtz_data.vort_2_errors.back() / Helmholtz_data.vort_2_norms.back(), 
                    Helmholtz_data.div_2_errors.back()  / Helmholtz_data.div_2_norms.back() 
                   );
        } else {
            fprintf( stdout,"\nSolver initializing with relative 2-norm errors of vels of %.2e\n", 
                    Helmholtz_data.vel_2_errors.back()  / Helmholtz_data.vel_2_norms.back()
                   );
        }

        // Set the solver
        if (constants::DO_TIMING) { clock_on = MPI_Wtime(); }
        Helmholtz_data.InitializeSolver();
        if (constants::DO_TIMING) { timing_records.add_to_record(MPI_Wtime() - clock_on, "InitializeSolver"); }
        Helmholtz_data.solver.setMaxIterations(iters_per_batch);

        bool keep_solving = true;
        unsigned int total_iters = 0, iter_cycle = 0;
        #if (DEBUG >= 0) and (DEBUG < 2)
        int perc_base = 5;
        int perc = perc_base, perc_count = 0;
        #endif

        while ( keep_solving ) {

            // Apply solver for base_iters number of iterations, and then test for convergence
            // using b - A*x as the rhs
            if (constants::DO_TIMING) { clock_on = MPI_Wtime(); }
            Helmholtz_data.Solve();
            if (constants::DO_TIMING) { timing_records.add_to_record(MPI_Wtime() - clock_on, "Solving"); }

            // Increment counter
            iter_cycle++;

            // Compute how much the solution has changed
            bool is_stagnated = ( iter_cycle > 0 ) ? Helmholtz_data.StagnationTestVector() : false;
            if ( is_stagnated ) { 
                keep_solving = false; 
                fprintf( stdout,"  Halting at solver cycle %d. Solver has stagnated.\n", iter_cycle );
            }

            // Increment the seed
            //      note that x0 stores the accumulated solution across all solver iterations
            //Helmholtz_data.x0 = Helmholtz_data.soln + Helmholtz_data.x0;
            Helmholtz_data.x0 = Helmholtz_data.soln;// + Helmholtz_data.x0; // solveWithGuess side-steps this

            // Compute errors
            Helmholtz_data.ComputeProjectionErrors( *coarsened_grid );

            // If sufficient convergence, stop
            if ( Helmholtz_data.IsConverged() ) { 
                keep_solving = false;
                fprintf( stdout,"  Halting at solver cycle %d. Solver has converged to desired tolerance.\n", iter_cycle );
            }

            // If too many iterations, stop
            total_iters += iters_per_batch;
            if ( total_iters >= max_iters ) { 
                keep_solving = false;
                fprintf( stdout,"  Halting at solver cycle %d. Iteration limit reached.\n", iter_cycle );
            }

            #if DEBUG >= 2
            if ( keep_solving ) {
                if (use_vort_div and use_vel) {
                    fprintf( stdout,"  Solver cycle %d complete with internal error %.2e. Relative 2-norm errors of vels, vort, and div are %.2e, %.2e, and %.2e\n", 
                            iter_cycle,
                            Helmholtz_data.solver.error(),
                            Helmholtz_data.vel_2_errors.back()  / Helmholtz_data.vel_2_norms.back(),
                            Helmholtz_data.vort_2_errors.back() / Helmholtz_data.vort_2_norms.back(), 
                            Helmholtz_data.div_2_errors.back()  / Helmholtz_data.div_2_norms.back() 
                           );
                } else if (use_vort_div) {
                    fprintf( stdout,"  Solver cycle %d complete with internal error %.2e. Relative 2-norm errors of vort and div are %.2e and %.2e\n", 
                            iter_cycle,
                            Helmholtz_data.solver.error(),
                            Helmholtz_data.vort_2_errors.back() / Helmholtz_data.vort_2_norms.back(), 
                            Helmholtz_data.div_2_errors.back()  / Helmholtz_data.div_2_norms.back() 
                           );
                } else {
                    fprintf( stdout,"  Solver cycle %d complete with internal error %.2e. Relative 2-norm errors of vels is %.2e\n", 
                            iter_cycle,
                            Helmholtz_data.solver.error(),
                            Helmholtz_data.vel_2_errors.back()  / Helmholtz_data.vel_2_norms.back()
                           );
                }
            }
            #elif DEBUG >= 0
            if ( keep_solving and ( total_iters / (double)max_iters * 100 ) >= perc ) {
                perc_count++;
                if (perc_count % 5 == 0) { fprintf(stdout, "|"); }
                else                     { fprintf(stdout, "."); }
                fflush(stdout);
                perc += perc_base;
            }
            #endif
        }

        if (use_vort_div and use_vel) {
            fprintf( stdout,"\nSolver completed with relative 2-norm errors of vels, vort, and div are %.2e, %.2e, and %.2e\n", 
                    Helmholtz_data.vel_2_errors.back()  / Helmholtz_data.vel_2_norms.back(),
                    Helmholtz_data.vort_2_errors.back() / Helmholtz_data.vort_2_norms.back(), 
                    Helmholtz_data.div_2_errors.back()  / Helmholtz_data.div_2_norms.back() 
                   );
        } else if (use_vort_div) {
            fprintf( stdout,"\nSolver completed with relative 2-norm errors of vort and div are %.2e and %.2e\n", 
                    Helmholtz_data.vort_2_errors.back() / Helmholtz_data.vort_2_norms.back(), 
                    Helmholtz_data.div_2_errors.back()  / Helmholtz_data.div_2_norms.back() 
                   );
        } else {
            fprintf( stdout,"\nSolver completed with relative 2-norm errors of vels is %.2e\n", 
                    Helmholtz_data.vel_2_errors.back()  / Helmholtz_data.vel_2_norms.back()
                   );
        }

        if (constants::DO_TIMING) { clock_on = MPI_Wtime(); }
        // Bring Psi, Phi onto the physical grid [not land-mapped]
        Helmholtz_data.Extract_PsiPhi( *coarsened_grid );

        // Now that solving is done, output the solution
        vars_to_add.clear();
        vars_to_add = { 
            "projected_u_lon", "projected_u_lat", //"residual_u_lon", "residual_u_lat", 
            "projected_vort", "projected_div", //"residual_div", "residual_vort",
            "proj_uo_tor", "proj_uo_pot", "proj_vo_tor", "proj_vo_pot" ,
            "proj_vort", "proj_div"
        };
        for ( size_t Ivar = 0; Ivar < vars_to_add.size(); Ivar++ ) {
            coarsened_grid->variables.insert( std::pair< std::string, std::vector<double> >(
                        vars_to_add[Ivar], std::vector<double>(Npts_coarse, 0.) ) );
        }

        // Compute velocity components
        toroidal_vel_from_F( coarsened_grid->variables.at("proj_uo_tor"),
                             coarsened_grid->variables.at("proj_vo_tor"),
                             coarsened_grid->variables.at("Psi"),
                             *coarsened_grid,
                             coarsened_grid->mask );

        potential_vel_from_F( coarsened_grid->variables.at("proj_uo_pot"),
                              coarsened_grid->variables.at("proj_vo_pot"),
                              coarsened_grid->variables.at("Phi"),
                              *coarsened_grid,
                              coarsened_grid->mask );
        
        // Compute vorticity and divergence
        toroidal_curl_u_dot_er( 
                coarsened_grid->variables.at("proj_vort"),
                coarsened_grid->variables.at("proj_uo_tor"),
                coarsened_grid->variables.at("proj_vo_tor"),
               *coarsened_grid,
                coarsened_grid->mask
                );

        toroidal_vel_div( 
                coarsened_grid->variables.at("proj_div"),
                coarsened_grid->variables.at("proj_uo_pot"),
                coarsened_grid->variables.at("proj_vo_pot"),
               *coarsened_grid,
                coarsened_grid->mask
                );


        //  we'll do this at each resolution for tracking convergence etc.
        Helmholtz_data.Extract_ProjectedVars( *coarsened_grid );

        if ( not(use_vort_div) ) {
            coarsened_grid->variables.at("proj_vort") = coarsened_grid->variables.at("projected_vort");
            coarsened_grid->variables.at("proj_div") = coarsened_grid->variables.at("projected_div");
        }
    
        if ( not(use_vel) ) {
            for ( size_t Ipt = 0; Ipt < Npts_coarse; Ipt++ ) {
                coarsened_grid->variables.at("projected_u_lon")[Ipt] = coarsened_grid->variables.at("proj_uo_tor")[Ipt] 
                    + coarsened_grid->variables.at("proj_uo_pot")[Ipt];
                coarsened_grid->variables.at("projected_u_lat")[Ipt] = coarsened_grid->variables.at("proj_vo_tor")[Ipt] 
                    + coarsened_grid->variables.at("proj_vo_pot")[Ipt];
            }
        }

        std::ostringstream name_stream;
        name_stream << "projection_R" << (refine_level) << ".nc";
        const std::string output_filename = name_stream.str();
        write_Helmholtz_output( output_filename, Helmholtz_data, *coarsened_grid, *coarsened_grid );
        if (constants::DO_TIMING) { timing_records.add_to_record(MPI_Wtime() - clock_on, "ProcessingSolution"); }

        // If we're doing timings, then print out and reset values now
        if (constants::DO_TIMING) { 
            timing_records.print();
            timing_records.reset();
            fflush(stdout);
        }

    }
}
