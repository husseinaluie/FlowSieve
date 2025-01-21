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

void ReconstructCoastalInterior( 
        HelmholtzDataClass & Helmholtz_data, 
        dataset & data
        ){

    const size_t num_islands = Helmholtz_data.Ncol - Helmholtz_data.Nrow + 1;
    const size_t Npts = data.mask.size();

    std::vector<double> coastal_Psi_sums( num_islands, 0. ), 
                        coastal_Phi_sums( num_islands, 0. ),
                        coastal_areas( num_islands, 0. );

    std::map< size_t, size_t > island_indices;

    size_t Ipt, curr_island;
    for ( curr_island = 0; curr_island < num_islands; curr_island++ ) {
        size_t island_rep = Helmholtz_data.island_reps[curr_island];

        island_indices.insert( std::pair< size_t, size_t >( island_rep, curr_island ) );

        size_t Npts_coast = Helmholtz_data.coastal_boundaries[island_rep].size();
        for ( Ipt = 0; Ipt < Npts_coast; Ipt++ ) {
            size_t index = Helmholtz_data.coastal_boundaries[island_rep][Ipt];
            double dA = data.areas[index];
            coastal_Psi_sums[curr_island] += data.variables.at("Psi")[index] * dA;
            coastal_Phi_sums[curr_island] += data.variables.at("Phi")[index] * dA;
            coastal_areas[curr_island] += dA;
        }
    }


    // Now loop over space and fill in any land points
    // with their appropriate coastal average
    #pragma omp parallel default(none) \
    shared( Helmholtz_data, data, coastal_Psi_sums, coastal_Phi_sums, coastal_areas, island_indices ) \
    private( Ipt ) firstprivate( Npts )
    {
        #pragma omp for collapse(1) schedule(guided)
        for ( Ipt = 0; Ipt < Npts; Ipt++ ) {
            if ( not( Helmholtz_data.all_land_neighbours[Ipt] ) ) { continue; } 
            size_t island_rep = Helmholtz_data.pt_maps_to[Ipt];

            size_t island_index = island_indices[island_rep];

            data.variables.at("Psi")[Ipt] = coastal_Psi_sums[island_index] / coastal_areas[island_index];
            data.variables.at("Phi")[Ipt] = coastal_Phi_sums[island_index] / coastal_areas[island_index];
        }
    }

}

void IdentifyCoastalBoundaries( 
        HelmholtzDataClass & Helmholtz_data, 
        const dataset * data
        ){

    const size_t num_islands = Helmholtz_data.Ncol - Helmholtz_data.Nrow + 1;
    const size_t Npts = data->mask.size();
    Helmholtz_data.island_reps.resize(num_islands);

    size_t Ipt, curr_island = 0;
    for ( Ipt = 0; Ipt < Npts; Ipt++ ) {
        if ( (Helmholtz_data.all_land_neighbours[Ipt]) and ( Helmholtz_data.pt_maps_to[Ipt] == Ipt ) ) {
            Helmholtz_data.island_reps[curr_island] = Ipt;
            curr_island++;
            if ( curr_island > num_islands ) {
                throw std::runtime_error("Failed to identify coastal boundaries. To many representatives.");
            }
        }
    }

    for ( curr_island = 0; curr_island < num_islands; curr_island++ ) {

        // Get the island representative
        Ipt = Helmholtz_data.island_reps[curr_island];

        // And set up the boundary list with an empty vector. We'll add to it as we find more.
        Helmholtz_data.coastal_boundaries.insert( std::pair< size_t, std::vector<size_t> >( Ipt, std::vector<size_t>(0) ) );


        std::vector<bool>   was_rejected(Npts, false), 
                            was_accepted(Npts, false), 
                            planned_for_testing(Npts, false);
        // intentionally using the secretly-a-bitset vector<bool>. bitset itself
        // doesn't allow dynamic sizing

        // Next, seed the 'points to test' with the adjacent points of the rep
        std::deque<size_t> points_to_test;
        for (size_t II = 0; II < data->num_neighbours; II++ ) {
            size_t neighbour_index = data->adjacency_indices[Ipt][II];
            points_to_test.push_back( neighbour_index );
            planned_for_testing[ neighbour_index ] = true;
        }

        // So long as we still have points to test, keep testing!
        while ( points_to_test.size() > 0 ) {

            // Pull out the most-recently-added point, and remove it from the 'to test' list,
            // since we're testing it now.
            // Since we're pulling out the most-recently-added, that effectively makes this a
            // depth-first search to build the kernel.
            size_t Jpt = points_to_test.front();
            points_to_test.pop_front();
            planned_for_testing[ Jpt ] = false;

            if ( data->mask[Jpt] ) { continue; } // somehow we hit water without going through the coast?
                                                 // mostly just a sanity-skip

            if ( not( Helmholtz_data.all_land_neighbours[Jpt] ) ) {
                // We've found a coastal point, so add it to the islands 'boundary'
                Helmholtz_data.coastal_boundaries[Ipt].push_back(Jpt);
                was_accepted[Jpt] = true;
            } else {
                // Otherwise, we're still in the island, so keep searching through the adjacency
                was_rejected[Jpt] = true;

                for ( size_t II = 0; II < data->num_neighbours; II++ ) {
                    size_t neighbour_index = data->adjacency_indices[Jpt][II];

                    // but first check if that neighbour has been rejected already
                    if ( was_rejected[neighbour_index] ) { continue; }

                    // then check if that neighbour is already accepted
                    if ( was_accepted[neighbour_index] ) { continue; }

                    // then check if that neighbour is already on the search list
                    if ( planned_for_testing[neighbour_index] ) { continue; }

                    // if it's not on any of those list already
                    // then add it to the 'points to test'
                    points_to_test.push_back( neighbour_index );
                    planned_for_testing[neighbour_index] = true;

                }
            }
        }
    }
}

void get_norms(
        double & vel_2_err,
        double & vort_2_err,
        double & div_2_err,
        double & vel_2_norm,
        double & vort_2_norm,
        double & div_2_norm,
        double & vel_inf_err,
        double & vort_inf_err,
        double & div_inf_err,
        double & vel_inf_norm,
        double & vort_inf_norm,
        double & div_inf_norm,
        double & vort_2_ener,
        double & div_2_ener,
        const dataset * coarsened_grid
        ) {

    double total_area = 0;
    size_t Ipt;
    bool use_vel = true, use_vort_div = true;
    const size_t & Npts_coarse = coarsened_grid->mask.size();

    #pragma omp parallel default(none) \
    shared( coarsened_grid ) \
    private( Ipt ) \
    firstprivate( Npts_coarse, use_vel, use_vort_div ) \
    reduction( +:total_area,vel_2_err,vort_2_err,div_2_err,vel_2_norm,vort_2_norm,div_2_norm ) \
    reduction( +:vort_2_ener,div_2_ener )\
    reduction( max:vel_inf_err,vort_inf_err,div_inf_err,vel_inf_norm,vort_inf_norm,div_inf_norm )
    {
        #pragma omp for collapse(1) schedule(static)
        for ( Ipt = 0; Ipt < Npts_coarse; Ipt++ ) {

            double area = coarsened_grid->areas[Ipt];
            double incr;

            total_area += area;

            if ( use_vel ) {
                // Velocity Errors
                incr = sqrt(
                        pow(
                               coarsened_grid->variables.at("proj_uo_tor")[Ipt]
                             + coarsened_grid->variables.at("proj_uo_pot")[Ipt]
                             - coarsened_grid->variables.at("u_lon")[Ipt]
                            , 2)
                        +  pow(
                               coarsened_grid->variables.at("proj_vo_tor")[Ipt]
                             + coarsened_grid->variables.at("proj_vo_pot")[Ipt]
                             - coarsened_grid->variables.at("u_lat")[Ipt]
                            , 2)
                        );
                vel_inf_err = std::fmax( incr, vel_inf_err );
                vel_2_err += area * pow( incr, 2 ); 
            }

            if ( use_vort_div ) {

                // Vorticity errors
                incr = std::fabs(
                          coarsened_grid->variables.at("proj_vort")[Ipt]
                        - coarsened_grid->variables.at("vort")[Ipt]
                        );
                vort_inf_err = std::fmax( incr, vort_inf_err );
                vort_2_err += area * pow( incr, 2 );

                incr = std::fabs(
                          coarsened_grid->variables.at("proj_vort")[Ipt]
                        );
                vort_2_ener += area * pow( incr, 2 );

                // Divergence errors
                incr = std::fabs(
                          coarsened_grid->variables.at("proj_div")[Ipt]
                        - coarsened_grid->variables.at("div")[Ipt]
                        );
                div_inf_err = std::fmax( incr, div_inf_err );
                div_2_err += area * pow( incr, 2 );
                incr = std::fabs(
                          coarsened_grid->variables.at("proj_div")[Ipt]
                        );
                div_2_ener += area * pow( incr, 2 );
            }

            // Reference norms
            if ( use_vel ) {
                incr = sqrt(   pow(  coarsened_grid->variables.at("u_lon")[Ipt], 2)
                             + pow(  coarsened_grid->variables.at("u_lat")[Ipt], 2)
                        );
                vel_inf_norm = std::fmax( incr, vel_inf_norm );
                vel_2_norm += area * pow( incr, 2 ); 
            }

            if ( use_vort_div ) {
                incr = std::fabs( coarsened_grid->variables.at("vort")[Ipt] );
                vort_inf_norm = std::fmax( incr, vort_inf_norm );
                vort_2_norm += area * pow( incr, 2 );

                incr = std::fabs( coarsened_grid->variables.at("div")[Ipt] );
                div_inf_norm = std::fmax( incr, div_inf_norm );
                div_2_norm += area * pow( incr, 2 );
            }
        }
    }

    vel_2_err  = sqrt( vel_2_err  / total_area );
    vort_2_err = sqrt( vort_2_err / total_area );
    div_2_err  = sqrt( div_2_err  / total_area );

    vort_2_ener = sqrt( vort_2_ener / total_area );
    div_2_ener  = sqrt( div_2_ener  / total_area );

    vel_2_norm  = sqrt( vel_2_norm  / total_area );
    vort_2_norm = sqrt( vort_2_norm / total_area );
    div_2_norm  = sqrt( div_2_norm  / total_area );
}


void filter_scalars(
       dataset * coarsened_grid,
       dataset & solution_grid,
       std::vector<std::string> & vars_to_filter,
       const double & filter_scale,
       const double & alpha = 1.,
       const bool & only_pole = false
        ) {

    const size_t Npts_coarse = solution_grid.mask.size();
    const size_t Nvars = vars_to_filter.size();
    size_t Ipt;

    // Filter at ~grid-scale for stability
    #pragma omp parallel default(none) \
    shared( coarsened_grid, solution_grid, vars_to_filter ) \
    private( Ipt ) \
    firstprivate( Npts_coarse, filter_scale, Nvars, only_pole )
    {
        #pragma omp for collapse(1) schedule(guided)
        for ( Ipt = 0; Ipt < Npts_coarse; Ipt++) {
   
            std::vector<bool>   was_rejected(Npts_coarse, false), 
                                was_accepted(Npts_coarse, false), 
                                planned_for_testing(Npts_coarse, false);
            // intentionally using the secretly-a-bitset vector<bool>. bitset itself
            // doesn't allow dynamic sizing
  
            std::vector<double> value_sums( Nvars, 0. );
 
            double target_lat = solution_grid.latitude[Ipt];
            double target_lon = solution_grid.longitude[Ipt];

            if ( only_pole and (std::fabs(target_lat) * 180. / M_PI < 80.) ) {
                for ( size_t Ivar = 0; Ivar < Nvars; Ivar++) {
                    coarsened_grid->variables.at(vars_to_filter[Ivar])[Ipt] = solution_grid.variables.at(vars_to_filter[Ivar])[Ipt];
                }
                    continue;
            }

            // Initialize kA
            double dArea = solution_grid.areas[Ipt],
                   kernel_normalization = dArea;

            for (size_t Ivar = 0; Ivar < Nvars; Ivar++) {
                double local_val = solution_grid.variables.at(vars_to_filter[Ivar])[Ipt];
                value_sums.at(Ivar) = dArea * local_val;
                // not +=, here we re-set the values to the local value
                // before looping to accumulate over the kernel
            }

            // Next, seed the 'points to test' with the adjacent points
            std::deque<size_t> points_to_test;
            for (size_t II = 0; II < solution_grid.num_neighbours; II++ ) {
                size_t neighbour_index = solution_grid.adjacency_indices[Ipt][II];
                points_to_test.push_back( neighbour_index );
                planned_for_testing[ neighbour_index ] = true;
            }

            // So long as we still have points to test, keep testing!
            // This implements a breadth-first search through the adjacency matrix to build
            //      filtering kernel. If a point is within distance, add it to the kernel,
            //      and then test it's neighbours. If those are in, test their neighbours,
            //      and so on. If a point is too far away, we do not test it's neighbours.
            // This then assumes that for any point X within distance L of Y, that there
            //      is an adjacency path from Y to X strictly using points within distance
            //      L of Y.
            // Along the way, test for double inclusing to make sure that we don't test 
            //      points repeatedly etc.
            while ( points_to_test.size() > 0 ) {

                // Pull out the most-recently-added point, and remove it from the 'to test' list,
                // since we're testing it now.
                // Since we're pulling out the most-recently-added, that effectively makes this a
                // depth-first search to build the kernel.
                size_t Jpt = points_to_test.front();
                points_to_test.pop_front();
                planned_for_testing[ Jpt ] = false;

                double Jlat = solution_grid.latitude[Jpt];
                double Jlon = solution_grid.longitude[Jpt];
                double local_dist = distance( target_lon, target_lat, Jlon, Jlat );

                double kern_val = kernel( local_dist, filter_scale );

                if ( ( kern_val > 1e-10 ) or ( local_dist <= 0.5 * filter_scale )  ) {
                    was_accepted[Jpt] = true;

                    // Accumulate the coarse values
                    dArea = solution_grid.areas[Jpt];
                    kernel_normalization += kern_val * dArea;
                    for ( size_t Ivar = 0; Ivar < Nvars; Ivar++) {
                        double local_val = solution_grid.variables.at(vars_to_filter[Ivar])[Jpt];
                        value_sums.at(Ivar) += kern_val * dArea * local_val;
                    }

                    // Since this point is in the kernel
                    // add its neighbours to the 'to search' list
                    for ( size_t II = 0; II < solution_grid.num_neighbours; II++ ) {
                        size_t neighbour_index = solution_grid.adjacency_indices[Jpt][II];

                        // but first check if that neighbour has been rejected already
                        if ( was_rejected[neighbour_index] ) { continue; }

                        // then check if that neighbour is already accepted
                        if ( was_accepted[neighbour_index] ) { continue; }

                        // then check if that neighbour is already on the search list
                        if ( planned_for_testing[neighbour_index] ) { continue; }

                        // if it's not on any of those list already
                        // then add it to the 'points to test'
                        points_to_test.push_back( neighbour_index );
                        planned_for_testing[neighbour_index] = true;

                    }
                } else {
                    // Otherwise, record this point as rejected, and move on
                    was_rejected[Jpt] = true;
                }

            } // loop through to make the kernel

            // Store the filtered values
            for ( size_t Ivar = 0; Ivar < Nvars; Ivar++ ) {
                coarsened_grid->variables.at(vars_to_filter[Ivar])[Ipt] = 
                    ( kernel_normalization == 0 )
                    ?
                    0.
                    :
                    value_sums[Ivar] / kernel_normalization;
            }
        }
    }

    #pragma omp parallel default(none) private(Ipt) \
    shared( coarsened_grid, solution_grid ) \
    firstprivate( Npts_coarse, alpha )
    {
        #pragma omp for collapse(1) schedule(static)
        for ( Ipt = 0; Ipt < Npts_coarse; Ipt++ ) {

            // alpha = 1 -> full filtering
            // alpha = 0 -> no filtering
            //const double alpha = dt * iters_per_batch / ( dt * 1e3 );
            //const double alpha = 1.;

            double small_scale_Psi = solution_grid.variables.at("Psi")[Ipt]
                - coarsened_grid->variables.at("Psi")[Ipt];
            solution_grid.variables.at("Psi")[Ipt] -= alpha * small_scale_Psi;

            double small_scale_Phi = solution_grid.variables.at("Phi")[Ipt]
                - coarsened_grid->variables.at("Phi")[Ipt];
        solution_grid.variables.at("Phi")[Ipt] -= alpha * small_scale_Phi;
        }
    }

}





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
        const double hyper_visc_coeff,
        const MPI_Comm comm
        ) {

    if ( not(use_vort_div and use_vel) ) {
        throw std::runtime_error("Must use both vortdiv and vels for time-step solver.");
    }

    const bool Dufort_Frankel = false;

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
    for ( int refine_level = num_refinements; refine_level >= 0; refine_level-- ) {

        #if DEBUG >= 0
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


            // Downsample the velocities onto the coarse grid
            std::vector<std::string> vars_to_map = { "u_lon", "u_lat", "vort", "div", "Lap4_Psi", "Lap4_Phi" };

            initialize_coarsened_grid( *coarsened_grid, source_data, Nlat_coarse, Nlon_coarse );
            Npts_coarse = coarsened_grid->longitude.size(); // The coarsened grid might have a different size
            for ( size_t Ivar = 0; Ivar < vars_to_map.size(); Ivar++ ) {
                coarsened_grid->variables.insert( std::pair< std::string, std::vector<double> >(
                            vars_to_map[Ivar], std::vector<double>(Npts_coarse, 0.) ) );
            }

            #if DEBUG >= 2
            fprintf(stdout, "Mapping velocities onto the coarse grid.\n");
            #endif
            vars_to_map = { "u_lon", "u_lat" };
            if (constants::DO_TIMING) { clock_on = MPI_Wtime(); }
            map_grid_to_grid( source_data, *coarsened_grid, vars_to_map );
            if (constants::DO_TIMING) { timing_records.add_to_record(MPI_Wtime() - clock_on, "Downsampling"); }

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
            if (constants::DO_TIMING) { timing_records.add_to_record(MPI_Wtime() - clock_on, "Computing Vort and Div"); }
        } else {
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
            if (constants::DO_TIMING) { timing_records.add_to_record(MPI_Wtime() - clock_on, "Downsampling"); }
        }

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


        // Add Psi and Phi variables to our coarsened grid
        coarsened_grid->variables.insert( std::pair< std::string, std::vector<double> >(
                    "Psi", std::vector<double>(Npts_coarse, 0.) ) );
        coarsened_grid->variables.insert( std::pair< std::string, std::vector<double> >(
                    "Phi", std::vector<double>(Npts_coarse, 0.) ) );

        // Now, up-sample the previous solutions onto the new grid
        if ( refine_level < num_refinements ) {

            std::vector<std::string> vars_to_map = { "Psi", "Phi" };

            if (constants::DO_TIMING) { clock_on = MPI_Wtime(); }
            map_grid_to_grid( solution_grid, *coarsened_grid, vars_to_map );
            if (constants::DO_TIMING) { timing_records.add_to_record(MPI_Wtime() - clock_on, "Upsampling"); }
        }

        // Re-initialize the solution grid on the new computational grid
        // To avoid having to re-build, just copy from the coarsened grid
        solution_grid.clear();
        solution_grid.copy_from_ptr( coarsened_grid );

        solution_grid.variables.insert( std::pair< std::string, std::vector<double> >(
                    "Psi", std::vector<double>(Npts_coarse, 0.) ) );
        solution_grid.variables.insert( std::pair< std::string, std::vector<double> >(
                    "Phi", std::vector<double>(Npts_coarse, 0.) ) );

        // Reset the Helmholtz data for this resolution
        Helmholtz_data.clear();
        Helmholtz_data.collapse_land = collapse_land;
        Helmholtz_data.use_vort_div = use_vort_div;
        Helmholtz_data.use_vel = use_vel;
        Helmholtz_data.weight_err = weight_err;
        Helmholtz_data.tolerance = rel_tol;
        Helmholtz_data.iterations_per_cycle = iters_per_batch;
        Helmholtz_data.iteration_max = max_iters;

        Helmholtz_data.IdentifyLandlockedPoints( *coarsened_grid, use_vort_div );
        Helmholtz_data.CreateLandCollapsingMap( *coarsened_grid );

        IdentifyCoastalBoundaries( Helmholtz_data, coarsened_grid );

        // We're done with the land mask, so remove it. This causes vort and
        // div to be computed on land cells to ensure agreement
        std::fill( coarsened_grid->mask.begin(), coarsened_grid->mask.end(), true );
        std::fill( solution_grid.mask.begin(), solution_grid.mask.end(), true ); 


        if ( refine_level == num_refinements ) {
            // Use spherical harmonics to set the seed on the coarsest grid
            if (constants::DO_TIMING) { clock_on = MPI_Wtime(); }
            SphericalHarmonicSolver(  coarsened_grid->variables.at("Psi"),
                                     *coarsened_grid,
                                      coarsened_grid->variables.at("vort") );
            SphericalHarmonicSolver(  coarsened_grid->variables.at("Phi"),
                                     *coarsened_grid,
                                      coarsened_grid->variables.at("div") );
            if (constants::DO_TIMING) { timing_records.add_to_record(MPI_Wtime() - clock_on, 
                                         "SphericalHarmonicSolver"); }
        }

        //
        //// Output the seed / starting guess
        //
        std::vector<std::string> vars_to_add = { 
            "projected_u_lon", "projected_u_lat", "projected_vort", "projected_div",
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

        coarsened_grid->variables.at("projected_vort") = coarsened_grid->variables.at("proj_vort");
        coarsened_grid->variables.at("projected_div") = coarsened_grid->variables.at("proj_div");

        for ( size_t Ipt = 0; Ipt < Npts_coarse; Ipt++ ) {
            coarsened_grid->variables.at("projected_u_lon")[Ipt] = 
                  coarsened_grid->variables.at("proj_uo_tor")[Ipt] 
                + coarsened_grid->variables.at("proj_uo_pot")[Ipt];
            coarsened_grid->variables.at("projected_u_lat")[Ipt] = 
                  coarsened_grid->variables.at("proj_vo_tor")[Ipt] 
                + coarsened_grid->variables.at("proj_vo_pot")[Ipt];
        }

        std::ostringstream seed_stream;
        seed_stream << "seed_R" << (refine_level) << ".nc";
        const std::string seed_filename = seed_stream.str();
        solution_grid.variables.at("Psi") = coarsened_grid->variables.at("Psi");
        solution_grid.variables.at("Phi") = coarsened_grid->variables.at("Phi");
        write_Helmholtz_output( seed_filename, Helmholtz_data, *coarsened_grid, solution_grid );


        // For DuFort-Frankel, we need the self-weights of the Laplacian
        std:: vector<double> w_n_n;
        double max_w_n_n = 0, min_da = 4 * M_PI * pow(6371e3, 2);
        if ( Dufort_Frankel) {
            w_n_n.resize( Npts_coarse );
            #pragma omp parallel default(none) private(Ipt) \
            shared( Helmholtz_data, coarsened_grid, w_n_n ) firstprivate( Npts_coarse ) \
            reduction( max:max_w_n_n )
            {
                #pragma omp for collapse(1) schedule(static)
                for (Ipt = 0; Ipt < Npts_coarse; Ipt++ ) {

                    // Only time-step on water and coastal points
                    if ( Helmholtz_data.all_land_neighbours[Ipt] ) { continue; }

                    double cos_lat = cos( coarsened_grid->latitude[Ipt] );
                    double tan_lat = tan( coarsened_grid->latitude[Ipt] );
                    w_n_n[Ipt] = 
                        (
                              coarsened_grid->adjacency_d2dlon2_weights[Ipt].back() / pow(cos_lat, 2)
                            + coarsened_grid->adjacency_d2dlat2_weights[Ipt].back()
                            - coarsened_grid->adjacency_ddlat_weights[  Ipt].back() * tan_lat
                        ) / pow( constants::R_earth, 2);

                    max_w_n_n = std::fmax( max_w_n_n, std::fabs( w_n_n[Ipt] ) );
                }
            }
        } else {
            #pragma omp parallel default(none) private(Ipt) \
            shared( Helmholtz_data, coarsened_grid ) firstprivate( Npts_coarse ) \
            reduction( min:min_da )
            {
                #pragma omp for collapse(1) schedule(static)
                for (Ipt = 0; Ipt < Npts_coarse; Ipt++ ) {
                    // Only time-step on water and coastal points
                    if ( Helmholtz_data.all_land_neighbours[Ipt] ) { continue; }
                    min_da = std::fmin( min_da, coarsened_grid->areas[Ipt] );
                }
            }
            fprintf( stdout,  "min(dA) = %.2e km^2\n", min_da / 1e6 );
        }
        const double typical_spacing = Dufort_Frankel ? sqrt(2. / max_w_n_n) : sqrt( min_da );

        bool keep_solving = true;
        unsigned long int iter_cycle = 0;

        std::vector<double> dPsidt_prev1, dPsidt_prev2, dPhidt_prev1, dPhidt_prev2, Psi_prev, Phi_prev;
        if ( Dufort_Frankel ) {
            Psi_prev.resize( Npts_coarse, 0 );
            Phi_prev.resize( Npts_coarse, 0 );
        } else {
            dPsidt_prev1.resize( Npts_coarse, 0 );
            dPsidt_prev2.resize( Npts_coarse, 0 );
            dPhidt_prev1.resize( Npts_coarse, 0 );
            dPhidt_prev2.resize( Npts_coarse, 0 );
        }
        const double hyper_visc = hyper_visc_coeff * pow( typical_spacing, 2. );

        // If we're not on the finest grid, filter to smooth out the poles a bit
        const double filter_scale = 2 * typical_spacing;
        //const double alpha = std::min( dt * iters_per_batch / ( 20 * dt / CFL ), 1. );
        const double alpha = 1.;
        const bool only_filter_pole = true;
        if ( refine_level > 0 ) {
            std::vector<std::string> init_vars_to_filter = { "vort", "div" };
            filter_scalars( coarsened_grid, *coarsened_grid, init_vars_to_filter, 
                    filter_scale, alpha, only_filter_pole );
        }

        while ( keep_solving ) {

            const double dt = ( iter_cycle < 5 ? CFL/10 : CFL) * pow( typical_spacing, 2. );
            if ( iter_cycle == 5 ) {
                fprintf( stdout, "Using time-step of %.2es, hyperviscosity of %.2e, and 'spacing' of %.2ekm\n", 
                        dt, hyper_visc, typical_spacing/1e3 );
            }

            if (constants::DO_TIMING) { clock_on = MPI_Wtime(); }

            #pragma omp parallel default(none) private(Ipt) \
            shared(Helmholtz_data, solution_grid, coarsened_grid, w_n_n, \
                    Psi_prev, Phi_prev, dPsidt_prev1, dPsidt_prev2, dPhidt_prev1, dPhidt_prev2 ) \
            firstprivate( dt, Npts_coarse, hyper_visc )
            {
                #pragma omp for collapse(1) schedule(static)
                for ( Ipt = 0; Ipt < Npts_coarse; Ipt++ ) {

                    // Only solve on coastal / water cells. Land will be filled in separately
                    if ( Helmholtz_data.all_land_neighbours[Ipt] ) { continue; }


                    if ( Dufort_Frankel ) {
                        // Dufort-Frankel version
                        // d(Psi)/dt = Lap(Psi) - vorticity
                        double del_Psi = 0.5 * ( solution_grid.variables.at("Psi")[Ipt] - Psi_prev[Ipt] );
                        double dPsi_dt = coarsened_grid->variables.at("proj_vort")[Ipt] 
                                         - coarsened_grid->variables.at("vort")[Ipt]
                                         - w_n_n[Ipt] * del_Psi;

                        // d(Phi)/dt = Lap(Phi) - divergence
                        double del_Phi = 0.5 * ( solution_grid.variables.at("Phi")[Ipt] - Phi_prev[Ipt] );
                        double dPhi_dt = coarsened_grid->variables.at("proj_div")[Ipt] 
                                        - coarsened_grid->variables.at("div")[Ipt]
                                            - w_n_n[Ipt] * del_Phi;

                        // Update previous values
                        Psi_prev[Ipt] = solution_grid.variables.at("Psi")[Ipt];
                        Phi_prev[Ipt] = solution_grid.variables.at("Phi")[Ipt];

                        // Update current values
                        const double DF_coef = dt / ( 1 - 0.5 * dt * w_n_n[Ipt] ); // Dufort-Frankel 'implicit' facctor
                        solution_grid.variables.at("Psi")[Ipt] += DF_coef * dPsi_dt;
                        solution_grid.variables.at("Phi")[Ipt] += DF_coef * dPhi_dt;
                    } else { 

                        // Otherwise, AB3

                        // d(Psi)/dt = Lap(Psi) - vorticity
                        double dPsi_dt = coarsened_grid->variables.at("proj_vort")[Ipt] 
                                         - coarsened_grid->variables.at("vort")[Ipt]
                                         - hyper_visc * coarsened_grid->variables.at("Lap4_Psi")[Ipt]; // hyperviscosity

                        // d(Phi)/dt = Lap(Phi) - divergence
                        double dPhi_dt = coarsened_grid->variables.at("proj_div")[Ipt] 
                                        - coarsened_grid->variables.at("div")[Ipt]
                                        - hyper_visc * coarsened_grid->variables.at("Lap4_Phi")[Ipt]; // hyperviscosity

                        double dPsidt_p  = dPsidt_prev1[Ipt],
                               dPsidt_pp = dPsidt_prev2[Ipt],
                               dPhidt_p  = dPhidt_prev1[Ipt],
                               dPhidt_pp = dPhidt_prev2[Ipt];

                               // Update the stored value
                               solution_grid.variables.at("Psi")[Ipt] += ( (23./12) * dPsi_dt - (16./12) * dPsidt_p + (5./12) * dPsidt_pp ) * dt;
                               solution_grid.variables.at("Phi")[Ipt] += ( (23./12) * dPhi_dt - (16./12) * dPhidt_p + (5./12) * dPhidt_pp ) * dt;

                               dPsidt_prev2[Ipt] = dPsidt_prev1[Ipt];
                               dPhidt_prev2[Ipt] = dPhidt_prev1[Ipt];

                               dPsidt_prev1[Ipt] = dPsi_dt;
                               dPhidt_prev1[Ipt] = dPhi_dt;
                    }
                }
            }

            // Now fill in the land interiours with their coastal means
            ReconstructCoastalInterior( Helmholtz_data, solution_grid );

            // Compute velocity components and vort,div
            // To measure errors and use for next iteration
            toroidal_vel_from_F( coarsened_grid->variables.at("proj_uo_tor"),
                    coarsened_grid->variables.at("proj_vo_tor"),
                    solution_grid.variables.at("Psi"),
                    *coarsened_grid,
                    coarsened_grid->mask );

            potential_vel_from_F( coarsened_grid->variables.at("proj_uo_pot"),
                    coarsened_grid->variables.at("proj_vo_pot"),
                    solution_grid.variables.at("Phi"),
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

            /*
            scalar_laplacian( 
                    coarsened_grid->variables.at("proj_vort"),
                    solution_grid.variables.at("Psi"),
                    *coarsened_grid,
                    coarsened_grid->mask
                    );

            scalar_laplacian( 
                    coarsened_grid->variables.at("proj_div"),
                    solution_grid.variables.at("Phi"),
                    *coarsened_grid,
                    coarsened_grid->mask
                    );
            */

            // Hyperviscosity terms
            if ( hyper_visc != 0 ) {
                scalar_laplacian( 
                        coarsened_grid->variables.at("Lap4_Psi"),
                        coarsened_grid->variables.at("proj_vort"),
                        *coarsened_grid,
                        coarsened_grid->mask
                        );

                scalar_laplacian( 
                        coarsened_grid->variables.at("Lap4_Phi"),
                        coarsened_grid->variables.at("proj_div"),
                        *coarsened_grid,
                        coarsened_grid->mask
                        );
            }

            iter_cycle++;
            if ( iter_cycle > max_iters ) { keep_solving = false; }
            if (constants::DO_TIMING) { timing_records.add_to_record(MPI_Wtime() - clock_on, "Solving"); }

            if ( iter_cycle % iters_per_batch == 0 ) {

                // Shift so Psi[0] = 0, Phi[0] = 0
                #pragma omp parallel default(none) private(Ipt) \
                shared( solution_grid ) firstprivate( Npts_coarse )
                {
                    #pragma omp for collapse(1) schedule(static)
                    for ( Ipt = 1; Ipt < Npts_coarse; Ipt++ ) {
                        solution_grid.variables.at("Psi")[Ipt] -= solution_grid.variables.at("Psi")[0];
                        solution_grid.variables.at("Phi")[Ipt] -= solution_grid.variables.at("Phi")[0];
                    }
                }
                solution_grid.variables.at("Psi")[0] = 0;
                solution_grid.variables.at("Phi")[0] = 0;
                
                // Apply a low-grade filter
                std::vector<std::string> vars_to_filter = { "Psi", "Phi" };
                filter_scalars( coarsened_grid, solution_grid, vars_to_filter, filter_scale, alpha, only_filter_pole );

                // compute errors and check convergence
                double vel_2_err  = 0, vel_inf_err  = 0,
                       vort_2_err = 0, vort_inf_err = 0,
                       div_2_err  = 0, div_inf_err  = 0;
                double vel_2_norm  = 0, vel_inf_norm  = 0,
                       vort_2_norm = 0, vort_inf_norm = 0,
                       div_2_norm  = 0, div_inf_norm  = 0;
                double vort_2_ener = 0, div_2_ener = 0;

                get_norms(  vel_2_err,  vort_2_err,  div_2_err,
                            vel_2_norm, vort_2_norm, div_2_norm,
                            vel_inf_err,  vort_inf_err,  div_inf_err,
                            vel_inf_norm, vort_inf_norm, div_inf_norm,
                            vort_2_ener, div_2_ener,
                            coarsened_grid
                          );

                // If things get way off the rails, stop.
                if ( (vel_2_err > 1e3 * vel_2_norm) 
                        or (vort_2_err > 1e3 * vort_2_norm) 
                        or (div_2_err > 1e3 * div_2_norm) ) {
                    throw std::runtime_error("Errors have exceeded 500%. Aborting.");
                }

                while ( false and ( (vort_2_ener > 1.05 * vort_2_norm ) or (div_2_ener > 1.05 * div_2_norm ) ) ) {
                    // So long as we exceed energy by more than 5%, keep filtering
                    filter_scalars( coarsened_grid, solution_grid, vars_to_filter, filter_scale );

                    toroidal_vel_from_F( coarsened_grid->variables.at("proj_uo_tor"),
                            coarsened_grid->variables.at("proj_vo_tor"),
                            solution_grid.variables.at("Psi"),
                            *coarsened_grid,
                            coarsened_grid->mask );

                    potential_vel_from_F( coarsened_grid->variables.at("proj_uo_pot"),
                            coarsened_grid->variables.at("proj_vo_pot"),
                            solution_grid.variables.at("Phi"),
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

                    /*
                    scalar_laplacian( 
                            coarsened_grid->variables.at("proj_vort"),
                            solution_grid.variables.at("Psi"),
                            *coarsened_grid,
                            coarsened_grid->mask
                            );

                    scalar_laplacian( 
                            coarsened_grid->variables.at("proj_div"),
                            solution_grid.variables.at("Phi"),
                            *coarsened_grid,
                            coarsened_grid->mask
                            );
                    */

                    // Re-compute the norms
                    get_norms( vel_2_err,  vort_2_err,  div_2_err,
                            vel_2_norm, vort_2_norm, div_2_norm,
                            vel_inf_err,  vort_inf_err,  div_inf_err,
                            vel_inf_norm, vort_inf_norm, div_inf_norm,
                            vort_2_ener, div_2_ener,
                               coarsened_grid
                             );

                }

                // Save the convergence records
                Helmholtz_data.vel_2_errors.push_back(  vel_2_err  );
                Helmholtz_data.vort_2_errors.push_back( vort_2_err );
                Helmholtz_data.div_2_errors.push_back(  div_2_err  );

                Helmholtz_data.vel_inf_errors.push_back(  vel_inf_err  );
                Helmholtz_data.vort_inf_errors.push_back( vort_inf_err );
                Helmholtz_data.div_inf_errors.push_back(  div_inf_err  );

                // Save the convergence reference values
                Helmholtz_data.vel_2_norms.push_back(  vel_2_norm  );
                Helmholtz_data.vort_2_norms.push_back( vort_2_norm );
                Helmholtz_data.div_2_norms.push_back(  div_2_norm  );

                Helmholtz_data.vel_inf_norms.push_back(  vel_inf_norm  );
                Helmholtz_data.vort_inf_norms.push_back( vort_inf_norm );
                Helmholtz_data.div_inf_norms.push_back(  div_inf_norm  );

                // If sufficient convergence, stop
                if ( Helmholtz_data.IsConverged() ) { 
                    keep_solving = false;
                    fprintf( stdout,"  Halting at solver cycle %'zu. Solver has converged to desired tolerance.\n", iter_cycle );
                }

                #if DEBUG >= 0
                if (use_vort_div and use_vel) {
                    fprintf( stdout,"  Solver cycle %'zu complete. Relative 2-norm errors of vels, vort, and div are %.2e, %.2e, and %.2e\n", 
                            iter_cycle,
                            Helmholtz_data.vel_2_errors.back()  / Helmholtz_data.vel_2_norms.back(),
                            Helmholtz_data.vort_2_errors.back() / Helmholtz_data.vort_2_norms.back(), 
                            Helmholtz_data.div_2_errors.back()  / Helmholtz_data.div_2_norms.back() 
                           );
                } else if (use_vort_div) {
                    fprintf( stdout,"  Solver cycle %'zu complete. Relative 2-norm errors of vort and div are %.2e and %.2e\n", 
                            iter_cycle,
                            Helmholtz_data.vort_2_errors.back() / Helmholtz_data.vort_2_norms.back(), 
                            Helmholtz_data.div_2_errors.back()  / Helmholtz_data.div_2_norms.back() 
                           );
                } else {
                    fprintf( stdout,"  Solver cycle %'zu complete. Relative 2-norm errors of vels is %.2e\n", 
                            iter_cycle,
                            Helmholtz_data.vel_2_errors.back()  / Helmholtz_data.vel_2_norms.back()
                           );
                }
                #endif
            }
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

        coarsened_grid->variables.at("projected_vort") = coarsened_grid->variables.at("proj_vort");
        coarsened_grid->variables.at("projected_div") = coarsened_grid->variables.at("proj_div");

        #pragma omp parallel default(none) private(Ipt) \
        shared( coarsened_grid ) firstprivate( Npts_coarse )
        {
            #pragma omp for collapse(1) schedule(static)
            for ( size_t Ipt = 0; Ipt < Npts_coarse; Ipt++ ) {
                coarsened_grid->variables.at("projected_u_lon")[Ipt] = 
                      coarsened_grid->variables.at("proj_uo_tor")[Ipt] 
                    + coarsened_grid->variables.at("proj_uo_pot")[Ipt];
                coarsened_grid->variables.at("projected_u_lat")[Ipt] = 
                      coarsened_grid->variables.at("proj_vo_tor")[Ipt] 
                    + coarsened_grid->variables.at("proj_vo_pot")[Ipt];
            }
        }

        std::ostringstream name_stream;
        name_stream << "projection_R" << (refine_level) << ".nc";
        const std::string output_filename = name_stream.str();
        write_Helmholtz_output( output_filename, Helmholtz_data, *coarsened_grid, solution_grid );
        add_attr_to_file( "CFL", CFL, output_filename );
        add_attr_to_file( "hyper_viscosity", hyper_visc, output_filename );
        add_attr_to_file( "dt", CFL * pow( typical_spacing, 2. ), output_filename );
        if (constants::DO_TIMING) { timing_records.add_to_record(MPI_Wtime() - clock_on, "ProcessingSolution"); }

        // If we're doing timings, then print out and reset values now
        if (constants::DO_TIMING) { 
            timing_records.print();
            timing_records.reset();
            fflush(stdout);
        }

    }
}
