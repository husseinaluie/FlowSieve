#include "../constants.hpp"
#include "../functions.hpp"
#include "../netcdf_io.hpp"
#include "../preprocess.hpp"
#include "../differentiation_tools.hpp"
#include <algorithm>
#include <vector>
#include <deque>
#include <omp.h>
#include <math.h>
//#include <eigen3/Eigen/Sparse>
//#include <eigen3/Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>


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
        const MPI_Comm comm
        ) {

    int wRank, wSize;
    MPI_Comm_rank( comm, &wRank );
    MPI_Comm_size( comm, &wSize );

    srand( time(NULL) );

    // Create some tidy names for variables
    const std::vector<double>   &latitude   = source_data.latitude,
                                &longitude  = source_data.longitude,
                                &dAreas     = source_data.areas;

    const std::vector<bool> &mask = (constants::FILTER_OVER_LAND) ? source_data.reference_mask : source_data.mask;

    const std::vector<int>  &myCounts = source_data.myCounts,
                            &myStarts = source_data.myStarts;

    std::vector<double>   &u_lat = source_data.variables.at("u_lat"),
                          &u_lon = source_data.variables.at("u_lon");

    // Create a 'no mask' mask variable
    //   we'll treat land values as zero velocity
    //   We do this because including land seems
    //   to introduce strong numerical issues
    const std::vector<bool> unmask(mask.size(), true);

    const int   Ntime   = myCounts.at(0),
                Ndepth  = myCounts.at(1);

    const size_t Npts = latitude.size();
    const size_t num_neighbours = source_data.num_neighbours;

    int Itime=0, Idepth=0;
    size_t Ipt, index, neighbour_ind, index_sub, iters_used = 0;

    // Fill in the land areas with zero velocity
    #pragma omp parallel default(none) shared( u_lon, u_lat, mask, stderr, wRank ) private( index )
    {
        #pragma omp for collapse(1) schedule(guided)
        for (index = 0; index < u_lon.size(); index++) {
            if (not(mask.at(index))) {
                u_lon.at(index) = 0.;
                u_lat.at(index) = 0.;
            } else if (    ( std::fabs( u_lon.at(index) ) > 30000.) 
                        or ( std::fabs( u_lat.at(index) ) > 30000.) 
                      ) {
                fprintf( stderr, "  Rank %d found a bad vel point at index %'zu! Setting to zero.\n", wRank, index );
                u_lon.at(index) = 0.;
                u_lat.at(index) = 0.;
            }
        }
    }

    // Identify which points have land-only neighbours
    std::vector<short int> all_land_neighbours(Npts, 0); // 1 = yes, 0 = no
    size_t num_coastal = 0;
    for (Ipt = 0; Ipt < Npts; Ipt++) {
        all_land_neighbours[Ipt] = 1;
        if (mask[Ipt]) { all_land_neighbours[Ipt] = 0; }
        for (neighbour_ind = 0; neighbour_ind < num_neighbours; neighbour_ind++ ) {
            if ( mask[source_data.adjacency_indices.at(Ipt)[neighbour_ind] ] ) {
                all_land_neighbours[Ipt] = 0;
            }
        }
        if ( ( all_land_neighbours[Ipt] == 0 ) and ( not(mask[Ipt]) ) ) {
            num_coastal++;
        }
    }


    // We're going to eliminate land from the solver by hard-enforcing that all points
    // over land have constants Phi,Psi
    // Here, we compute the contiguous land blocks
    std::vector<size_t> pt_maps_to( Npts, Npts ); // size Npts, starting value Npts
    std::deque<size_t> points_to_test;
    size_t num_mapped_points = 0, num_mapped_onto = 0, num_mapped_coastal = 0, Ineighbour, Itest;
    for (Ipt = 0; Ipt < Npts; Ipt++) {
        //if (all_land_neighbours[Ipt] == 0) { 
            // Water+coast points map to themselves
        if (mask[Ipt]) { 
            // Water points map to themselves
            pt_maps_to[Ipt] = Ipt; 
        } else {
            // Check if we've already mapped this point
            if (pt_maps_to[Ipt] == Npts) {
                // Begin a depth-first search through the adjacency
                // matrix, mapping all land-only neighbours to this one
                // This point will become the representative for this 'island'
                points_to_test.clear();
                pt_maps_to[Ipt] = Ipt; 
                num_mapped_onto++;
                for (neighbour_ind = 0; neighbour_ind < num_neighbours; neighbour_ind++ ) {
                    points_to_test.push_back( source_data.adjacency_indices.at(Ipt)[neighbour_ind] );
                }
                // So long as we still have points to test, keep testing!
                while ( points_to_test.size() > 0 ) {
                    Itest = points_to_test.front();
                    points_to_test.pop_front();
                    if ( pt_maps_to[Itest] < Npts ) { continue; } // already mapped, skip
                    //else if ( all_land_neighbours[Itest] == 0) { continue; } // water/coast, so skip
                    else if (mask[Itest]) { continue; } // water, so skip
                    else {
                        // Otherwise, map it to Ipt, and add its neighbours to the test list
                        pt_maps_to[Itest] = Ipt;
                        num_mapped_points++;
                        for (neighbour_ind = 0; neighbour_ind < num_neighbours; neighbour_ind++ ) {
                            Ineighbour = source_data.adjacency_indices.at(Itest)[neighbour_ind];
                            if ( pt_maps_to[Ineighbour] == Npts ) {
                                // Only add the neighbour if we haven't already looked at it
                                points_to_test.push_back( Ineighbour );
                            }
                        }
                    }
                }
            }
        }
    }
    fprintf( stdout, "Mapping %'zu land points onto %'zu 'islands' (contiguous land masses). %'zu points are coastal.\n",
          num_mapped_points+num_mapped_onto, num_mapped_onto, num_coastal );
    const size_t Npts_mapped = Npts - num_mapped_points,
          Ncol = Npts - num_mapped_points,
          Nrow = Npts - (num_mapped_points+num_mapped_onto) + num_coastal;

    std::vector<size_t> num_mapped_before( Npts, 0 );
    std::vector<size_t> num_mapped_before_noncoastal( Npts, 0 );
    size_t counter = 0, noncoastal_counter = 0;
    for (Ipt = 1; Ipt < Npts; Ipt++) {
        if ( pt_maps_to[Ipt-1] != (Ipt-1) ) { counter++; }
        num_mapped_before[Ipt] = counter;

        //if ( ( pt_maps_to[Ipt-1] != (Ipt-1) ) and ( all_land_neighbours[Ipt-1] == 1) ) { 
        if ( all_land_neighbours[Ipt-1] == 1 ) { 
            noncoastal_counter++; 
        }
        num_mapped_before_noncoastal[Ipt] = noncoastal_counter;
    }

    // Storage vectors
    std::vector<double> 
        full_Psi(        u_lon.size(), 0. ),
        full_Phi(        u_lon.size(), 0. ),
        full_u_lon_tor(  u_lon.size(), 0. ),
        full_u_lat_tor(  u_lon.size(), 0. ),
        full_u_lon_pot(  u_lon.size(), 0. ),
        full_u_lat_pot(  u_lon.size(), 0. ),
        u_lon_tor_seed(  Npts, 0. ),
        u_lat_tor_seed(  Npts, 0. ),
        u_lon_pot_seed(  Npts, 0. ),
        u_lat_pot_seed(  Npts, 0. );

    // alglib variables
    std::vector<double> 
        RHS_vector( 4*Nrow+2 , 0. ),
        Psi_seed(     Npts, 0. ),
        Phi_seed(     Npts, 0. ),
        work_arr(     Npts, 0. ),
        div_term(     Npts, 0. ),
        vort_term(    Npts, 0. ),
        u_lon_rem(    Npts, 0. ),
        u_lat_rem(    Npts, 0. );
    

    // Copy the starting seed.
    if (single_seed) {
        #pragma omp parallel \
        default(none) \
        shared(Psi_seed, Phi_seed, seed_tor, seed_pot) \
        private( index ) \
        firstprivate( Npts )
        {
            #pragma omp for collapse(1) schedule(static)
            for (index = 0; index < Npts; ++index) {
                Psi_seed.at(index) = seed_tor.at(index);
                Phi_seed.at(index) = seed_pot.at(index);
            }
        }
    }

    // Get a magnitude for the derivatives, to help normalize the rows of the 
    //  Laplace entries to have similar magnitude to the others.
    long double deriv_ref_1 = 0;
    #pragma omp parallel default(none) \
    private( Ineighbour, index ) \
    shared( source_data ) \
    firstprivate( num_neighbours, Npts ) \
    reduction( +:deriv_ref_1 )
    {
        #pragma omp for collapse(1) schedule(static)
        for (index = 0; index < Npts; ++index) {
            for ( Ineighbour = 0; Ineighbour < num_neighbours + 1; Ineighbour++ ) {
                deriv_ref_1 += std::fabs( source_data.adjacency_ddlat_weights.at(index).at(Ineighbour) ) / (Npts*num_neighbours);
            }
        }
    }
    const double deriv_scale_factor = deriv_ref_1;
    fprintf( stdout, "deriv-scale-factor: %g\n", deriv_scale_factor );

    //
    //// Build the LHS part of the problem
    //      Ordering is: [  u_from_psi      u_from_phi   ] *  [ psi ]   =    [  u   ]
    //                   [  v_from_psi      v_from_phi   ]    [ phi ]        [  v   ]
    //
    //      Ordering is: [           - ddlat   sec(phi) * ddlon   ] *  [ psi ]   =    [     u     ]
    //                   [  sec(phi) * ddlon              ddlat   ]    [ phi ]        [     v     ]
    //      followed by entries forcing psi and phi to be constant over land
    

    #if DEBUG >= 0
    if (wRank == 0) {
        fprintf(stdout, "Building the LHS of the least squares problem.\n");
        fflush(stdout);
    }
    #endif

    double val;
    size_t column_skip, row_skip, land_counter = 0;

    const bool USE_TRUE_2ND_DERIV = false;
    double *F_array;
    typedef Eigen::Triplet<double> T;
    const int pts_per_1st_deriv = (constants::ADJACENCY_SIZE + 1),
              pts_per_2nd_deriv = pow(constants::ADJACENCY_SIZE + 1, USE_TRUE_2ND_DERIV ? 1 : 2);
    std::vector<T> Aij_triplets( 
            2 * Nrow * ( 2 * pts_per_1st_deriv + 2 * pts_per_2nd_deriv + pts_per_1st_deriv) + 2,
            T(0,0,0) );
    Aij_triplets[Aij_triplets.size() - 2] = T( 4 * Nrow + 0, 0,    Tikhov_Laplace ); // Force Psi = 0 on first point
    Aij_triplets[Aij_triplets.size() - 1] = T( 4 * Nrow + 1, Ncol, Tikhov_Laplace ); // Force Phi = 0 on first point
    size_t Itriplet, Ipt_mapped, neighbour_mapped;
    double weight_val, cos_lat_inv, R_inv, rand_val;
    bool is_pole;
    #pragma omp parallel default(none) \
    shared( mask, dAreas, latitude, source_data, Aij_triplets, \
            pt_maps_to, num_mapped_before, num_mapped_before_noncoastal, all_land_neighbours ) \
    private( Ipt, Ineighbour, neighbour_ind, Itriplet, row_skip, column_skip, is_pole, val, \
             weight_val, cos_lat_inv, R_inv, rand_val, counter, \
             Ipt_mapped, neighbour_mapped ) \
    firstprivate( Npts, Nrow, Ncol, weight_err, num_neighbours, Tikhov_Laplace, \
                  stdout, pts_per_1st_deriv, Npts_mapped, deriv_scale_factor )
    { 
        #pragma omp for collapse(1) schedule(static)
        for ( Ipt = 0; Ipt < Npts; Ipt++ ) {

            //if ( ( pt_maps_to[Ipt] != Ipt ) and (all_land_neighbours[Ipt] == 1) ) { continue; } // Skip points that were mapped
            if ( all_land_neighbours[Ipt] == 1 ) { continue; } // Skip points that were mapped

            weight_val = weight_err ? dAreas.at(Ipt) : 1.;
            cos_lat_inv = 1. / cos(latitude.at(Ipt));
            R_inv = 1. / constants::R_earth;

            for ( Ineighbour = 0; Ineighbour < num_neighbours + 1; Ineighbour++ ) {

                neighbour_ind = (Ineighbour < num_neighbours) ? 
                                        source_data.adjacency_indices.at(Ipt).at(Ineighbour) :
                                        Ipt;
                neighbour_ind = pt_maps_to[neighbour_ind]; // convert to mapped coordinated
                neighbour_mapped = neighbour_ind - num_mapped_before[neighbour_ind]; // all land removed from columns
                //Ipt_mapped = pt_maps_to[Ipt] - num_mapped_before_noncoastal[pt_maps_to[Ipt]]; // coast included in rows
                Ipt_mapped = Ipt - num_mapped_before_noncoastal[Ipt]; // coast included in rows
                if ( ( Ipt_mapped < 0 ) or (Ipt_mapped > Nrow) ) { 
                    fprintf( stdout, "BAD POINT! %'zu -> %'zu\n", Ipt, Ipt_mapped );
                    assert(false);
                }

                is_pole = std::fabs( std::fabs( latitude.at(Ipt) * 180.0 / M_PI ) - 90 ) < 1e-6;
                if ( is_pole ) { fprintf(stdout, "SKIPPING POLE POINT!\n"); continue; }

                // LON first derivative

                val  = source_data.adjacency_ddlon_weights.at(Ipt).at(Ineighbour);
                val *= weight_val * cos_lat_inv * R_inv;

                // Psi part (of u_lat)
                column_skip = 0 * Ncol + neighbour_mapped;
                row_skip    = 1 * Nrow + Ipt_mapped;
                Itriplet = Ipt_mapped * pts_per_1st_deriv + Ineighbour;
                Aij_triplets[Itriplet] = T( row_skip, column_skip, val );

                // Phi part (of u_lon)
                column_skip = 1 * Ncol + neighbour_mapped;
                row_skip    = 0 * Nrow + Ipt_mapped;
                Itriplet = (2*Nrow + Ipt_mapped) * pts_per_1st_deriv + Ineighbour;
                Aij_triplets[Itriplet] = T( row_skip, column_skip, val );


                // LAT first derivative

                val  = source_data.adjacency_ddlat_weights.at(Ipt).at(Ineighbour);
                val *= weight_val * R_inv;

                // Psi part (of u_lon)
                column_skip = 0 * Ncol + neighbour_mapped;
                row_skip    = 0 * Nrow + Ipt_mapped;
                Itriplet = (Nrow + Ipt_mapped) * pts_per_1st_deriv + Ineighbour;
                Aij_triplets[Itriplet] = T( row_skip, column_skip, -val );

                // Phi part (of u_lat)
                column_skip = 1 * Ncol + neighbour_mapped;
                row_skip    = 1 * Nrow + Ipt_mapped;
                Itriplet = (3*Nrow + Ipt_mapped) * pts_per_1st_deriv + Ineighbour;
                Aij_triplets[Itriplet] = T( row_skip, column_skip, val );

                // and do Lap for both Phi and Psi
                for ( counter = 0; counter < 2; counter++ ) {

                    // Second LON derivative
                    if ( USE_TRUE_2ND_DERIV ) {
                        val  = source_data.adjacency_d2dlon2_weights.at(Ipt).at(Ineighbour);
                        val *= weight_val * pow(cos_lat_inv * R_inv, 2.);
                        val *= Tikhov_Laplace / deriv_scale_factor;

                        column_skip = neighbour_mapped;
                        column_skip += counter * Ncol;
                        row_skip    = Ipt_mapped;
                        row_skip    += 2 * Nrow;
                        row_skip    += counter * Nrow;
                        Itriplet = Ipt_mapped * pts_per_2nd_deriv + Ineighbour;
                        Itriplet += 4 * Nrow * pts_per_1st_deriv;
                        Aij_triplets[Itriplet] = T( row_skip, column_skip, val );
                    } else {
                        for ( size_t D2_ind = 0; D2_ind < num_neighbours+1; D2_ind++ ) {
                            val  =   source_data.adjacency_ddlon_weights.at(Ipt).at(Ineighbour)
                                * source_data.adjacency_ddlon_weights.at(neighbour_ind).at(D2_ind);
                            val *= weight_val * pow(R_inv, 2.) * cos_lat_inv / cos(latitude.at(neighbour_ind));
                            val *= Tikhov_Laplace / deriv_scale_factor;

                            column_skip = source_data.adjacency_indices.at(neighbour_ind).at(D2_ind);
                            column_skip = pt_maps_to[column_skip];
                            column_skip = column_skip - num_mapped_before[column_skip];
                            column_skip += counter * Ncol;
                            row_skip    = Ipt_mapped;
                            row_skip    += 2 * Nrow;
                            row_skip    += counter * Nrow;
                            Itriplet = Ipt_mapped * pts_per_2nd_deriv + Ineighbour*(num_neighbours+1)  + D2_ind;
                            Itriplet += 4 * Nrow * pts_per_1st_deriv;
                            Aij_triplets[Itriplet] = T( row_skip, column_skip, val );
                        }
                    }

                    // Second LAT derivative
                    if ( USE_TRUE_2ND_DERIV ) {
                        val = source_data.adjacency_d2dlat2_weights.at(Ipt).at(Ineighbour);
                        val *= weight_val * pow(R_inv, 2.);
                        val *= Tikhov_Laplace / deriv_scale_factor;

                        column_skip = neighbour_mapped;
                        column_skip += counter * Ncol;
                        row_skip    = Ipt_mapped;
                        row_skip    += 2 * Nrow;
                        row_skip    += counter * Nrow;
                        Itriplet    = (Nrow + Ipt_mapped) * pts_per_2nd_deriv + Ineighbour;
                        Itriplet += 4 * Nrow * pts_per_1st_deriv;
                        Aij_triplets[Itriplet] = T( row_skip, column_skip, val );
                    } else {
                        for ( size_t D2_ind = 0; D2_ind < num_neighbours+1; D2_ind++ ) {
                            val  =   source_data.adjacency_ddlat_weights.at(Ipt).at(Ineighbour)
                                * source_data.adjacency_ddlat_weights.at(neighbour_ind).at(D2_ind);
                            val *= weight_val * pow(R_inv, 2.);
                            val *= Tikhov_Laplace / deriv_scale_factor;

                            column_skip = source_data.adjacency_indices.at(neighbour_ind).at(D2_ind);
                            column_skip = pt_maps_to[column_skip];
                            column_skip = column_skip - num_mapped_before[column_skip];
                            column_skip += counter * Ncol;
                            row_skip    = Ipt_mapped;
                            row_skip    += 2 * Nrow;
                            row_skip    += counter * Nrow;
                            Itriplet = (Nrow + Ipt_mapped) * pts_per_2nd_deriv + Ineighbour*(num_neighbours+1)  + D2_ind;
                            Itriplet += 4 * Nrow * pts_per_1st_deriv;
                            Aij_triplets[Itriplet] = T( row_skip, column_skip, val );
                        }
                    }

                    // First LAT derivative

                    val = - source_data.adjacency_ddlat_weights.at(Ipt).at(Ineighbour) * tan( latitude.at(Ipt) );
                    val *= weight_val * pow(R_inv, 2.);
                    val *= Tikhov_Laplace / deriv_scale_factor;

                    column_skip = neighbour_mapped;
                    column_skip += counter * Ncol;
                    row_skip    = Ipt_mapped;
                    row_skip    += 2 * Nrow;
                    row_skip    += counter * Nrow;
                    Itriplet    = 2 * Nrow * pts_per_2nd_deriv + Ipt_mapped * pts_per_1st_deriv + Ineighbour;
                    Itriplet += 4 * Nrow * pts_per_1st_deriv;
                    Aij_triplets[Itriplet] = T( row_skip, column_skip, val );
                }
            }
        }
    }

    Eigen::SparseMatrix<double> LHS_matr( 4*Nrow+2, 2*Ncol );
    fprintf( stdout, "%'zu, %'zu\n", 4*Nrow+2, 2*Ncol );
    LHS_matr.setFromTriplets( Aij_triplets.begin(), Aij_triplets.end() );

    #if DEBUG >= 1
    if (wRank == 0) {
        fprintf(stdout, "Declaring the least squares problem and computing.\n");
        fflush(stdout);
    }
    #endif

    LHS_matr.makeCompressed();
    Eigen::LeastSquaresConjugateGradient< Eigen::SparseMatrix<double> > solver;
    solver.setMaxIterations(max_iters);
    solver.setTolerance(rel_tol);
    //Eigen::SparseLU< Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > solver;
    //Eigen::SparseQR< Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > solver;
    solver.compute( LHS_matr );
    if ( solver.info() == Eigen::NumericalIssue ) {
        fprintf( stderr, "The provided data did not satisfy the prerequisites..\n" );
        return;
    } else if ( solver.info() == Eigen::NoConvergence ) {
        fprintf( stderr, "Iterative procedure did not converge.\n" );
        return;
    } else if ( solver.info() == Eigen::InvalidInput ) {
        fprintf( stderr, "The inputs are invalid, or the algorithm has been improperly called.\n" );
        return;
    } else if ( solver.info() != Eigen::Success ) {
        fprintf( stderr, "Eigen decomposition failed in an unknown way.\n" );
        return;
    }

    // Counters to track termination types
    int terminate_count_abs_tol = 0,
        terminate_count_rel_tol = 0,
        terminate_count_max_iter = 0,
        terminate_count_rounding = 0,
        terminate_count_other = 0;

    // Now do the solve!
    for (int Itime = 0; Itime < Ntime; ++Itime) {
        for (int Idepth = 0; Idepth < Ndepth; ++Idepth) {

            if (not(single_seed)) {
                #if DEBUG >= 2
                fprintf( stdout, "Extracting seed.\n" );
                fflush(stdout);
                #endif
                // If single_seed == false, then we were provided seed values, pull out the appropriate values here
                #pragma omp parallel \
                default(none) \
                shared( Psi_seed, Phi_seed, seed_tor, seed_pot, Itime, Idepth, stdout ) \
                private( index, index_sub ) \
                firstprivate( Ntime, Ndepth, Npts )
                {
                    #pragma omp for collapse(1) schedule(static)
                    for (index = 0; index < Npts; ++index) {
                        Psi_seed.at(index) = seed_tor.at(index + Npts*(Itime*Ndepth + Idepth));
                        Phi_seed.at(index) = seed_pot.at(index + Npts*(Itime*Ndepth + Idepth));
                    }
                }
            }

            // Get velocity from seed
            #if DEBUG >= 3
            fprintf( stdout, "Getting velocities from seed.\n" );
            fflush(stdout);
            #endif
            toroidal_vel_from_F(  u_lon_tor_seed, u_lat_tor_seed, Psi_seed, source_data, use_mask ? mask : unmask);
            potential_vel_from_F( u_lon_pot_seed, u_lat_pot_seed, Phi_seed, source_data, use_mask ? mask : unmask);

            #if DEBUG >= 3
            fprintf( stdout, "Subtracting seed velocity to get remaining.\n" );
            fflush(stdout);
            #endif
            #pragma omp parallel default(none) \
            shared( Itime, Idepth, stdout, \
                    u_lon, u_lon_tor_seed, u_lon_pot_seed, u_lon_rem, \
                    u_lat, u_lat_tor_seed, u_lat_pot_seed, u_lat_rem ) \
            private( index, index_sub ) \
            firstprivate( Ntime, Ndepth, Npts )
            {
                #pragma omp for collapse(1) schedule(static)
                for (index_sub = 0; index_sub < Npts; ++index_sub) {
                    index = index_sub + Npts*(Itime*Ndepth + Idepth);
                    u_lon_rem.at( index_sub ) = u_lon.at(index) - u_lon_tor_seed.at(index_sub) - u_lon_pot_seed.at(index_sub);
                    u_lat_rem.at( index_sub ) = u_lat.at(index) - u_lat_tor_seed.at(index_sub) - u_lat_pot_seed.at(index_sub);
                }
            }

            #if DEBUG >= 3
            fprintf( stdout, "Getting divergence and vorticity from remaining velocity.\n" );
            fflush(stdout);
            #endif
            toroidal_vel_div(        div_term, u_lon_rem, u_lat_rem, source_data, use_mask ? mask : unmask );
            toroidal_curl_u_dot_er( vort_term, u_lon_rem, u_lat_rem, source_data, use_mask ? mask : unmask );

            //
            //// Set up the RHS_vector
            //

            #if DEBUG >= 2
            if ( wRank == 0 ) {
                fprintf(stdout, "Building the RHS of the least squares problem.\n");
                fflush(stdout);
            }
            #endif
            
            double is_pole;
            std::fill( RHS_vector.begin(), RHS_vector.end(), 0. );
            #pragma omp parallel default(none) \
            shared( dAreas, RHS_vector, u_lon_rem, u_lat_rem, vort_term, div_term, \
                    num_mapped_before_noncoastal, pt_maps_to ) \
            private( Ipt ) \
            firstprivate( weight_err, Npts, Nrow, Tikhov_Laplace, deriv_scale_factor )
            {
                #pragma omp for collapse(1) schedule(static)
                for ( Ipt = 0; Ipt < Npts; ++Ipt) {
                    if ( pt_maps_to[Ipt] != Ipt ) { continue; }
                    RHS_vector.at(         Ipt - num_mapped_before_noncoastal[Ipt]) = 
                        u_lon_rem.at(Ipt) * ( weight_err ? dAreas.at(Ipt) : 1. );
                    RHS_vector.at(  Nrow + Ipt - num_mapped_before_noncoastal[Ipt]) = 
                        u_lat_rem.at(Ipt) * ( weight_err ? dAreas.at(Ipt) : 1. );
                    RHS_vector.at(2*Nrow + Ipt - num_mapped_before_noncoastal[Ipt]) = 
                        (Tikhov_Laplace / deriv_scale_factor) * vort_term.at(Ipt) * ( weight_err ? dAreas.at(Ipt) : 1. );
                    RHS_vector.at(3*Nrow + Ipt - num_mapped_before_noncoastal[Ipt]) = 
                        (Tikhov_Laplace / deriv_scale_factor) * div_term.at(Ipt) * ( weight_err ? dAreas.at(Ipt) : 1. );
                }
            }
            Eigen::VectorXd RHS = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(RHS_vector.data(), RHS_vector.size());

            //
            //// Now apply the least-squares solver
            //
            #if DEBUG >= 0
            if ( wRank == 0 ) {
                fprintf(stdout, "Solving the least squares problem.\n");
                fflush(stdout);
            }
            #endif
            Eigen::VectorXd F_Eigen = solver.solve( RHS );
            #if DEBUG >= 0
            if ( wRank == 0 ) {
                fprintf( stdout, "    Solver converged after %ld iterations to error %g.\n", 
                        solver.iterations(), solver.error() );
                fflush(stdout);
            }
            #endif
            std::vector<double> Psi_vector(Npts, 0), Phi_vector(Npts, 0);
            for (size_t ii = 0; ii < Npts; ++ii) {
                Psi_vector[ii] = F_Eigen[              pt_maps_to[ii] - num_mapped_before[pt_maps_to[ii]]];
                Phi_vector[ii] = F_Eigen[Npts_mapped + pt_maps_to[ii] - num_mapped_before[pt_maps_to[ii]]];
            }
            //std::vector<double> Psi_vector(F_Eigen.data(),        F_Eigen.data() +   Npts);
            //std::vector<double> Phi_vector(F_Eigen.data() + Npts, F_Eigen.data() + 2*Npts);

            // Add the seed back in
            for (size_t ii = 0; ii < Npts; ++ii) {
                Psi_vector.at(ii) += Psi_seed.at(ii);
                Phi_vector.at(ii) += Phi_seed.at(ii);
            }

            // Get velocity associated to computed F field
            #if DEBUG >= 2
            if ( wRank == 0 ) {
                fprintf(stdout, " Extracting velocities and divergence from toroidal field.\n");
                fflush(stdout);
            }
            #endif

            std::vector<double> u_lon_tor(Npts, 0.), u_lat_tor(Npts, 0.), u_lon_pot(Npts, 0.), u_lat_pot(Npts, 0.);
            toroidal_vel_from_F(  u_lon_tor, u_lat_tor, Psi_vector, source_data, use_mask ? mask : unmask);
            potential_vel_from_F( u_lon_pot, u_lat_pot, Phi_vector, source_data, use_mask ? mask : unmask);

            //
            //// Store into the full arrays
            //
            #if DEBUG >= 2
            if ( wRank == 0 ) {
                fprintf(stdout, " Storing values into output arrays\n");
                fflush(stdout);
            }
            #endif
            #pragma omp parallel \
            default(none) \
            shared( full_u_lon_tor, u_lon_tor, full_u_lat_tor, u_lat_tor, \
                    full_u_lon_pot, u_lon_pot, full_u_lat_pot, u_lat_pot, \
                    full_Psi, full_Phi, Psi_vector, Phi_vector, \
                    Phi_seed, Psi_seed, \
                    Itime, Idepth ) \
            private( index, index_sub ) \
            firstprivate( Ndepth, Ntime, single_seed, Npts )
            {
                #pragma omp for collapse(1) schedule(static)
                for (index_sub = 0; index_sub < Npts; ++index_sub) {
                    index = index_sub + Npts*(Itime*Ndepth + Idepth);

                    full_u_lon_tor.at(index) = u_lon_tor.at(index_sub) ;
                    full_u_lat_tor.at(index) = u_lat_tor.at(index_sub) ;

                    full_u_lon_pot.at(index) = u_lon_pot.at(index_sub) ;
                    full_u_lat_pot.at(index) = u_lat_pot.at(index_sub) ;

                    full_Psi.at(index) = Psi_vector.at( index_sub );
                    full_Phi.at(index) = Phi_vector.at( index_sub );

                    // If we don't have a seed for the next iteration, use this solution as the seed
                    if (single_seed) {
                        Psi_seed.at(index_sub) = Psi_vector.at(index_sub);
                        Phi_seed.at(index_sub) = Phi_vector.at(index_sub);
                    }
                }
            }

            #if DEBUG >= 0
            if ( source_data.full_Ndepth > 1 ) {
                fprintf(stdout, "  --  --  Rank %d done depth %d after %'zu iterations\n", wRank, Idepth + myStarts.at(1), iters_used );
                fflush(stdout);
            }
            #endif

        }

        #if DEBUG >= 0
        if ( source_data.full_Ntime > 1 ) {
            fprintf(stdout, " -- Rank %d done time %d after %'zu iterations\n", wRank, Itime + myStarts.at(0), iters_used );
            fflush(stdout);
        }
        #endif

        #if DEBUG >= 0
        if ( ( source_data.full_Ntime = 1 ) and ( source_data.full_Ndepth = 1 ) ) {
            fprintf(stdout, " -- Rank done after %'zu iterations\n", iters_used );
            fflush(stdout);
        }
        #endif
    }

    //
    //// Print termination counts
    //

    int total_count_abs_tol, total_count_rel_tol, total_count_max_iter, total_count_rounding, total_count_other;

    MPI_Reduce( &terminate_count_abs_tol,  &total_count_abs_tol,  1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD );
    MPI_Reduce( &terminate_count_rel_tol,  &total_count_rel_tol,  1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD );
    MPI_Reduce( &terminate_count_max_iter, &total_count_max_iter, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD );
    MPI_Reduce( &terminate_count_rounding, &total_count_rounding, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD );
    MPI_Reduce( &terminate_count_other,    &total_count_other,    1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD );

    #if DEBUG >= 0
    if (wRank == 0) {
        fprintf( stdout, "\n" );
        fprintf( stdout, "Termination counts: %'d from absolute tolerance\n", total_count_abs_tol );
        fprintf( stdout, "                    %'d from relative tolerance\n", total_count_rel_tol );
        fprintf( stdout, "                    %'d from iteration maximum\n", total_count_max_iter );
        fprintf( stdout, "                    %'d from rounding errors \n", total_count_rounding );
        fprintf( stdout, "                    %'d from other causes \n", total_count_other );
        fprintf( stdout, "\n" );
    }
    #endif


    //
    //// Write the output
    //

    const int ndims = 3;
    size_t starts[ndims] = {
        size_t(myStarts.at(0)), size_t(myStarts.at(1)), 0
    };
    size_t counts[ndims] = { size_t(Ntime), size_t(Ndepth), Npts };

    std::vector<std::string> vars_to_write;
    if (not(constants::MINIMAL_OUTPUT)) {
        vars_to_write.push_back("u_lon_tor");
        vars_to_write.push_back("u_lat_tor");

        vars_to_write.push_back("u_lon_pot");
        vars_to_write.push_back("u_lat_pot");

        vars_to_write.push_back("vorticity");
        vars_to_write.push_back("divergence");

        vars_to_write.push_back("proj_vorticity");
        vars_to_write.push_back("proj_divergence");
    }

    vars_to_write.push_back("Psi");
    vars_to_write.push_back("Phi");

    initialize_output_file( source_data, vars_to_write, output_fname.c_str(), -1);

    if (not(constants::MINIMAL_OUTPUT)) {
        write_field_to_output(full_u_lon_tor,  "u_lon_tor",  starts, counts, output_fname.c_str(), &unmask);
        write_field_to_output(full_u_lat_tor,  "u_lat_tor",  starts, counts, output_fname.c_str(), &unmask);

        write_field_to_output(full_u_lon_pot,  "u_lon_pot",  starts, counts, output_fname.c_str(), &unmask);
        write_field_to_output(full_u_lat_pot,  "u_lat_pot",  starts, counts, output_fname.c_str(), &unmask);

        write_field_to_output(vort_term,  "vorticity",  starts, counts, output_fname.c_str(), &unmask);
        write_field_to_output(div_term,   "divergence", starts, counts, output_fname.c_str(), &unmask);

        toroidal_vel_div(        div_term, full_u_lon_pot, full_u_lat_pot, source_data, use_mask ? mask : unmask );
        toroidal_curl_u_dot_er( vort_term, full_u_lon_tor, full_u_lat_tor, source_data, use_mask ? mask : unmask );

        write_field_to_output(vort_term,  "proj_vorticity",  starts, counts, output_fname.c_str(), &unmask);
        write_field_to_output(div_term,   "proj_divergence", starts, counts, output_fname.c_str(), &unmask);
    }

    write_field_to_output(full_Psi, "Psi", starts, counts, output_fname.c_str(), &unmask);
    write_field_to_output(full_Phi, "Phi", starts, counts, output_fname.c_str(), &unmask);

    // Store some solver information
    add_attr_to_file("rel_tol",         rel_tol,                        output_fname.c_str());
    add_attr_to_file("max_iters",       (double) max_iters,             output_fname.c_str());
    add_attr_to_file("diff_order",      (double) constants::DiffOrd,    output_fname.c_str());
    add_attr_to_file("use_mask",        (double) use_mask,              output_fname.c_str());
    add_attr_to_file("weight_err",      (double) weight_err,            output_fname.c_str());
    add_attr_to_file("Tikhov_Laplace",  Tikhov_Laplace,                 output_fname.c_str());


    //
    //// At the very end, compute the L2 and LInf error for each time/depth
    //

    #if DEBUG >= 1
    if (wRank == 0) {
        fprintf(stdout, "Computing the error of the projection.\n");
    }
    #endif

    std::vector<double> projection_2error(      Ntime * Ndepth, 0. ),
                        projection_Inferror(    Ntime * Ndepth, 0. ),
                        velocity_Infnorm(       Ntime * Ndepth, 0. ),
                        projection_KE(          Ntime * Ndepth, 0. ),
                        toroidal_KE(            Ntime * Ndepth, 0. ),
                        potential_KE(           Ntime * Ndepth, 0. ),
                        velocity_2norm(         Ntime * Ndepth, 0. ),
                        tot_areas(              Ntime * Ndepth, 0. );
    double total_area, error2, errorInf, velInf, tor_KE, pot_KE, proj_KE, orig_KE;
    for (int Itime = 0; Itime < Ntime; ++Itime) {
        for (int Idepth = 0; Idepth < Ndepth; ++Idepth) {

            total_area = 0.;
            error2 = 0.;
            tor_KE = 0.;
            pot_KE = 0.;
            proj_KE = 0.;
            orig_KE = 0.;
            errorInf = 0.;
            velInf = 0.;

            #pragma omp parallel \
            default(none) \
            shared( full_u_lon_tor, full_u_lat_tor, full_u_lon_pot, full_u_lat_pot, \
                    u_lon, u_lat, Itime, Idepth, dAreas, latitude ) \
            reduction(+ : total_area, error2, tor_KE, pot_KE, proj_KE, orig_KE) \
            reduction( max : errorInf, velInf )\
            private( index, index_sub ) \
            firstprivate( Ndepth, Ntime, Npts )
            {
                #pragma omp for collapse(1) schedule(static)
                for (index_sub = 0; index_sub < Npts; ++index_sub) {
                    index = index_sub + Npts*(Itime*Ndepth + Idepth);

                    total_area += dAreas.at(index_sub);

                    error2 += dAreas.at(index_sub) * (
                                    pow( u_lon.at(index) - full_u_lon_tor.at(index) - full_u_lon_pot.at(index) , 2.)
                                 +  pow( u_lat.at(index) - full_u_lat_tor.at(index) - full_u_lat_pot.at(index) , 2.)
                            );

                    errorInf = std::fmax( 
                                    errorInf,
                                    sqrt(     pow( u_lon.at(index) - full_u_lon_tor.at(index) - full_u_lon_pot.at(index) , 2.)
                                           +  pow( u_lat.at(index) - full_u_lat_tor.at(index) - full_u_lat_pot.at(index) , 2.)
                                         )
                                    );

                    velInf = std::fmax( velInf,  std::fabs( sqrt( pow( u_lon.at(index) , 2.) +  pow( u_lat.at(index) , 2.) ) )  );

                    tor_KE += dAreas.at(index_sub) * ( pow( full_u_lon_tor.at(index), 2.) + pow( full_u_lat_tor.at(index), 2.) );
                    pot_KE += dAreas.at(index_sub) * ( pow( full_u_lon_pot.at(index), 2.) + pow( full_u_lat_pot.at(index), 2.) );

                    proj_KE += dAreas.at(index_sub) * (
                                    pow( full_u_lon_tor.at(index) + full_u_lon_pot.at(index) , 2.)
                                 +  pow( full_u_lat_tor.at(index) + full_u_lat_pot.at(index) , 2.)
                            );

                    orig_KE += dAreas.at(index_sub) * ( pow( u_lon.at(index), 2.) + pow( u_lat.at(index), 2.) );
                }
            }
            size_t int_index = Index( Itime, Idepth, 0, 0, Ntime, Ndepth, 1, 1);

            tot_areas.at(int_index) = total_area;

            projection_2error.at(   int_index ) = sqrt( error2   / total_area );
            projection_Inferror.at( int_index ) = errorInf;

            velocity_2norm.at(   int_index ) = sqrt( orig_KE  / total_area );
            velocity_Infnorm.at( int_index ) = velInf;

            projection_KE.at( int_index ) = sqrt( proj_KE  / total_area );
            toroidal_KE.at(   int_index ) = sqrt( tor_KE   / total_area );
            potential_KE.at(  int_index ) = sqrt( pot_KE   / total_area );
        }
    }

    const char* dim_names[] = {"time", "depth"};
    const int ndims_error = 2;
    if (wRank == 0) {
        add_var_to_file( "total_area",    dim_names, ndims_error, output_fname.c_str() );

        add_var_to_file( "projection_2error",    dim_names, ndims_error, output_fname.c_str() );
        add_var_to_file( "projection_Inferror",  dim_names, ndims_error, output_fname.c_str() );

        add_var_to_file( "velocity_2norm",   dim_names, ndims_error, output_fname.c_str() );
        add_var_to_file( "velocity_Infnorm", dim_names, ndims_error, output_fname.c_str() );

        add_var_to_file( "projection_KE",  dim_names, ndims_error, output_fname.c_str() );
        add_var_to_file( "toroidal_KE",    dim_names, ndims_error, output_fname.c_str() );
        add_var_to_file( "potential_KE",   dim_names, ndims_error, output_fname.c_str() );
    }
    MPI_Barrier(MPI_COMM_WORLD);

    size_t starts_error[ndims_error] = { size_t(myStarts.at(0)), size_t(myStarts.at(1)) };
    size_t counts_error[ndims_error] = { size_t(Ntime), size_t(Ndepth) };

    write_field_to_output( tot_areas,   "total_area",   starts_error, counts_error, output_fname.c_str() );

    write_field_to_output( projection_2error,   "projection_2error",   starts_error, counts_error, output_fname.c_str() );
    write_field_to_output( projection_Inferror, "projection_Inferror", starts_error, counts_error, output_fname.c_str() );

    write_field_to_output( velocity_2norm,   "velocity_2norm",   starts_error, counts_error, output_fname.c_str() );
    write_field_to_output( velocity_Infnorm, "velocity_Infnorm", starts_error, counts_error, output_fname.c_str() );

    write_field_to_output( projection_KE, "projection_KE", starts_error, counts_error, output_fname.c_str() );
    write_field_to_output( toroidal_KE,   "toroidal_KE",   starts_error, counts_error, output_fname.c_str() );
    write_field_to_output( potential_KE,  "potential_KE",  starts_error, counts_error, output_fname.c_str() );

}
