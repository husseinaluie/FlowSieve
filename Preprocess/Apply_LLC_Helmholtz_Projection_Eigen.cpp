#include "../constants.hpp"
#include "../functions.hpp"
#include "../netcdf_io.hpp"
#include "../preprocess.hpp"
#include "../differentiation_tools.hpp"
#include <algorithm>
#include <vector>
#include <omp.h>
#include <math.h>
//#include <eigen3/Eigen/Sparse>
//#include <eigen3/Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>


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
        const MPI_Comm comm
        ) {

    int wRank, wSize;
    MPI_Comm_rank( comm, &wRank );
    MPI_Comm_size( comm, &wSize );

    // Create some tidy names for variables
    const std::vector<double>   &latitude   = source_data.latitude,
                                &longitude  = source_data.longitude,
                                &dAreas     = source_data.areas;

    const std::vector<short int> &mask = (constants::FILTER_OVER_LAND) ? source_data.reference_mask : source_data.mask;

    const std::vector<int>  &myCounts = source_data.myCounts,
                            &myStarts = source_data.myStarts;

    std::vector<double>   &u_lat = source_data.variables.at("u_lat"),
                          &u_lon = source_data.variables.at("u_lon");

    // Create a 'no mask' mask variable
    //   we'll treat land values as zero velocity
    //   We do this because including land seems
    //   to introduce strong numerical issues
    const std::vector<short int> unmask(mask.size(), true);

    const int   Ntime   = myCounts.at(0),
                Ndepth  = myCounts.at(1);

    const size_t Npts = latitude.size();

    int Itime=0, Idepth=0;
    size_t index, index_sub, iters_used = 0;

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

    size_t num_land_points = 0;
    if (constants::FILTER_OVER_LAND) {
        for (index = 0; index < u_lon.size(); index++) {
            if (not(mask.at(index))) { num_land_points++; }
        }
    }
    #if DEBUG>=0
    if (wRank == 0) { fprintf(stdout, " Identified %'zu land points.\n", num_land_points); }
    #endif

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
    assert( Tikhov_Laplace < 0 );
    const size_t Nboxrows = ( Tikhov_Laplace > 0 ) ? 4 : 2;
    std::vector<double> 
        //RHS_vector( Npts + 2*num_land_points, 0. ),
        RHS_vector( Npts + constants::ADJACENCY_SIZE*num_land_points, 0. ),
        Psi_seed(       Npts, 0. ),
        Phi_seed(       Npts, 0. ),
        work_arr(       Npts, 0. ),
        div_term(       Npts, 0. ),
        vort_term(      Npts, 0. ),
        u_lon_rem(      Npts, 0. ),
        u_lat_rem(      Npts, 0. );
    

    fprintf( stdout, "Copy seed\n" );
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

    const size_t num_neighbours = source_data.num_neighbours;
    // Get a magnitude for the derivatives, to help normalize the rows of the 
    //  Laplace entries to have similar magnitude to the others.
    long double deriv_ref_1 = 0, deriv_ref_2 = 0;
    size_t Ineighbour;
    #pragma omp parallel default(none) \
    private( Ineighbour, index ) \
    shared( source_data ) \
    firstprivate( num_neighbours, Npts ) \
    reduction( +:deriv_ref_1,deriv_ref_2 )
    {
        #pragma omp for collapse(1) schedule(static)
        for (index = 0; index < Npts; ++index) {
            for ( Ineighbour = 0; Ineighbour < num_neighbours + 1; Ineighbour++ ) {
                deriv_ref_1 += std::fabs( source_data.adjacency_ddlat_weights.at(index).at(Ineighbour) ) / (Npts*num_neighbours);
                //deriv_ref_2 += std::fabs( source_data.adjacency_d2dlat2_weights.at(index).at(Ineighbour) ) / (Npts*num_neighbours);
            }
        }
    }
    //const double deriv_scale_factor = deriv_ref_2 / deriv_ref_1;
    const double deriv_scale_factor = deriv_ref_1;
    fprintf( stdout, "deriv-scale-factor: %g\n", deriv_scale_factor );

    //
    //// Build the LHS part of the problem
    //      Ordering is: [  u_from_psi      u_from_phi   ] *  [ psi ]   =    [  u   ]
    //                   [  v_from_psi      v_from_phi   ]    [ phi ]        [  v   ]
    //
    //      Ordering is: [           - ddlat   sec(phi) * ddlon   ] *  [ psi ]   =    [     u     ]
    //                   [  sec(phi) * ddlon              ddlat   ]    [ phi ]        [     v     ]
    //                   [           Laplace                  0   ]                   [ vort(u,v) ]
    //                   [                 0            Laplace   ]                   [  div(u,v) ]
    
    #if DEBUG >= 1
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
    //std::vector<T> Aij_triplets( Npts * (2 * pts_per_2nd_deriv + pts_per_1st_deriv) + 2 * num_land_points * pts_per_1st_deriv );
    std::vector<T> Aij_triplets( Npts * (2 * pts_per_2nd_deriv + pts_per_1st_deriv) 
                                  + 2 * num_land_points * constants::ADJACENCY_SIZE );
    size_t Itriplet, Ipt, neighbour_ind;
    double weight_val, cos_lat_inv, R_inv;
    bool is_pole, all_land_neighbours;
    #pragma omp parallel default(none) \
    shared( mask, dAreas, latitude, source_data, Aij_triplets ) \
    private( Ipt, Ineighbour, neighbour_ind, Itriplet, row_skip, column_skip, is_pole, val, \
             weight_val, cos_lat_inv, R_inv, all_land_neighbours ) \
    firstprivate( Npts, weight_err, num_neighbours, Tikhov_Laplace, \
                  stdout, deriv_scale_factor, pts_per_1st_deriv, pts_per_2nd_deriv )
    { 
        #pragma omp for collapse(1) schedule(static)
        for ( Ipt = 0; Ipt < Npts; Ipt++ ) {

            weight_val = weight_err ? dAreas.at(Ipt) : 1.;
            cos_lat_inv = 1. / cos(latitude.at(Ipt));
            R_inv = 1. / constants::R_earth;

            if ( not(mask[Ipt]) ) {
                // If it's land, instead of solving the Lapacian, just do the
                // equal-to-neighbours bit

                // Check that all neighbours are land
                all_land_neighbours = true;
                for ( Ineighbour = 0; Ineighbour < num_neighbours; Ineighbour++ ) {
                    neighbour_ind = source_data.adjacency_indices.at(Ipt).at(Ineighbour);
                    if ( mask[neighbour_ind] ) { all_land_neighbours = false; }
                }

                if ( all_land_neighbours ) {
                    // all neighbours are land, so just take the first
                    neighbour_ind = source_data.adjacency_indices.at(Ipt).at(0);
                    column_skip = neighbour_ind;
                    row_skip    = Ipt;
                    if ( USE_TRUE_2ND_DERIV ) {
                        Itriplet = Ipt * pts_per_2nd_deriv + 0;
                    } else {
                        Itriplet = Ipt * pts_per_2nd_deriv + 0*(num_neighbours+1) + 0;
                    }
                    Aij_triplets[Itriplet  ] = T( row_skip, Ipt,            weight_val * Tikhov_Laplace );
                    Aij_triplets[Itriplet+1] = T( row_skip, neighbour_ind, -weight_val * Tikhov_Laplace );
                    continue;
                }
            }

            for ( Ineighbour = 0; Ineighbour < num_neighbours + 1; Ineighbour++ ) {

                neighbour_ind = (Ineighbour < num_neighbours) ? 
                                        source_data.adjacency_indices.at(Ipt).at(Ineighbour) :
                                        Ipt;

                is_pole = std::fabs( std::fabs( latitude.at(Ipt) * 180.0 / M_PI ) - 90 ) < 1e-6;
                if ( is_pole ) { fprintf(stdout, "SKIPPING POLE POINT!\n"); continue; }

                //
                //// Second LON derivative
                //

                if ( USE_TRUE_2ND_DERIV ) {
                    val  = source_data.adjacency_d2dlon2_weights.at(Ipt).at(Ineighbour);
                    val *= weight_val * pow(cos_lat_inv * R_inv, 2.);

                    column_skip = neighbour_ind;
                    row_skip    = Ipt;
                    Itriplet = Ipt * pts_per_2nd_deriv + Ineighbour;
                    Aij_triplets[Itriplet] = T( row_skip, column_skip, val );
                } else {
                    for ( size_t D2_ind = 0; D2_ind < num_neighbours+1; D2_ind++ ) {
                        val  =   source_data.adjacency_ddlon_weights.at(Ipt).at(Ineighbour)
                               * source_data.adjacency_ddlon_weights.at(neighbour_ind).at(D2_ind);
                        val *= weight_val * pow(R_inv, 2.) * cos_lat_inv / cos(latitude.at(neighbour_ind));

                        column_skip = source_data.adjacency_indices.at(neighbour_ind).at(D2_ind);
                        row_skip    = Ipt;
                        Itriplet = Ipt * pts_per_2nd_deriv + Ineighbour*(num_neighbours+1)  + D2_ind;
                        Aij_triplets[Itriplet] = T( row_skip, column_skip, val );
                    }
                }
    
                //
                //// Second LAT derivative
                //
    
                if ( USE_TRUE_2ND_DERIV ) {
                    val = source_data.adjacency_d2dlat2_weights.at(Ipt).at(Ineighbour);
                    val *= weight_val * pow(R_inv, 2.);

                    column_skip = neighbour_ind;
                    row_skip    = Ipt;
                    Itriplet    = (Npts + Ipt) * pts_per_2nd_deriv + Ineighbour;
                    Aij_triplets[Itriplet] = T( row_skip, column_skip, val );
                } else {
                    for ( size_t D2_ind = 0; D2_ind < num_neighbours+1; D2_ind++ ) {
                        val  =   source_data.adjacency_ddlat_weights.at(Ipt).at(Ineighbour)
                               * source_data.adjacency_ddlat_weights.at(neighbour_ind).at(D2_ind);
                        val *= weight_val * pow(R_inv, 2.);

                        column_skip = source_data.adjacency_indices.at(neighbour_ind).at(D2_ind);
                        row_skip    = Ipt;
                        Itriplet = (Npts + Ipt) * pts_per_2nd_deriv + Ineighbour*(num_neighbours+1)  + D2_ind;
                        Aij_triplets[Itriplet] = T( row_skip, column_skip, val );
                    }
                }

                //
                //// First LAT derivative
                //
    
                val = - source_data.adjacency_ddlat_weights.at(Ipt).at(Ineighbour) * tan( latitude.at(Ipt) );
                val *= weight_val * pow(R_inv, 2.);
    
                column_skip = neighbour_ind;
                row_skip    = Ipt;
                Itriplet    = 2 * Npts * pts_per_2nd_deriv + Ipt * pts_per_1st_deriv + Ineighbour;
                Aij_triplets[Itriplet] = T( row_skip, column_skip, val );

            }
        }
    }

    size_t Iland, count, row, column;

    #pragma omp parallel default(none) \
    shared( source_data, Aij_triplets, mask, dAreas, latitude ) \
    private( count, Iland, Ipt, Ineighbour, row, column, val, \
             weight_val, cos_lat_inv, R_inv, neighbour_ind, Itriplet, \
             all_land_neighbours ) \
    firstprivate( num_land_points, Npts, num_neighbours, weight_err,\
                  Tikhov_Laplace )
    {
        count = 0;
        #pragma omp for collapse(1) schedule(guided)
        for ( Iland = 0; Iland < num_land_points; Iland++ ) {
            for ( Ipt = 0; Ipt < Npts; Ipt++ ) {
                if ( not(mask[Ipt]) ) {
                    if ( Iland == count ) {

                        weight_val = weight_err ? dAreas.at(Ipt) : 1.;
                        cos_lat_inv = 1. / cos(latitude.at(Ipt));
                        R_inv = 1. / constants::R_earth;

                        // Check that all neighbours are land
                        all_land_neighbours = true;
                        for ( Ineighbour = 0; Ineighbour < num_neighbours; Ineighbour++ ) {
                            neighbour_ind = source_data.adjacency_indices.at(Ipt).at(Ineighbour);
                            if ( mask[neighbour_ind] ) { all_land_neighbours = false; }
                        }
                        if ( not(all_land_neighbours) ) { continue; }

                        // If they are, enforce equality between all neighbours (and centre)
                        for ( Ineighbour = 0; Ineighbour < num_neighbours; Ineighbour++ ) {
                            neighbour_ind = source_data.adjacency_indices.at(Ipt).at(Ineighbour);
                            row = Npts + Iland * constants::ADJACENCY_SIZE + Ineighbour;
                            Itriplet = Npts * (2 * pts_per_2nd_deriv + pts_per_1st_deriv) 
                                        + 2 * (Iland * constants::ADJACENCY_SIZE + Ineighbour);
                            Aij_triplets[Itriplet  ] = T( row, Ipt,            weight_val * Tikhov_Laplace );
                            Aij_triplets[Itriplet+1] = T( row, neighbour_ind, -weight_val * Tikhov_Laplace );
                        }
                    }
                    count++;
                }
            }
        }
    }





    //Eigen::SparseMatrix<double> LHS_matr( Npts + 2*num_land_points, Npts );
    Eigen::SparseMatrix<double> LHS_matr( Npts + constants::ADJACENCY_SIZE*num_land_points, Npts );
    LHS_matr.setFromTriplets( Aij_triplets.begin(), Aij_triplets.end() );

    fprintf(stdout, "  Land counter: %'zu\n", land_counter);

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
    solver.compute( LHS_matr );
    if ( solver.info() != Eigen::Success ) {
        // decomposition failed
        fprintf( stderr, "!!Eigen decomposition failed.\n" );
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

            #if DEBUG >= 2
            if ( wRank == 0 ) {
                fprintf(stdout, "Building the RHS of the least squares problem.\n");
                fflush(stdout);
            }
            #endif

            //
            //// Set up the RHS_vector
            //
            
            double is_pole;
            std::fill( RHS_vector.begin(), RHS_vector.end(), 0. );
            #pragma omp parallel default(none) \
            shared( dAreas, RHS_vector, vort_term ) \
            private( index_sub ) firstprivate( weight_err, Npts )
            {
                #pragma omp for collapse(1) schedule(static)
                for (index_sub = 0; index_sub < Npts; ++index_sub) {
                    RHS_vector.at( index_sub) = vort_term.at(index_sub) * ( weight_err ? dAreas.at(index_sub) : 1. );
                }
            }

            /*
            fprintf( stdout, "Also adding the seed info over land, %zu\n", Nboxrows );
            land_counter = 0;
            for (index_sub = 0; index_sub < Npts; ++index_sub) {
                if (not(mask.at(index_sub))) {
                    double weight_val = weight_err ? dAreas.at(index_sub) : 1.;

                    row_skip    = Nboxrows*Npts + 0*num_land_points + land_counter;
                    //RHS_vector.at( row_skip ) = -dPsi_dlon;
                    RHS_vector.at( row_skip ) = -u_lat_tor_seed[index_sub] * weight_val;

                    row_skip    = Nboxrows*Npts + 1*num_land_points + land_counter;
                    //RHS_vector.at( row_skip ) = -dPhi_dlon;
                    RHS_vector.at( row_skip ) = -u_lon_pot_seed[index_sub] * weight_val;

                    row_skip    = Nboxrows*Npts + 2*num_land_points + land_counter;
                    //RHS_vector.at( row_skip ) = -dPsi_dlat;
                    RHS_vector.at( row_skip ) = -u_lon_tor_seed[index_sub] * weight_val;

                    row_skip    = Nboxrows*Npts + 3*num_land_points + land_counter;
                    //RHS_vector.at( row_skip ) = -dPhi_dlat;
                    RHS_vector.at( row_skip ) = -u_lat_pot_seed[index_sub] * weight_val;

                    land_counter++;
                }
            }
            fprintf( stdout, "%zu land points\n", land_counter );
            */
            Eigen::VectorXd RHS = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(RHS_vector.data(), RHS_vector.size());

            //
            //// Now apply the least-squares solver
            //
            #if DEBUG >= 2
            if ( wRank == 0 ) {
                fprintf(stdout, "Solving the least squares problem for Psi.\n");
                fflush(stdout);
            }
            #endif
            Eigen::VectorXd F_Eigen = solver.solve( RHS );
            #if DEBUG >= 2
            if ( wRank == 0 ) {
                fprintf( stdout, "    Solver converged after %ld iterations to error %g.\n", 
                        solver.iterations(), solver.error() );
                fflush(stdout);
            }
            #endif
            std::vector<double> Psi_vector(F_Eigen.data(), F_Eigen.data() + Npts);

            //
            //// And repeat for divergence!
            //
            std::fill( RHS_vector.begin(), RHS_vector.end(), 0. );
            #pragma omp parallel default(none) \
            shared( dAreas, RHS_vector, div_term ) \
            private( index_sub ) firstprivate( weight_err, Npts )
            {
                #pragma omp for collapse(1) schedule(static)
                for (index_sub = 0; index_sub < Npts; ++index_sub) {
                    RHS_vector.at( index_sub) = div_term.at(index_sub) * ( weight_err ? dAreas.at(index_sub) : 1. );
                }
            }
            RHS = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(RHS_vector.data(), RHS_vector.size());

            #if DEBUG >= 2
            if ( wRank == 0 ) {
                fprintf(stdout, "Solving the least squares problem for Phi.\n");
                fflush(stdout);
            }
            #endif
            F_Eigen = solver.solve( RHS );
            #if DEBUG >= 2
            if ( wRank == 0 ) {
                fprintf( stdout, "    Solver converged after %ld iterations to error %g.\n", 
                        solver.iterations(), solver.error() );
                fflush(stdout);
            }
            #endif
            std::vector<double> Phi_vector(F_Eigen.data(), F_Eigen.data() + Npts);
            

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
