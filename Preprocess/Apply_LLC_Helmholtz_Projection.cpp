#include "../constants.hpp"
#include "../functions.hpp"
#include "../netcdf_io.hpp"
#include "../preprocess.hpp"
#include "../differentiation_tools.hpp"
#include <algorithm>
#include <vector>
#include <omp.h>
#include <math.h>
#include "../ALGLIB/stdafx.h"
#include "../ALGLIB/linalg.h"
#include "../ALGLIB/solvers.h"

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
        const MPI_Comm comm
        ) {

    int wRank, wSize;
    MPI_Comm_rank( comm, &wRank );
    MPI_Comm_size( comm, &wSize );

    // Create some tidy names for variables
    const std::vector<double>   &latitude   = source_data.latitude,
                                &longitude  = source_data.longitude,
                                &dAreas     = source_data.areas;

    const std::vector<bool> &mask = source_data.mask;

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
    const size_t Nboxrows = ( Tikhov_Laplace > 0 ) ? 4 : 2;
    alglib::real_1d_array rhs;//, proj_vals;
    std::vector<double> 
        RHS_vector( Nboxrows * Npts, 0. ),
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
    double deriv_ref_1 = 0, deriv_ref_2 = 0;
    for ( size_t Ineighbour = 0; Ineighbour < num_neighbours + 1; Ineighbour++ ) {
        deriv_ref_1 += std::fabs( source_data.adjacency_ddlat_weights.at(0).at(Ineighbour) );
        deriv_ref_2 += std::fabs( source_data.adjacency_d2dlat2_weights.at(0).at(Ineighbour) );
    }
    //const double deriv_scale_factor = deriv_ref_2 / deriv_ref_1;
    const double deriv_scale_factor = deriv_ref_1;
    fprintf( stdout, "deriv-scale-factor: %g\n", deriv_scale_factor );


    rhs.attach_to_ptr( Nboxrows * Npts, &RHS_vector[0] );

    alglib::linlsqrstate state;
    alglib::linlsqrreport report;

    alglib::real_1d_array F_alglib;

    double *F_array;

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

    alglib::sparsematrix LHS_matr;
    alglib::sparsecreate(Nboxrows*Npts, 2*Npts, LHS_matr);
    for ( size_t Ipt = 0; Ipt < Npts; Ipt++ ) {
        //fprintf( stdout, " - %'zu \n", Ipt );

        double  weight_val = weight_err ? dAreas.at(Ipt) : 1.,
                cos_lat_inv = 1. / cos(latitude.at(Ipt)),
                R_inv = 1. / constants::R_earth;

        for ( size_t Ineighbour = 0; Ineighbour < num_neighbours + 1; Ineighbour++ ) {
            //fprintf( stdout, " - - %'zu \n", Ineighbour );

            size_t neighbour_ind = (Ineighbour < num_neighbours) ? 
                                    source_data.adjacency_indices.at(Ipt).at(Ineighbour) :
                                    Ipt;

            bool is_pole = std::fabs( std::fabs( latitude.at(Ipt) * 180.0 / M_PI ) - 90 ) < 0.001;
            if ( is_pole ) { continue; }

            //
            //// LON first derivative
            //

            double val  = source_data.adjacency_ddlon_weights.at(Ipt).at(Ineighbour);
            val *= weight_val * cos_lat_inv * R_inv;
            //fprintf( stdout, " - - %g \n", val );

            // Psi part
            size_t  column_skip = 0 * Npts + neighbour_ind,
                    row_skip    = 1 * Npts + Ipt;
            alglib::sparseadd( LHS_matr, row_skip, column_skip, val );

            // Phi part
            column_skip = 1 * Npts + neighbour_ind,
            row_skip    = 0 * Npts + Ipt;
            alglib::sparseadd( LHS_matr, row_skip, column_skip, val );


            //
            //// LAT first derivative
            //

            val  = source_data.adjacency_ddlat_weights.at(Ipt).at(Ineighbour);
            val *= weight_val * R_inv;
            //fprintf( stdout, " - - %g \n", val );

            // Psi part
            column_skip = 0 * Npts + neighbour_ind,
            row_skip    = 0 * Npts + Ipt;
            alglib::sparseadd( LHS_matr, row_skip, column_skip, -val );

            // Phi part
            column_skip = 1 * Npts + neighbour_ind,
            row_skip    = 1 * Npts + Ipt;
            alglib::sparseadd( LHS_matr, row_skip, column_skip, val );


            if (Tikhov_Laplace > 0) {

                //
                //// Second LON derivative
                //
                //  NOTE:   We're doing this as the first derivative of the first
                //          derivative. Typically, this is not great, since it
                //          tends to have a higher error. However, in our case
                //          the broader stencil is helpful. It's also more internally
                //          consistent since later we take derivs( vel ) = deriv( deriv ( Helm ) )

                for ( size_t D2_ind = 0; D2_ind < num_neighbours+1; D2_ind++ ) {
                    val  =   source_data.adjacency_ddlon_weights.at(Ipt).at(Ineighbour)
                           * source_data.adjacency_ddlon_weights.at(neighbour_ind).at(D2_ind);
                    val *= weight_val * pow(R_inv, 2.) * cos_lat_inv / cos(latitude.at(neighbour_ind));
                    val *= Tikhov_Laplace / deriv_scale_factor;

                    column_skip = 0 * Npts + source_data.adjacency_indices.at(neighbour_ind).at(D2_ind);
                    row_skip    = 2 * Npts + Ipt;
                    alglib::sparseadd( LHS_matr, row_skip, column_skip, val );

                    column_skip = 1 * Npts + source_data.adjacency_indices.at(neighbour_ind).at(D2_ind);
                    row_skip    = 3 * Npts + Ipt;
                    alglib::sparseadd( LHS_matr, row_skip, column_skip, val );
                }

                /*
                // Version using actual second derivative
                val  = source_data.adjacency_d2dlon2_weights.at(Ipt).at(Ineighbour);
                val *= weight_val * pow(cos_lat_inv * R_inv, 2.);
                val *= Tikhov_Laplace / deriv_scale_factor;

                column_skip = 0 * Npts + neighbour_ind;
                row_skip    = 2 * Npts + Ipt;
                alglib::sparseadd( LHS_matr, row_skip, column_skip, val );

                column_skip = 1 * Npts + neighbour_ind;
                row_skip    = 3 * Npts + Ipt;
                alglib::sparseadd( LHS_matr, row_skip, column_skip, val );
                */


                //
                //// Second LAT derivative
                //
                for ( size_t D2_ind = 0; D2_ind < num_neighbours+1; D2_ind++ ) {
                    val  =   source_data.adjacency_ddlat_weights.at(Ipt).at(Ineighbour)
                           * source_data.adjacency_ddlat_weights.at(neighbour_ind).at(D2_ind);
                    val *= weight_val * pow(R_inv, 2.);
                    val *= Tikhov_Laplace / deriv_scale_factor;

                    column_skip = 0 * Npts + source_data.adjacency_indices.at(neighbour_ind).at(D2_ind);
                    row_skip    = 2 * Npts + Ipt;
                    alglib::sparseadd( LHS_matr, row_skip, column_skip, val );

                    column_skip = 1 * Npts + source_data.adjacency_indices.at(neighbour_ind).at(D2_ind);
                    row_skip    = 3 * Npts + Ipt;
                    alglib::sparseadd( LHS_matr, row_skip, column_skip, val );
                }

                /*
                // Version using actual second derivative
                val  = source_data.adjacency_d2dlat2_weights.at(Ipt).at(Ineighbour);
                val *= weight_val * pow(R_inv, 2.);
                val *= Tikhov_Laplace / deriv_scale_factor;

                column_skip = 0 * Npts + neighbour_ind;
                row_skip    = 2 * Npts + Ipt;
                alglib::sparseadd( LHS_matr, row_skip, column_skip, val );

                column_skip = 1 * Npts + neighbour_ind;
                row_skip    = 3 * Npts + Ipt;
                alglib::sparseadd( LHS_matr, row_skip, column_skip, val );
                */

                //
                //// LAT first derivative
                //
                val  = source_data.adjacency_ddlat_weights.at(Ipt).at(Ineighbour);
                val *= weight_val * pow(R_inv, 2.) * tan( latitude.at(Ipt) );
                val *= Tikhov_Laplace / deriv_scale_factor;

                // Psi part
                column_skip = 0 * Npts + neighbour_ind,
                row_skip    = 2 * Npts + Ipt;
                alglib::sparseadd( LHS_matr, row_skip, column_skip, -val );

                // Phi part
                column_skip = 1 * Npts + neighbour_ind,
                row_skip    = 3 * Npts + Ipt;
                alglib::sparseadd( LHS_matr, row_skip, column_skip, -val );

            }

        }
    }

    alglib::sparseconverttocrs(LHS_matr);

    #if DEBUG >= 1
    if (wRank == 0) {
        fprintf(stdout, "Declaring the least squares problem.\n");
        fflush(stdout);
    }
    #endif
    alglib::linlsqrcreate(Nboxrows*Npts, 2*Npts, state);
    alglib::linlsqrsetcond(state, rel_tol, rel_tol, max_iters);

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
            shared( Itime, Idepth, RHS_vector, stdout, \
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
            #pragma omp parallel default(none) \
            shared( dAreas, latitude, Itime, Idepth, RHS_vector, div_term, vort_term, u_lon_rem, u_lat_rem ) \
            private( index, index_sub, is_pole ) \
            firstprivate( Ndepth, Ntime, Npts, Tikhov_Laplace, weight_err, deriv_scale_factor, stdout )
            {
                #pragma omp for collapse(1) schedule(static)
                for (index_sub = 0; index_sub < Npts; ++index_sub) {
                    index = index_sub + Npts*(Itime*Ndepth + Idepth);

                    is_pole = std::fabs( std::fabs( latitude.at(index_sub) * 180.0 / M_PI ) - 90 ) < 0.001;

                    RHS_vector.at( 0*Npts + index_sub) = u_lon_rem.at(index_sub);
                    RHS_vector.at( 1*Npts + index_sub) = u_lat_rem.at(index_sub);

                    if (Tikhov_Laplace > 0) {
                        if ( is_pole ) {
                            RHS_vector.at( 2*Npts + index_sub) = 0.;
                            RHS_vector.at( 3*Npts + index_sub) = 0.;
                        } else {
                            RHS_vector.at( 2*Npts + index_sub) = vort_term.at(index_sub) * Tikhov_Laplace / deriv_scale_factor;
                            RHS_vector.at( 3*Npts + index_sub) = div_term.at( index_sub) * Tikhov_Laplace / deriv_scale_factor;
                        }
                    }

                    if ( weight_err ) {
                        RHS_vector.at( 0*Npts + index_sub) *= dAreas.at(index_sub);
                        RHS_vector.at( 1*Npts + index_sub) *= dAreas.at(index_sub);
                        if (Tikhov_Laplace > 0) {
                            RHS_vector.at( 2*Npts + index_sub) *= dAreas.at(index_sub);
                            RHS_vector.at( 3*Npts + index_sub) *= dAreas.at(index_sub);
                        }
                    }

                    //fprintf( stdout, " %'zu : %.3g %.3g %.3g %.3g\n", index, 
                    //        RHS_vector.at( 0*Npts + index_sub), 
                    //        RHS_vector.at( 1*Npts + index_sub), 
                    //        RHS_vector.at( 2*Npts + index_sub), 
                    //        RHS_vector.at( 3*Npts + index_sub) );
                }
            }

            //
            //// Now apply the least-squares solver
            //
            #if DEBUG >= 2
            if ( wRank == 0 ) {
                fprintf(stdout, "Solving the least squares problem.\n");
                fflush(stdout);
            }
            #endif
            alglib::linlsqrsolvesparse(state, LHS_matr, rhs);
            alglib::linlsqrresults(state, F_alglib, report);

            /*    Rep     -   optimization report:
                * Rep.TerminationType completetion code:
                    *  1    ||Rk||<=EpsB*||B||
                    *  4    ||A^T*Rk||/(||A||*||Rk||)<=EpsA
                    *  5    MaxIts steps was taken
                    *  7    rounding errors prevent further progress,
                            X contains best point found so far.
                            (sometimes returned on singular systems)
                    *  8    user requested termination via calling
                            linlsqrrequesttermination()
                * Rep.IterationsCount contains iterations count
                * NMV countains number of matrix-vector calculations
            */

            #if DEBUG >= 1
            if      (report.terminationtype == 1) { fprintf(stdout, "Termination type: absolulte tolerance reached.\n"); }
            else if (report.terminationtype == 4) { fprintf(stdout, "Termination type: relative tolerance reached.\n"); }
            else if (report.terminationtype == 5) { fprintf(stdout, "Termination type: maximum number of iterations reached.\n"); }
            else if (report.terminationtype == 7) { fprintf(stdout, "Termination type: round-off errors prevent further progress.\n"); }
            else if (report.terminationtype == 8) { fprintf(stdout, "Termination type: user requested (?)\n"); }
            else                                  { fprintf(stdout, "Termination type: unknown\n"); }
            #endif
            if      (report.terminationtype == 1) { terminate_count_abs_tol++; }
            else if (report.terminationtype == 4) { terminate_count_rel_tol++; }
            else if (report.terminationtype == 5) { terminate_count_max_iter++; }
            else if (report.terminationtype == 7) { terminate_count_rounding++; }
            else if (report.terminationtype == 8) { terminate_count_other++; }
            else                                  { terminate_count_other++; }

            iters_used = linlsqrpeekiterationscount( state );

            #if DEBUG >= 2
            if ( wRank == 0 ) {
                fprintf(stdout, " Done solving the least squares problem.\n");
                fflush(stdout);
            }
            #endif

            // Extract the solution and add the seed back in
            F_array = F_alglib.getcontent();
            std::vector<double> Psi_vector(F_array,        F_array +     Npts),
                                Phi_vector(F_array + Npts, F_array + 2 * Npts);
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
            //alglib::sparsemv( LHS_matr, F_alglib, proj_vals );
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
    }

    vars_to_write.push_back("Psi");
    vars_to_write.push_back("Phi");

    initialize_output_file( source_data, vars_to_write, output_fname.c_str(), -1);

    if (not(constants::MINIMAL_OUTPUT)) {
        write_field_to_output(full_u_lon_tor,  "u_lon_tor",  starts, counts, output_fname.c_str(), &unmask);
        write_field_to_output(full_u_lat_tor,  "u_lat_tor",  starts, counts, output_fname.c_str(), &unmask);

        write_field_to_output(full_u_lon_pot,  "u_lon_pot",  starts, counts, output_fname.c_str(), &unmask);
        write_field_to_output(full_u_lat_pot,  "u_lat_pot",  starts, counts, output_fname.c_str(), &unmask);
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
