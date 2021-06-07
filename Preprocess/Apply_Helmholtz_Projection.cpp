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
        const MPI_Comm comm
        ) {

    int wRank, wSize;
    MPI_Comm_rank( comm, &wRank );
    MPI_Comm_size( comm, &wSize );

    // Create some tidy names for variables
    const std::vector<double>   &time       = source_data.time,
                                &depth      = source_data.depth,
                                &latitude   = source_data.latitude,
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
                Ndepth  = myCounts.at(1),
                Nlat    = myCounts.at(2),
                Nlon    = myCounts.at(3);

    const size_t Npts = Nlat * Nlon;

    int Itime, Idepth, Ilat, Ilon;
    size_t index, index_sub;

    // Fill in the land areas with zero velocity
    #pragma omp parallel default(none) shared( u_lon, u_lat, mask ) private( index )
    {
        #pragma omp for collapse(1) schedule(guided)
        for (index = 0; index < u_lon.size(); index++) {
            if (not(mask.at(index))) {
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
        u_lon_tor_seed(  u_lon.size(), 0. ),
        u_lat_tor_seed(  u_lon.size(), 0. ),
        u_lon_pot_seed(  u_lon.size(), 0. ),
        u_lat_pot_seed(  u_lon.size(), 0. );

    // alglib variables
    alglib::real_1d_array rhs;
    std::vector<double> 
        RHS_vector( 2 * Npts, 0.),
        Psi_seed(       Npts, 0.),
        Phi_seed(       Npts, 0.);
    

    // Copy the starting seed.
    if (single_seed) {
        #pragma omp parallel default(none) shared(Psi_seed, Phi_seed, seed_tor, seed_pot) private( Ilat, Ilon, index )
        {
            #pragma omp for collapse(2) schedule(static)
            for (Ilat = 0; Ilat < Nlat; ++Ilat) {
                for (Ilon = 0; Ilon < Nlon; ++Ilon) {
                    index = Index(0, 0, Ilat, Ilon, 1, 1, Nlat, Nlon);
                    Psi_seed.at(index) = seed_tor.at(index);
                    Phi_seed.at(index) = seed_pot.at(index);
                }
            }
        }
    }

    rhs.attach_to_ptr( 2 * Npts, &RHS_vector[0] );

    alglib::linlsqrstate state;
    alglib::linlsqrreport report;

    alglib::real_1d_array F_alglib;

    double *F_array;

    //
    //// Build the LHS part of the problem
    //      Ordering is: [  u_from_psi      u_from_phi   ] *  [ psi ]   =    [  u   ]
    //                   [  v_from_psi      v_from_phi   ]    [ phi ]        [  v   ]
    //
    #if DEBUG >= 1
    if (wRank == 0) {
        fprintf(stdout, "Building the LHS of the least squares problem.\n");
        fflush(stdout);
    }
    #endif

    alglib::sparsematrix LHS_matr;
    alglib::sparsecreate(2*Npts, 2*Npts, LHS_matr);

    // Put in {u,v}_from_{psi,phi} bits
    sparse_vel_from_PsiPhi( LHS_matr, source_data, Itime, Idepth, use_mask ? mask : unmask, weight_err );

    alglib::sparseconverttocrs(LHS_matr);

    #if DEBUG >= 1
    if (wRank == 0) {
        fprintf(stdout, "Declaring the least squares problem.\n");
        fflush(stdout);
    }
    #endif
    alglib::linlsqrcreate(2*Npts, 2*Npts, state);
    alglib::linlsqrsetcond(state, rel_tol, rel_tol, max_iters);

    // Now do the solve!
    for (int Itime = 0; Itime < Ntime; ++Itime) {
        for (int Idepth = 0; Idepth < Ndepth; ++Idepth) {

            if (not(single_seed)) {
                // If single_seed == false, then we were provided seed values, pull out the appropriate values here
                #pragma omp parallel \
                default(none) \
                shared( Psi_seed, Phi_seed, seed_tor, seed_pot, Itime, Idepth ) \
                private( Ilat, Ilon, index, index_sub )
                {
                    #pragma omp for collapse(2) schedule(static)
                    for (Ilat = 0; Ilat < Nlat; ++Ilat) {
                        for (Ilon = 0; Ilon < Nlon; ++Ilon) {
                            index = Index(Itime, Idepth, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon);
                            index_sub = Index(0, 0, Ilat, Ilon, 1, 1, Nlat, Nlon);
                            Psi_seed.at(index_sub) = seed_tor.at(index);
                            Phi_seed.at(index_sub) = seed_pot.at(index);
                        }
                    }
                }
            }

            // Get velocity from seed
            toroidal_vel_from_F(  u_lon_tor_seed, u_lat_tor_seed, Psi_seed, longitude, latitude, Ntime, Ndepth, Nlat, Nlon, use_mask ? mask : unmask);
            potential_vel_from_F( u_lon_pot_seed, u_lat_pot_seed, Phi_seed, longitude, latitude, Ntime, Ndepth, Nlat, Nlon, use_mask ? mask : unmask);

            #if DEBUG >= 1
            if ( (wRank == 0) and (Itime == 0) ) {
                fprintf(stdout, "Building the RHS of the least squares problem.\n");
                fflush(stdout);
            }
            #endif

            //
            //// Set up the RHS_vector
            //
            //  Subtract off the velocity from the seed
            //
            
            #pragma omp parallel default(none) \
            shared( dAreas, Itime, Idepth, RHS_vector, u_lon, u_lon_tor_seed, u_lon_pot_seed, u_lat, u_lat_tor_seed, u_lat_pot_seed ) \
            private( Ilat, Ilon, index, index_sub )
            {
                #pragma omp for collapse(2) schedule(static)
                for (Ilat = 0; Ilat < Nlat; ++Ilat) {
                    for (Ilon = 0; Ilon < Nlon; ++Ilon) {
                        index_sub = Index( 0,     0,      Ilat, Ilon, 1,     1,      Nlat, Nlon);
                        index     = Index( Itime, Idepth, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon);

                        RHS_vector.at(       index_sub) = u_lon.at(index) - u_lon_tor_seed.at(index_sub) - u_lon_pot_seed.at(index_sub);
                        RHS_vector.at(Npts + index_sub) = u_lat.at(index) - u_lat_tor_seed.at(index_sub) - u_lat_pot_seed.at(index_sub);

                        if ( weight_err ) {
                            RHS_vector.at(       index_sub) *= dAreas.at(index_sub);
                            RHS_vector.at(Npts + index_sub) *= dAreas.at(index_sub);
                        }
                    }
                }
            }

            //
            //// Now apply the least-squares solver
            //
            #if DEBUG >= 1
            if ( (wRank == 0) and (Itime == 0) ) {
                fprintf(stdout, "Solving the least squares problem.\n");
                fflush(stdout);
            }
            #endif
            alglib::linlsqrsolvesparse(state, LHS_matr, rhs);
            alglib::linlsqrresults(state, F_alglib, report);

            #if DEBUG >= 1
            if      (report.terminationtype == 1) { fprintf(stdout, "Termination type: absolulte tolerance reached.\n"); }
            else if (report.terminationtype == 4) { fprintf(stdout, "Termination type: relative tolerance reached.\n"); }
            else if (report.terminationtype == 5) { fprintf(stdout, "Termination type: maximum number of iterations reached.\n"); }
            else if (report.terminationtype == 7) { fprintf(stdout, "Termination type: round-off errors prevent further progress.\n"); }
            else if (report.terminationtype == 8) { fprintf(stdout, "Termination type: user requested (?)\n"); }
            else                                  { fprintf(stdout, "Termination type: unknown\n"); }
            #endif

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
            if ( (wRank == 0) and (Itime == 0) ) {
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
            #if DEBUG >= 1
            if ( (wRank == 0) and (Itime == 0) ) {
                fprintf(stdout, " Extracting velocities and divergence from toroidal field.\n");
                fflush(stdout);
            }
            #endif

            std::vector<double> u_lon_tor(Npts, 0.), u_lat_tor(Npts, 0.), u_lon_pot(Npts, 0.), u_lat_pot(Npts, 0.);
            toroidal_vel_from_F(  u_lon_tor, u_lat_tor, Psi_vector, longitude, latitude, Ntime, Ndepth, Nlat, Nlon, use_mask ? mask : unmask);
            potential_vel_from_F( u_lon_pot, u_lat_pot, Phi_vector, longitude, latitude, Ntime, Ndepth, Nlat, Nlon, use_mask ? mask : unmask);

            //
            //// Store into the full arrays
            //
            #if DEBUG >= 1
            if ( (wRank == 0) and (Itime == 0) ) {
                fprintf(stdout, " Storing values into output arrays\n");
                fflush(stdout);
            }
            #endif
            #pragma omp parallel \
            default(none) \
            shared( full_u_lon_tor, u_lon_tor, full_u_lat_tor, u_lat_tor, \
                    full_u_lon_pot, u_lon_pot, full_u_lat_pot, u_lat_pot, \
                    full_Psi, full_Phi, Psi_vector, Phi_vector, \
                    Itime, Idepth ) \
            private( Ilat, Ilon, index, index_sub )
            {
                #pragma omp for collapse(2) schedule(static)
                for (Ilat = 0; Ilat < Nlat; ++Ilat) {
                    for (Ilon = 0; Ilon < Nlon; ++Ilon) {
                        index = Index(Itime, Idepth, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon);

                        index_sub = Index(0, 0, Ilat, Ilon, 1, 1, Nlat, Nlon);

                        // add the mean velocity back in
                        full_u_lon_tor.at(index) = u_lon_tor.at(index_sub) ;
                        full_u_lat_tor.at(index) = u_lat_tor.at(index_sub) ;

                        full_u_lon_pot.at(index) = u_lon_pot.at(index_sub) ;
                        full_u_lat_pot.at(index) = u_lat_pot.at(index_sub) ;

                        full_Psi.at(index) = Psi_vector.at( index_sub );
                        full_Phi.at(index) = Phi_vector.at( index_sub );
                    }
                }
            }

            // If we don't have a seed for the next iteration, use this solution as the seed
            if (single_seed) {
                Psi_seed = Psi_vector;
                Phi_seed = Phi_vector;
            }

            #if DEBUG >= 1
            fprintf(stdout, "  --  --  Rank %d done depth %d\n", wRank, Idepth + myStarts.at(1) );
            fflush(stdout);
            #endif

        }

        #if DEBUG >= 1
        fprintf(stdout, " -- Rank %d done time %d\n", wRank, Itime + myStarts.at(0) );
        fflush(stdout);
        #endif
    }

    //
    //// Write the output
    //

    const int ndims = 4;
    size_t starts[ndims] = {
        size_t(myStarts.at(0)), size_t(myStarts.at(1)), 
        size_t(myStarts.at(2)), size_t(myStarts.at(3))
    };
    size_t counts[ndims] = {
        size_t(Ntime), size_t(Ndepth), 
        size_t(Nlat),  size_t(Nlon)
    };

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
    add_attr_to_file("rel_tol",    rel_tol,                     output_fname.c_str());
    add_attr_to_file("max_iters",  (double) max_iters,          output_fname.c_str());
    add_attr_to_file("diff_order", (double) constants::DiffOrd, output_fname.c_str());
    add_attr_to_file("use_mask",   (double) use_mask,           output_fname.c_str());
    add_attr_to_file("weight_err", (double) weight_err,         output_fname.c_str());


    //
    //// At the very end, compute the L2 error for each time/depth
    //

    std::vector<double> projection_error(   Ntime * Ndepth, 0. ),
                        projection_KE(      Ntime * Ndepth, 0. ),
                        toroidal_KE(        Ntime * Ndepth, 0. ),
                        potential_KE(       Ntime * Ndepth, 0. ),
                        original_KE(        Ntime * Ndepth, 0. );
    double total_area, error, tor_KE, pot_KE, proj_KE, orig_KE;
    for (int Itime = 0; Itime < Ntime; ++Itime) {
        for (int Idepth = 0; Idepth < Ndepth; ++Idepth) {
            #pragma omp parallel \
            default(none) \
            shared( full_u_lon_tor, full_u_lat_tor, full_u_lon_pot, full_u_lat_pot, \
                    u_lon, u_lat, Itime, Idepth, dAreas ) \
            reduction(+ : total_area, error, tor_KE, pot_KE, proj_KE, orig_KE) \
            private( Ilat, Ilon, index, index_sub )
            {
                #pragma omp for collapse(2) schedule(static)
                for (Ilat = 0; Ilat < Nlat; ++Ilat) {
                    for (Ilon = 0; Ilon < Nlon; ++Ilon) {
                        index_sub = Index( 0,     0,      Ilat, Ilon, 1,     1,      Nlat, Nlon);
                        index     = Index( Itime, Idepth, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon);

                        total_area += dAreas.at(index_sub);

                        error += dAreas.at(index_sub) * (
                                        pow( u_lon.at(index) - full_u_lon_tor.at(index) - full_u_lon_pot.at(index) , 2.)
                                     +  pow( u_lat.at(index) - full_u_lat_tor.at(index) - full_u_lat_pot.at(index) , 2.)
                                );

                        tor_KE += dAreas.at(index_sub) * ( pow( full_u_lon_tor.at(index), 2.) + pow( full_u_lat_tor.at(index), 2.) );
                        pot_KE += dAreas.at(index_sub) * ( pow( full_u_lon_pot.at(index), 2.) + pow( full_u_lat_pot.at(index), 2.) );

                        proj_KE += dAreas.at(index_sub) * (
                                        pow( full_u_lon_tor.at(index) + full_u_lon_pot.at(index) , 2.)
                                     +  pow( full_u_lat_tor.at(index) + full_u_lat_pot.at(index) , 2.)
                                );

                        orig_KE += dAreas.at(index_sub) * ( pow( u_lon.at(index), 2.) + pow( u_lat.at(index), 2.) );
                    }
                }
            }
            size_t int_index = Index( Itime, Idepth, 0, 0, Ntime, Ndepth, 1, 1);
            projection_error.at( int_index ) = error    / total_area;
            projection_KE.at(    int_index ) = proj_KE  / total_area;
            toroidal_KE.at(      int_index ) = tor_KE   / total_area;
            potential_KE.at(     int_index ) = pot_KE   / total_area;
            original_KE.at(      int_index ) = orig_KE  / total_area;
        }
    }

    const char* dim_names[] = {"time", "depth"};
    const int ndims_error = 2;
    add_var_to_file( "projection_error",    dim_names, ndims_error, output_fname.c_str() );
    add_var_to_file( "projection_KE",       dim_names, ndims_error, output_fname.c_str() );
    add_var_to_file( "toroidal_KE",         dim_names, ndims_error, output_fname.c_str() );
    add_var_to_file( "potential_KE",        dim_names, ndims_error, output_fname.c_str() );
    add_var_to_file( "original_KE",         dim_names, ndims_error, output_fname.c_str() );

    size_t starts_error[ndims_error] = { size_t(myStarts.at(0)), size_t(myStarts.at(1)) };
    size_t counts_error[ndims_error] = { size_t(Ntime), size_t(Ndepth) };

    write_field_to_output( projection_error, "projection_error", starts_error, counts_error, output_fname.c_str() );
    write_field_to_output( projection_KE,    "projection_KE",    starts_error, counts_error, output_fname.c_str() );
    write_field_to_output( toroidal_KE,      "toroidal_KE",      starts_error, counts_error, output_fname.c_str() );
    write_field_to_output( potential_KE,     "potential_KE",     starts_error, counts_error, output_fname.c_str() );
    write_field_to_output( original_KE,      "original_KE",      starts_error, counts_error, output_fname.c_str() );

}
