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


void Apply_Potential_Projection(
        const std::string output_fname,
        dataset & source_data,
        const std::vector<double> & seed,
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
                Ndepth  = myCounts.at(1),
                Nlat    = myCounts.at(2),
                Nlon    = myCounts.at(3);

    const int Npts = Nlat * Nlon;

    int Itime = 0, Idepth = 0, Ilat, Ilon; 
    size_t index, index_sub;

    // Get the velocity means ( will be stored in output file for reference )
    //      mean meridional flow will be captured by the potential (ulat * sin_lat is non-zero)
    std::vector<double> u_lon_means, u_lat_means;
    compute_spatial_average(u_lon_means, u_lon, dAreas, Ntime, Ndepth, Nlat, Nlon, mask);
    compute_spatial_average(u_lat_means, u_lat, dAreas, Ntime, Ndepth, Nlat, Nlon, mask);

    // Fill in the land areas with zero velocity
    #pragma omp parallel \
    default(none) \
    shared( u_lon, u_lat, mask ) \
    private( index )
    {
        #pragma omp for collapse(1) schedule(guided)
        for (index = 0; index < u_lon.size(); index++) {
            if (not(mask.at(index))) {
                u_lon.at(index) = 0.;
                u_lat.at(index) = 0.;
            }
        }
    }

    // Get the divergence of the original reference field, for comparison
    std::vector<double> full_div_orig(u_lon.size(), 0.);
    toroidal_vel_div(full_div_orig, u_lon, u_lat, longitude, latitude,
            Ntime, Ndepth, Nlat, Nlon, use_mask ? mask : unmask);

    // Storage vectors
    std::vector<double> 
        full_F(         u_lon.size(), 0. ),
        full_u_lon_pot( u_lon.size(), 0. ),
        full_u_lat_pot( u_lon.size(), 0. ),
        full_RHS(       u_lon.size(), 0. ),
        full_div_pot(   u_lon.size(), 0. ),
        full_seed(      u_lon.size(), 0. );

    // alglib variables
    alglib::real_1d_array rhs;
    std::vector<double> 
        div_term(   Npts, 0. ), 
        Lap_F_pot(  Npts, 0. ),
        F_seed(     Npts, 0. ),
        F_seed_Lap( Npts, 0. );

    // Copy the starting seed.
    if (single_seed) {
        #pragma omp parallel \
        default(none) \
        shared(F_seed, seed) \
        private( Ilat, Ilon, index ) \
        firstprivate( Nlon, Nlat )
        {
            #pragma omp for collapse(2) schedule(static)
            for (Ilat = 0; Ilat < Nlat; ++Ilat) {
                for (Ilon = 0; Ilon < Nlon; ++Ilon) {
                    index = Index(0, 0, Ilat, Ilon,
                                  1, 1, Nlat, Nlon);
                    F_seed.at(index) = seed.at(index);
                }
            }
        }
    }

    rhs.attach_to_ptr(Npts, &div_term[0]);

    alglib::linlsqrstate state;
    alglib::linlsqrreport report;

    alglib::real_1d_array F_alglib;

    double *F_array;

    //
    //// Build the LHS part of the problem (Lap)
    //
    if (wRank == 0) {
        fprintf(stdout, "Building the LHS of the least squares problem.\n");
        fflush(stdout);
    }

    alglib::sparsematrix Lap;
    alglib::sparsecreate(Npts, Npts, Lap);

    toroidal_sparse_Lap(Lap, source_data, Itime, Idepth, use_mask ? mask : unmask, weight_err);
    alglib::sparseconverttocrs(Lap);

    if (wRank == 0) {
        fprintf(stdout, "Declaring the least squares problem.\n");
        fflush(stdout);
    }
    alglib::linlsqrcreate(Npts, Npts, state);
    alglib::linlsqrsetcond(state, rel_tol, rel_tol, max_iters);

    // Now do the solve!
    for (int Itime = 0; Itime < Ntime; ++Itime) {
        for (int Idepth = 0; Idepth < Ndepth; ++Idepth) {

            // If single_seed == false, then we were provided seed values
            if (not(single_seed)) {
                #pragma omp parallel \
                default(none) \
                shared( F_seed, seed, Itime, Idepth ) \
                private( Ilat, Ilon, index, index_sub ) \
                firstprivate( Nlon, Nlat, Ndepth, Ntime )
                {
                    #pragma omp for collapse(2) schedule(static)
                    for (Ilat = 0; Ilat < Nlat; ++Ilat) {
                        for (Ilon = 0; Ilon < Nlon; ++Ilon) {
                            index = Index(Itime, Idepth, Ilat, Ilon,
                                          Ntime, Ndepth, Nlat, Nlon);
                            index_sub = Index(0, 0, Ilat, Ilon,
                                              1, 1, Nlat, Nlon);
                            F_seed.at(index_sub) = seed.at(index);
                        }
                    }
                }
            }

            //
            //// Build the RHS of the problem (curl term)
            //
            if ( (wRank == 0) and (Itime == 0) ) {
                fprintf(stdout, "Building the RHS of the least squares problem.\n");
                fflush(stdout);
            }

            // Get Lap of computed F term
            //
            // The 'seed' approach is a bit different than usual.
            //   There doesn't seem to be a way to provide an actual
            //   seed to the solver, so instead I'm going to modify
            //   the RHS of the problem. That is, if we're solving
            //   Ax = b, then we write x = x' + x0, where x0 is the
            //   seed / initial guess. We then instead solve the 
            //   problem Ax' = b - Ax0
            //
            toroidal_Lap_F(F_seed_Lap, F_seed, longitude, latitude,
                    Ntime, Ndepth, Nlat, Nlon, use_mask ? mask : unmask);

            #pragma omp parallel \
            default(none) \
            shared( div_term, full_div_orig, F_seed_Lap, dAreas, Itime, Idepth ) \
            private( Ilat, Ilon, index_sub, index ) \
            firstprivate( Nlon, Nlat, Ndepth, Ntime, weight_err )
            {
                #pragma omp for collapse(2) schedule(static)
                for (Ilat = 0; Ilat < Nlat; ++Ilat) {
                    for (Ilon = 0; Ilon < Nlon; ++Ilon) {
                        index = Index(Itime, Idepth, Ilat, Ilon,
                                      Ntime, Ndepth, Nlat, Nlon);
                        index_sub = Index(0, 0, Ilat, Ilon,
                                          1, 1, Nlat, Nlon);
                        div_term.at(index_sub) = full_div_orig.at(index) - F_seed_Lap.at(index_sub);

                        if (weight_err) { div_term.at(index_sub) *= dAreas.at(index_sub); }
                    }
                }
            }

            //
            //// Now apply the least-squares solver
            //
            if ( (wRank == 0) and (Itime == 0) ) {
                fprintf(stdout, "Solving the least squares problem.\n");
                fflush(stdout);
            }
            alglib::linlsqrsolvesparse(state, Lap, rhs);
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

            // Extract the solution and add the seed back in
            F_array = F_alglib.getcontent();
            std::vector<double> F_vector(F_array, F_array + Npts);
            for (size_t ii = 0; ii < F_vector.size(); ++ii) {
                F_vector.at(ii) += F_seed.at(ii);
            }

            // Get velocity associated to computed F field
            std::vector<double> 
                u_lon_pot(Npts, 0.), u_lat_pot(Npts, 0.), div_pot(Npts, 0.);
            potential_vel_from_F(u_lon_pot, u_lat_pot, F_vector, longitude, latitude,
                                Ntime, Ndepth, Nlat, Nlon, use_mask ? mask : unmask);
            toroidal_vel_div(div_pot, u_lon_pot, u_lat_pot, longitude, latitude,
                                Ntime, Ndepth, Nlat, Nlon, use_mask ? mask : unmask);

            //
            //// Store into the full arrays
            //
            #pragma omp parallel \
            default(none) \
            shared( full_u_lon_pot, u_lon_pot, full_u_lat_pot, u_lat_pot, \
                    full_div_pot, div_pot, full_F, F_vector, full_seed, F_seed, full_RHS, div_term, \
                    Itime, Idepth ) \
            private( Ilat, Ilon, index, index_sub ) \
            firstprivate( Nlon, Nlat, Ndepth, Ntime )
            {
                #pragma omp for collapse(2) schedule(static)
                for (Ilat = 0; Ilat < Nlat; ++Ilat) {
                    for (Ilon = 0; Ilon < Nlon; ++Ilon) {
                        index = Index(Itime, Idepth, Ilat, Ilon,
                                Ntime, Ndepth, Nlat, Nlon);

                        index_sub = Index(0, 0, Ilat, Ilon,
                                          1, 1, Nlat, Nlon);

                        // add the mean velocity back in
                        full_u_lon_pot.at(index) = u_lon_pot.at(index_sub);
                        full_u_lat_pot.at(index) = u_lat_pot.at(index_sub);

                        full_div_pot.at(  index) = div_pot.at(  index_sub);
                        full_F.at(        index) = F_vector.at( index_sub);

                        full_seed.at(index) = F_seed.at(   index_sub);
                        full_RHS.at( index) = div_term.at( index_sub);
                    }
                }
            }

            // Copy this solution into the seed for the next iteration
            if (single_seed) {
                F_seed = F_vector;
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
    vars_to_write.push_back("u_lon");
    vars_to_write.push_back("u_lat");

    vars_to_write.push_back("F");
    vars_to_write.push_back("F_seed");

    if (not(constants::MINIMAL_OUTPUT)) {
        vars_to_write.push_back("RHS");
    }

    initialize_output_file( source_data, vars_to_write, output_fname.c_str(), -1);
    
    write_field_to_output(full_u_lon_pot,  "u_lon",    starts, counts, output_fname.c_str(), &mask);
    write_field_to_output(full_u_lat_pot,  "u_lat",    starts, counts, output_fname.c_str(), &mask);

    write_field_to_output(full_F,          "F",        starts, counts, output_fname.c_str(), &unmask);
    write_field_to_output(full_seed,       "F_seed",   starts, counts, output_fname.c_str(), &unmask);

    if (not(constants::MINIMAL_OUTPUT)) {
        write_field_to_output(full_RHS,        "RHS",      starts, counts, output_fname.c_str(), &mask);
    }

    // Also write the mean velocity
    if (wRank == 0) {
        const char* dim_names[] = {"time", "depth"};
        const int ndims = 2;
        add_var_to_file("u_lon_mean", dim_names, ndims, output_fname.c_str());
        add_var_to_file("u_lat_mean", dim_names, ndims, output_fname.c_str());
    }
    MPI_Barrier(comm);
    size_t mean_starts[2] = { size_t(myStarts.at(0)),   size_t(myStarts.at(1)) };
    size_t mean_counts[2] = { size_t(Ntime),            size_t(Ndepth) };

    write_field_to_output( u_lon_means, "u_lon_mean", mean_starts, mean_counts, output_fname.c_str(), NULL );
    write_field_to_output( u_lat_means, "u_lat_mean", mean_starts, mean_counts, output_fname.c_str(), NULL );

    // Store some solver information
    add_attr_to_file("rel_tol",    rel_tol,                     output_fname.c_str());
    add_attr_to_file("max_iters",  (double) max_iters,          output_fname.c_str());
    add_attr_to_file("diff_order", (double) constants::DiffOrd, output_fname.c_str());
    add_attr_to_file("use_mask",   (double) use_mask,           output_fname.c_str());
    add_attr_to_file("weight_err", (double) weight_err,         output_fname.c_str());

}
