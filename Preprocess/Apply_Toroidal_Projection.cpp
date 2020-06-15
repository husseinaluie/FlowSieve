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


void Apply_Toroidal_Projection(
        std::vector<double> & u_lon,
        std::vector<double> & u_lat,
        const std::vector<double> & time,
        const std::vector<double> & depth,
        const std::vector<double> & latitude,
        const std::vector<double> & longitude,
        const std::vector<double> & dAreas,
        const std::vector<double> & mask,
        const std::vector<int>    & myCounts,
        const std::vector<int>    & myStarts,
        const std::vector<double> & seed,
        const bool single_seed,
        const MPI_Comm comm
        ) {

    int wRank, wSize;
    MPI_Comm_rank( comm, &wRank );
    MPI_Comm_size( comm, &wSize );

    const double rel_tol    = 1e-7;
    const int    max_iters  = 1e2;
    const bool   weight_err = true; // weight err by cell size

    // Create a 'no mask' mask variable
    //   we'll treat land values as zero velocity
    //   We do this because including land seems
    //   to introduce strong numerical issues
    std::vector<double> unmask(mask.size(), 1.);

    const int Ntime   = myCounts.at(0);
    const int Ndepth  = myCounts.at(1);
    const int Nlat    = myCounts.at(2);
    const int Nlon    = myCounts.at(3);

    const int Npts = Nlat * Nlon;

    int Itime, Idepth, Ilat, Ilon, index, index_sub, mean_ind;

    // Get the velocity means
    std::vector<double> u_lon_means, u_lat_means;
    compute_spatial_average(u_lon_means, u_lon, dAreas, Ntime, Ndepth, Nlat, Nlon, mask);
    compute_spatial_average(u_lat_means, u_lat, dAreas, Ntime, Ndepth, Nlat, Nlon, mask);

    // Fill in the land areas with zero velocity
    //   also subtract the mean off (will be added back later)
    for (index = 0; index < (int)u_lon.size(); index++) {
        if (mask.at(index) == 0) {
            u_lon.at(index) = 0.;
            u_lat.at(index) = 0.;
        } else {
            Index1to4(index, Itime, Idepth, Ilat, Ilon,
                             Ntime, Ndepth, Nlat, Nlon);
            mean_ind  = Index(0, 0, Itime, Idepth,
                              1, 1, Ntime, Ndepth);

            u_lon.at(index) -= u_lon_means.at(mean_ind);
            u_lat.at(index) -= u_lat_means.at(mean_ind);
        }
    }

    // Get the divergence of the original reference field, for comparison
    std::vector<double> full_div_orig(u_lon.size(), 0.);
    toroidal_vel_div(full_div_orig, u_lon, u_lat, longitude, latitude,
            Ntime, Ndepth, Nlat, Nlon, unmask);

    // Storage vectors
    std::vector<double> 
        full_F(         u_lon.size(), 0. ),
        full_u_lon_tor( u_lon.size(), 0. ),
        full_u_lat_tor( u_lon.size(), 0. ),
        full_div_tor(   u_lon.size(), 0. ),
        full_RHS(       u_lon.size(), 0. ),
        full_seed(      u_lon.size(), 0. );

    // alglib variables
    alglib::real_1d_array rhs;
    std::vector<double> 
        curl_term( Npts, 0.), 
        Lap_F_tor( Npts, 0.),
        F_seed(    Npts, 0.),
        F_seed_Lap(Npts, 0.);
    

    // Copy the starting seed.
    if (single_seed) {
        for (Ilat = 0; Ilat < Nlat; ++Ilat) {
            for (Ilon = 0; Ilon < Nlon; ++Ilon) {
                index = Index(0, 0, Ilat, Ilon,
                              1, 1, Nlat, Nlon);
                F_seed.at(index) = seed.at(index);
            }
        }
    }

    rhs.attach_to_ptr(Npts, &curl_term[0]);

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

    toroidal_sparse_Lap(Lap, latitude, longitude, Itime, Idepth,
            Ntime, Ndepth, Nlat, Nlon, mask, dAreas, weight_err);
            //Ntime, Ndepth, Nlat, Nlon, unmask, dAreas, weight_err);
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
                    Ntime, Ndepth, Nlat, Nlon, unmask);

            toroidal_curl_u_dot_er(curl_term, u_lon, u_lat, longitude, latitude, 
                    Itime, Idepth, Ntime, Ndepth, Nlat, Nlon, unmask, &F_seed_Lap);

            if (weight_err) {
                // Weight by area if requested
                for (Ilat = 0; Ilat < Nlat; ++Ilat) {
                    for (Ilon = 0; Ilon < Nlon; ++Ilon) {
                        index_sub = Index(0, 0, Ilat, Ilon,
                                          1, 1, Nlat, Nlon);

                        curl_term.at(index_sub) *= dAreas.at(index_sub);
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

            // Extract the solution and add the seed back in
            F_array = F_alglib.getcontent();
            std::vector<double> F_vector(F_array, F_array + Npts);
            for (size_t ii = 0; ii < F_vector.size(); ++ii) {
                F_vector.at(ii) += F_seed.at(ii);
            }

            // Get velocity associated to computed F field
            std::vector<double> 
                u_lon_tor(Npts, 0.), u_lat_tor(Npts, 0.), div_tor(Npts, 0.);
            toroidal_vel_from_F(u_lon_tor, u_lat_tor, F_vector, longitude, latitude,
                                Ntime, Ndepth, Nlat, Nlon, unmask);
            toroidal_vel_div(div_tor, u_lon_tor, u_lat_tor, longitude, latitude,
                                Ntime, Ndepth, Nlat, Nlon, unmask);

            // Get velocity associated with seed F field
            std::vector<double> 
                u_lon_seed(Npts, 0.), u_lat_seed(Npts, 0.);
            toroidal_vel_from_F(u_lon_seed, u_lat_seed, F_seed, longitude, latitude,
                                Ntime, Ndepth, Nlat, Nlon, unmask);

            // Get Lap of computed F term
            toroidal_Lap_F(Lap_F_tor, F_vector, longitude, latitude,
                    Ntime, Ndepth, Nlat, Nlon, unmask);

            //
            //// Store into the full arrays
            //
            for (Ilat = 0; Ilat < Nlat; ++Ilat) {
                for (Ilon = 0; Ilon < Nlon; ++Ilon) {
                    index = Index(Itime, Idepth, Ilat, Ilon,
                                  Ntime, Ndepth, Nlat, Nlon);

                    index_sub = Index(0, 0, Ilat, Ilon,
                                      1, 1, Nlat, Nlon);
                    mean_ind  = Index(0, 0, Itime, Idepth,
                                      1, 1, Ntime, Ndepth);

                    // add the mean velocity back in
                    full_u_lon_tor.at(index) = u_lon_tor.at(index_sub) 
                        + u_lon_means.at(mean_ind);
                    full_u_lat_tor.at(index) = u_lat_tor.at(index_sub)
                        + u_lat_means.at(mean_ind);

                    full_div_tor.at(  index) = div_tor.at(  index_sub);
                    full_F.at(        index) = F_vector.at( index_sub);

                    full_seed.at(index) = F_seed.at(   index_sub);
                    full_RHS.at( index) = curl_term.at(index_sub);

                }
            }

            // Copy this solution into the seed for the next iteration
            if (single_seed) {
                F_seed = F_vector;
            }

            // 
            #if DEBUG >= 0
            fprintf(stdout, "    Rank %d done time %d\n", wRank, Itime + myStarts.at(0) );
            fflush(stdout);
            #endif

        }
    }

    //
    //// Write the output
    //

    const int ndims = 4;
    size_t starts[ndims] = {
        size_t(myStarts.at(0)), size_t(myStarts.at(1)), 
        size_t(myStarts.at(2)), size_t(myStarts.at(3))};
    size_t counts[ndims] = {
        size_t(Ntime), size_t(Ndepth), 
        size_t(Nlat), size_t(Nlon)};

    std::vector<std::string> vars_to_write;
    vars_to_write.push_back("u_lon");
    vars_to_write.push_back("u_lat");

    vars_to_write.push_back("F");
    vars_to_write.push_back("F_seed");

    vars_to_write.push_back("RHS");

    vars_to_write.push_back("div_tor");
    vars_to_write.push_back("div_orig");

    char fname [50];
    snprintf(fname, 50, "toroidal_projection.nc");

    initialize_output_file(time, depth, longitude, latitude,
            mask, vars_to_write, fname, 0);
    
    write_field_to_output(full_u_lon_tor,  "u_lon",    starts, counts, fname, &mask);
    write_field_to_output(full_u_lat_tor,  "u_lat",    starts, counts, fname, &mask);

    write_field_to_output(full_RHS,        "RHS",      starts, counts, fname, &mask);

    write_field_to_output(full_F,          "F",        starts, counts, fname, &unmask);
    write_field_to_output(full_seed,       "F_seed",   starts, counts, fname, &mask);

    write_field_to_output(full_div_tor,    "div_tor",  starts, counts, fname, &mask);
    write_field_to_output(full_div_orig,   "div_orig", starts, counts, fname, &mask);

    add_attr_to_file("rel_tol",    rel_tol,                     fname);
    add_attr_to_file("max_iters",  (double) max_iters,          fname);
    add_attr_to_file("diff_order", (double) constants::DiffOrd, fname);

}
