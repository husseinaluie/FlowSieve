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
        const std::vector<double> & mask,
        const std::vector<int>    & myCounts,
        const std::vector<int>    & myStarts,
        const MPI_Comm comm
        ) {

    int wRank, wSize;
    MPI_Comm_rank( comm, &wRank );
    MPI_Comm_size( comm, &wSize );

    const double rel_tol   = 1e-5;
    const int    max_iters = 1e5;

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

    int Ilat, Ilon, index, index_sub;

    // Fill in the land areas with zero velocity
    for (index = 0; index < (int)u_lon.size(); index++) {
        if (mask.at(index) == 0) {
            u_lon.at(index) = 0.;
            u_lat.at(index) = 0.;
        }
    }

    // Storage vectors
    std::vector<double> full_F(u_lon.size(), 0.),
        full_u_lon_tor(u_lon.size(), 0.),
        full_u_lat_tor(u_lon.size(), 0.),
        full_curl_term(u_lon.size(), 0.);

    // alglib variables
    alglib::real_1d_array rhs;
    std::vector<double> curl_term(Npts, 0.);
    rhs.attach_to_ptr(Npts, &curl_term[0]);

    alglib::linlsqrstate state;
    alglib::linlsqrreport report;

    alglib::real_1d_array F_alglib;

    double *F_array;

    //
    //// Build the LHS part of the problem (Lap)
    //
    fprintf(stdout, "Building the LHS of the least squares problem.\n");
    fflush(stdout);

    alglib::sparsematrix Lap;
    alglib::sparsecreate(Npts, Npts, Lap);

    toroidal_sparse_Lap(Lap, latitude, longitude, Nlat, Nlon, unmask);
    alglib::sparseconverttocrs(Lap);

    fprintf(stdout, "Declaring the least squares problem.\n");
    fflush(stdout);
    alglib::linlsqrcreate(Npts, Npts, state);
    alglib::linlsqrsetcond(state, rel_tol, rel_tol, max_iters);

    // Now do the solve!
    for (int Itime = 0; Itime < Ntime; ++Itime) {
        for (int Idepth = 0; Idepth < Ndepth; ++Idepth) {

            //
            //// Build the RHS of the problem (curl term)
            //
            fprintf(stdout, "Building the RHS of the least squares problem.\n");
            fflush(stdout);

            toroidal_curl_u_dot_er(curl_term, u_lon, u_lat, longitude, latitude, 
                    Itime, Idepth, Ntime, Ndepth, Nlat, Nlon, unmask);

            //
            //// Now apply the least-squares solver
            //
            fprintf(stdout, "Solving the least squares problem.\n");
            fflush(stdout);
            alglib::linlsqrsolvesparse(state, Lap, rhs);
            alglib::linlsqrresults(state, F_alglib, report);

            F_array = F_alglib.getcontent();
            std::vector<double> F_vector(F_array, F_array + Npts);

            // Get associated velocity
            std::vector<double> u_lon_tor(Npts, 0.), u_lat_tor(Npts, 0.);
            toroidal_vel_from_F(u_lon_tor, u_lat_tor, F_vector, longitude, latitude,
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

                    full_F.at(index) = F_vector.at(index_sub);

                    full_u_lon_tor.at(index) = u_lon_tor.at(index_sub);
                    full_u_lat_tor.at(index) = u_lat_tor.at(index_sub);

                    full_curl_term.at(index) = curl_term.at(index_sub);

                }
            }
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
    vars_to_write.push_back("F");
    vars_to_write.push_back("u_lon");
    vars_to_write.push_back("u_lat");
    vars_to_write.push_back("curl_term");

    char fname [50];
    snprintf(fname, 50, "toroidal_projection.nc");

    initialize_output_file(time, depth, longitude, latitude,
            mask, vars_to_write, fname, 0);
    
    write_field_to_output(full_u_lon_tor, "u_lon",     starts, counts, fname, &mask);
    write_field_to_output(full_u_lat_tor, "u_lat",     starts, counts, fname, &mask);
    write_field_to_output(full_F,         "F",         starts, counts, fname, &mask);
    write_field_to_output(full_curl_term, "curl_term", starts, counts, fname, &mask);

    add_attr_to_file("rel_tol", rel_tol, fname);
    add_attr_to_file("max_iters", (double) max_iters, fname);

}
