#include "../../constants.hpp"
#include "../../functions.hpp"
#include "../../netcdf_io.hpp"
#include "../../differentiation_tools.hpp"
#include <algorithm>
#include <vector>
#include <omp.h>
#include <math.h>


double get_dt_tor(
        const std::vector<double> & u_lon,
        const std::vector<double> & u_lat,
        const std::vector<double> & longitude,
        const std::vector<double> & latitude,
        const std::vector<double> & mask
        ) {

    double max_u = 0, max_v = 0;
    size_t index;
    const int Nlat = latitude.size();

    #pragma omp parallel \
    default(none) shared( u_lon, u_lat, mask ) \
    private(index) reduction(max : max_u, max_v)
    {
        #pragma omp for collapse(1) schedule(guided)
        for (index = 0; index < u_lon.size(); index++) {
            if (mask.at(index) == 1) {
                max_u = fmax( max_u, fabs(u_lon.at(index)) );
                max_v = fmax( max_v, fabs(u_lat.at(index)) );
            }
        }
    }

    const double dlon = longitude.at(1) - longitude.at(0);

    double min_cos_lat = 1.;
    for (int Ilat = 0; Ilat < Nlat; ++Ilat) {
        min_cos_lat = fmin( min_cos_lat, cos(latitude.at(Ilat)) );
    }

    double dlat = 1e5;
    for (int Ilat = 0; Ilat < Nlat-1; ++Ilat) {
        dlat = fmin( dlat, latitude.at(Ilat+1) - latitude.at(Ilat) );
    }
    
    double dt_1, dt_2;
    dt_1 = pow(  min_cos_lat * dlon * constants::R_earth, 2. );
    dt_2 = fmin( min_cos_lat * dlon, dlat ) * constants::R_earth / fmax(max_u, max_v);

    return 0.004 * fmin(dt_1, dt_2);
}


void apply_toroidal_projection(
        std::vector<double> & u_lon_tor,
        std::vector<double> & u_lat_tor,
        std::vector<double> & u_lon, // not const because we fill in
        std::vector<double> & u_lat, // masked areas with zeros
        const std::vector<double> & areas,
        const std::vector<double> & longitude,
        const std::vector<double> & latitude,
        const std::vector<double> & depth,
        const std::vector<double> & time,
        const std::vector<double> & mask,
        const std::vector<int>    & myCounts,
        const std::vector<int>    & myStarts,
        const MPI_Comm comm
        ) {

    int wRank, wSize;
    MPI_Comm_rank( comm, &wRank );
    MPI_Comm_size( comm, &wSize );

    // Create a 'no mask' mask variable
    //   we'll treat land values as zero velocity
    //   We do this because including land seems
    //   to introduce strong numerical issues
    std::vector<double> unmasked(mask.size(), 1.);

    int index;
    #pragma omp parallel \
    default(none) shared( mask, u_lon, u_lat ) private(index)
    {
        #pragma omp for collapse(1) schedule(guided)
        for (index = 0; index < (int)mask.size(); index++) {
            if (mask.at(index) == 0) { // Only land areas
                u_lon.at(index) = 0.;
                u_lat.at(index) = 0.;
            }
        }
    }

    const int Ntime   = myCounts.at(0);
    const int Ndepth  = myCounts.at(1);
    const int Nlat    = myCounts.at(2);
    const int Nlon    = myCounts.at(3);

    std::vector<double> curl_term(u_lon.size(), 0.), RHS(u_lon.size(), 0.), 
        work_arr_1(u_lon.size(), 0.), work_arr_2(u_lon.size(), 0.),
        Lap_F(u_lon.size(), 0.), F(u_lon.size(), 0.),
        RHS_old(u_lon.size(), 0.);

    double RHS_norm = 1.;

    //
    //// Compute curl term and its norm 
    //

    if (wRank == 0) {
        fprintf(stdout, "Computing curl term.\n");
        fflush(stdout);
    }
    toroidal_curl_u_dot_er(curl_term, u_lon, u_lat, longitude, latitude, 
            Ntime, Ndepth, Nlat, Nlon, unmasked);


    double curl_norm = 0., net_area = 0., err;
    int area_index, Itime, Idepth, Ilat, Ilon;

    #pragma omp parallel \
    default(none) \
    shared( latitude, longitude, F, unmasked, areas, curl_term ) \
    private(Itime, Idepth, Ilat, Ilon, index, area_index) \
    reduction(+ : curl_norm, net_area)
    {
        #pragma omp for collapse(1) schedule(guided)
        for (index = 0; index < (int)F.size(); index++) {

            Index1to4(index, Itime, Idepth, Ilat, Ilon,
                             Ntime, Ndepth, Nlat, Nlon);

            area_index = Index(0, 0, Ilat, Ilon,
                               1, 1, Nlat, Nlon);

            if (unmasked.at(index) == 1) { // Skip land areas
                curl_norm += pow(curl_term.at(index), 2) * areas.at(area_index);
                net_area += areas.at(area_index);
            }
        }
    }

    const double ref_norm = sqrt(curl_norm / net_area);
    if (wRank == 0) {
        fprintf(stdout, "Reference norm: %.3g\n", ref_norm);
        fflush(stdout);
    }

    //
    //// Get dt
    //
    const double dt = get_dt_tor(u_lon, u_lat, longitude, latitude, mask);

    if (wRank == 0) {
        fprintf(stdout, "Beginning iterations (using dt = %.3g).\n", dt);
        fflush(stdout);
    }

    //
    //// Initialize output file
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
    vars_to_write.push_back("Lap_F");
    vars_to_write.push_back("curl_term");
    vars_to_write.push_back("RHS");

    char fname [50];
    snprintf(fname, 50, "toroidal_projection.nc");

    initialize_output_file(time, depth, longitude, latitude,
            mask, vars_to_write, fname, 0);

    // Iterate until converged
    bool converged = false, write = false;

    const double conv_test = 1e-2; // 1e-3
    double check = 1.0;

    const int max_iters = 1e7; // 1e8
    const int iter_base = 5e4;
    int iter_check = 0;
    int iter = 0;

    while ( not(converged) ) {

        write = false;

        // Compute RHS
        RHS_norm = toroidal_comp_RHS( RHS, work_arr_1, work_arr_2, Lap_F,
                F, curl_term, longitude, latitude, 
                Ntime, Ndepth, Nlat, Nlon, unmasked, areas);

        // Update RHS_old (for first iteration)
        if (iter == 0) {
            #pragma omp parallel default(none) shared( unmasked, RHS, RHS_old ) private(index)
            {
                #pragma omp for collapse(1) schedule(guided)
                for (index = 0; index < (int)RHS.size(); ++index) {
                    if (unmasked.at(index) == 1) {
                        RHS_old.at(index) = RHS.at(index);
                    }
                }
            }
        }
    
        // Update F
        #pragma omp parallel default(none) shared( unmasked, F, RHS, RHS_old ) private(index)
        {
            #pragma omp for collapse(1) schedule(guided)
            for (index = 0; index < (int)F.size(); ++index) {
                if (unmasked.at(index) == 1) {
                    F.at(index) += dt * ( 1.5 * RHS.at(index) - 0.5 * RHS_old.at(index) );
                }
            }
        }

        // Update RHS_old
        #pragma omp parallel default(none) shared( unmasked, RHS, RHS_old ) private(index)
        {
            #pragma omp for collapse(1) schedule(guided)
            for (index = 0; index < (int)RHS.size(); ++index) {
                if (unmasked.at(index) == 1) {
                    RHS_old.at(index) = RHS.at(index);
                }
            }
        }

        // Test convergence
        err = RHS_norm / ref_norm;
        if ( err < check ) {
            if (wRank == 0) {
                fprintf(stdout, " .. err down to %.5e (%.3g / %.3g) after %d iters\n",
                        err, RHS_norm, ref_norm, iter);
                fflush(stdout);
            }
            write = true;
            check = check / pow(10, 1./4);
        }
        converged = err < conv_test;

        iter++;
        // Occasionally write out some information
        if (iter >= iter_check) {
            if (wRank == 0) {
                fprintf(stdout, " .. (1-err) down to %.5e (%.3g / %.3g) after %d iters\n",
                        1 - err, RHS_norm, ref_norm, iter);
                fflush(stdout);
            }
            write = true;
            iter_check += iter_base;
        }
        
        // If we've passed the iteration cap, then halt
        if (iter >= max_iters) {
            if (wRank == 0) {
                fprintf(stdout, 
                        " Exceeded iteration cap (%d iters),"
                        " halting with err %.5e (%.3g / %.3g).\n", 
                        iter, err, RHS_norm, ref_norm);
                fflush(stdout);
            }
            write = true;
            converged = true;
        }

        if (write) {
            // Extract the velocity
            if (wRank == 0) {
                fprintf(stdout, "Computing velocity from toroidal F.\n");
                fflush(stdout);
            }
            std::vector<double> vel_lon(F.size(), 1.), vel_lat(F.size(), 1.);
            toroidal_vel_from_F(vel_lon, vel_lat, F, longitude, latitude,
                    Ntime, Ndepth, Nlat, Nlon, unmasked);

            write_field_to_output(vel_lon,   "u_lon",     starts, counts, fname, &mask);
            write_field_to_output(vel_lat,   "u_lat",     starts, counts, fname, &mask);
            write_field_to_output(F,         "F",         starts, counts, fname, &mask);
            write_field_to_output(Lap_F,     "Lap_F",     starts, counts, fname, &mask);
            write_field_to_output(curl_term, "curl_term", starts, counts, fname, &mask);
            write_field_to_output(RHS,       "RHS",       starts, counts, fname, &mask);
        }
    }
} 
