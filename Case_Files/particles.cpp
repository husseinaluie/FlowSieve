#include <fenv.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <algorithm>
#include <math.h>
#include <vector>
#include <mpi.h>
#include <omp.h>
#include <cassert>

#include "../netcdf_io.hpp"
#include "../functions.hpp"
#include "../constants.hpp"
#include "../particles.hpp"

int main(int argc, char *argv[]) {
    
    // PERIODIC_Y implies UNIFORM_LAT_GRID
    static_assert( (constants::UNIFORM_LAT_GRID) or (not(constants::PERIODIC_Y)),
            "PERIODIC_Y requires UNIFORM_LAT_GRID.\n"
            "Please update constants.hpp accordingly.\n");

    static_assert( ( (constants::PERIODIC_X) and (not(constants::PERIODIC_Y)) ),
            "The particles routine currently requires globe-like periodicity.\n"
            "Please update constants.hpp accordingly.\n");


    // Specify the number of OpenMP threads
    //   and initialize the MPI world
    int thread_safety_provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &thread_safety_provided);
    //MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI::ERRORS_THROW_EXCEPTIONS);

    int wRank=-1, wSize=-1;
    MPI_Comm_rank( MPI_COMM_WORLD, &wRank );
    MPI_Comm_size( MPI_COMM_WORLD, &wSize );

    // Default input file name
    const int str_len = 50;
    char    input_fname[   str_len], out_name[      str_len],
            flag_version[  str_len], flag_input[    str_len],
            zonal_vel_name[str_len], flag_zonal_vel[str_len], 
            merid_vel_name[str_len], flag_merid_vel[str_len];

    snprintf(zonal_vel_name,    str_len, "uo");
    snprintf(merid_vel_name,    str_len, "vo");

    snprintf(input_fname,       str_len, "input.nc");
    snprintf(flag_version,      str_len, "--version");
    snprintf(flag_input,        str_len, "--input_file");
    snprintf(flag_zonal_vel,    str_len, "--zonal_vel");
    snprintf(flag_merid_vel,    str_len, "--merid_vel");

    // Parse command-line flags
    int ii = 1;
    while(ii < argc) {  
        if (wRank == 0) {
            fprintf(stdout, "Argument %d : %s\n", ii, argv[ii]);
        }

        if (strcmp(argv[ii], flag_version) == 0) {
            // check if the flag is 'version'
            if (wRank == 0) {
                //print_compile_info(filter_scales);
            }
            return 0;
        } else if (strcmp(argv[ii], flag_input) == 0) {
            // check if the flag is 'input_file' and, if it is
            //   the next input is then used as the filename
            //   of the input
            snprintf(input_fname, 50, argv[ii+1]);
            ++ii;
        } else if (strcmp(argv[ii], flag_zonal_vel) == 0) {
            // check if we're given the name of the zonal velocity var
            snprintf(zonal_vel_name, str_len, argv[ii+1]);
            ++ii;
        } else if (strcmp(argv[ii], flag_merid_vel) == 0) {
            // check if we're given the name of the meridional velocity var
            snprintf(merid_vel_name, str_len, argv[ii+1]);
            ++ii;
        } else {
            // Otherwise, the flag is unrecognized
            if (wRank == 0) {
                fprintf(stderr, "Flag %s not recognized.\n", argv[ii]);
            }
            return -1;
        }
        ++ii;
    }

    #if DEBUG >= 0
    if (wRank == 0) {
        fprintf(stdout, "\n\n");
        fprintf(stdout, "Compiled at %s on %s.\n", __TIME__, __DATE__);
        fprintf(stdout, "  Version %d.%d.%d \n", 
                MAJOR_VERSION, MINOR_VERSION, PATCH_VERSION);
        fprintf(stdout, "\n");
        fflush(stdout);
    }
    #endif

    #if DEBUG >= 0
    if (wRank == 0) {
        if (constants::CARTESIAN) { 
            fprintf(stdout, "Using Cartesian coordinates.\n");
        } else {
            fprintf(stdout, "Using spherical coordinates.\n");
        }
        fflush(stdout);
    }
    #endif
    MPI_Barrier(MPI_COMM_WORLD);

    // Print processor assignments
    const int max_threads = omp_get_max_threads();
    omp_set_num_threads( max_threads );
    #if DEBUG >= 2
    int tid, nthreads;
    #pragma omp parallel default(none) private(tid, nthreads) \
        shared(stdout) firstprivate(wRank, wSize)
    {
        tid = omp_get_thread_num();
        nthreads = omp_get_num_threads();
        fprintf(stdout, "Hello from thread %d of %d on processor %d of %d.\n", 
                tid+1, nthreads, wRank+1, wSize);
    }
    fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD);
    #endif

    std::vector<double> longitude, latitude, time, depth;
    std::vector<double> u_lon, u_lat;
    std::vector<double> mask;
    std::vector<int> myCounts, myStarts;
    size_t II;

    // Read in source data / get size information
    #if DEBUG >= 1
    if (wRank == 0) { fprintf(stdout, "Reading in source data.\n\n"); }
    #endif

    // Read in the grid coordinates
    read_var_from_file(longitude, "longitude", input_fname);
    read_var_from_file(latitude,  "latitude",  input_fname);
    read_var_from_file(time,      "time",      input_fname);
    read_var_from_file(depth,     "depth",     input_fname);

    // Read in the velocity fields
    read_var_from_file(u_lon, zonal_vel_name, input_fname, &mask, &myCounts, &myStarts);
    read_var_from_file(u_lat, merid_vel_name, input_fname, &mask, &myCounts, &myStarts);

    const int Ntime  = time.size();
    #if DEBUG >= 1
    const int Ndepth = depth.size();
    #endif
    const int Nlon   = longitude.size();
    const int Nlat   = latitude.size();

    #if DEBUG >= 1
    fprintf(stdout, "Processor %d has (%d, %d, %d, %d) from (%d, %d, %d, %d)\n", 
            wRank, 
            myCounts[0], myCounts[1], myCounts[2], myCounts[3],
            Ntime, Ndepth, Nlat, Nlon);
    fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD);
    #endif

    if (not(constants::CARTESIAN)) {
        // Convert coordinate to radians
        if (wRank == 0) { fprintf(stdout, "Converting to radians.\n\n"); }
        int ii;
        const double D2R = M_PI / 180.;
        #pragma omp parallel default(none) private(ii) shared(longitude, latitude)
        { 
            #pragma omp for collapse(1) schedule(static)
            for (ii = 0; ii < Nlon; ii++) { longitude.at(ii) = longitude.at(ii) * D2R; }

            #pragma omp for collapse(1) schedule(static)
            for (ii = 0; ii < Nlat; ii++) { latitude.at(ii) = latitude.at(ii) * D2R; }
        }
    }

    // Set the output times
    const int Nouts = Ntime * 4;
    const double start_time = time.front(),
                 final_time = time.at(time.size()-2);

    std::vector<double> target_times(Nouts);
    for ( II = 0; II < target_times.size(); ++II ) {
        target_times.at(II) = start_time 
            + (double)II * ( final_time - start_time ) / Nouts;
    }

    // Get particle positions
    const int Npts = 20;
    std::vector<double> starting_lat(Npts), starting_lon(Npts);
    particles_initial_positions(starting_lat, starting_lon, Npts,
            latitude, longitude, mask);

    // Trajectories dimension (essentially just a numbering)
    std::vector<double> trajectories(Npts);
    for (II = 0; II < Npts; ++II) { trajectories.at(II) = (double) II; };

    // Initialize particle output file
    snprintf(out_name, str_len, "particles.nc");
    std::vector<std::string> vars_to_write;
    initialize_particle_file(target_times, trajectories, vars_to_write, out_name);

    std::vector<double> part_lon_hist(Npts * Nouts, constants::fill_value), 
                        part_lat_hist(Npts * Nouts, constants::fill_value);

    for (II = 0; II < Npts; ++II) {
        part_lon_hist.at(II) = starting_lon.at(II);
        part_lat_hist.at(II) = starting_lat.at(II);
    }

    size_t starts[2], counts[2];
    starts[0] = 0;
    counts[0] = Nouts;

    starts[1] = 0;
    counts[1] = Npts;

    write_field_to_output(part_lon_hist, "longitude", starts, counts, out_name);
    write_field_to_output(part_lat_hist, "latitude",  starts, counts, out_name);

    // Convert to seconds
    for (II = 0; II < time.size();         ++II) { time.at(II)         *= 60 * 60; }
    for (II = 0; II < target_times.size(); ++II) { target_times.at(II) *= 60 * 60; }

    // Now do the particle routine
    particles_evolve_trajectories(
        part_lon_hist, part_lat_hist,
        starting_lat,  starting_lon,
        target_times,
        u_lon, u_lat,
        time, latitude, longitude,        
        mask);

    fprintf(stdout, "\n\nDone evolving particles, now writing outputs.\n\n");

    std::vector<double> out_mask(part_lon_hist.size());
    for (II = 0; II < out_mask.size(); ++II) {
        out_mask.at(II) = part_lon_hist.at(II) == constants::fill_value ? 0 : 1;
    }

    write_field_to_output(part_lon_hist, "longitude", 
            starts, counts, out_name, &out_mask);
    write_field_to_output(part_lat_hist, "latitude",  
            starts, counts, out_name, &out_mask);

    MPI_Finalize();
    return 0;
}
