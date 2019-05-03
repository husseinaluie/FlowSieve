#include <fenv.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <algorithm>
#include <math.h>
#include <vector>
#include <mpi.h>
#include <omp.h>

#include "netcdf_io.hpp"
#include "functions.hpp"
#include "constants.hpp"

int main(int argc, char *argv[]) {

    // Enable all floating point exceptions but FE_INEXACT
    //feenableexcept(FE_ALL_EXCEPT & ~FE_INEXACT);

    // Specify the number of OpenMP threads
    //   and initialize the MPI world
    int thread_safety_provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &thread_safety_provided);
    MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI::ERRORS_THROW_EXCEPTIONS);
    //MPI_Status status;
    const double start_time = MPI_Wtime();

    int wRank=-1, wSize=-1;
    MPI_Comm_rank( MPI_COMM_WORLD, &wRank );
    MPI_Comm_size( MPI_COMM_WORLD, &wSize );

    // For the time being, hard-code the filter scales
    //   include scales as a comma-separated list
    //   scales are given in metres
    // A zero scale will cause everything to nan out
    std::vector<double> filter_scales {150e3};

    // Parse command-line flags
    char buffer [50];
    if (wRank == 0) {
        for(int ii = 1; ii < argc; ++ii) {  
            fprintf(stdout, "Argument %d : %s\n", ii, argv[ii]);
            snprintf(buffer, 50, "--version");
            if (strcmp(argv[ii], buffer) == 0) {
                print_compile_info(filter_scales);
                return 0;
            } else {
                fprintf(stderr, "Flag %s not recognized.\n", argv[ii]);
                return -1;
            }
        }
    }

    #if DEBUG >= 0
    if (wRank == 0) {
        fprintf(stdout, "\n\n");
        fprintf(stdout, "Compiled at %s on %s.\n", __TIME__, __DATE__);
        fprintf(stdout, "  Version %d.%d.%d \n", 
                MAJOR_VERSION, MINOR_VERSION, PATCH_VERSION);
        fprintf(stdout, "\n");
    }
    #endif

    #if DEBUG >= 0
    if (wRank == 0) {
        if (constants::CARTESIAN) { 
            fprintf(stdout, "Using spherical coordinates.\n");
        } else {
            fprintf(stdout, "Using Cartesian coordinates.\n");
        }
    }
    #endif
    MPI_Barrier(MPI_COMM_WORLD);

    // Print processor assignments
    #if DEBUG >= 2
    int tid, nthreads;
    const int max_threads = omp_get_max_threads();
    omp_set_num_threads( max_threads );
    #pragma omp parallel default(none) private(tid, nthreads) \
        shared(stdout) firstprivate(wRank, wSize)
    {
        tid = omp_get_thread_num();
        nthreads = omp_get_num_threads();
        fprintf(stdout, "Hello from thread %d of %d on processor %d of %d.\n", 
                tid+1, nthreads, wRank+1, wSize);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    #endif

    std::vector<double> longitude;
    std::vector<double> latitude;
    std::vector<double> time;
    std::vector<double> depth;

    std::vector<double> u_r;
    std::vector<double> u_lon;
    std::vector<double> u_lat;

    std::vector<double> rho;
    std::vector<double> p;

    std::vector<double> mask;
    std::vector<int> myCounts, myStarts;
    size_t II;

    // Read in source data / get size information
    #if DEBUG >= 1
    if (wRank == 0) { fprintf(stdout, "Reading in source data.\n\n"); }
    #endif

    // Read in the grid coordinates
    read_var_from_file(longitude, "longitude", "input.nc");
    read_var_from_file(latitude,  "latitude",  "input.nc");
    read_var_from_file(time,      "time",      "input.nc");
    read_var_from_file(depth,     "depth",     "input.nc");

    // Read in the velocity fields
    read_var_from_file(u_lon, "uo", "input.nc", &mask, &myCounts, &myStarts);
    read_var_from_file(u_lat, "vo", "input.nc", &mask, &myCounts, &myStarts);

    // No u_r in inputs, so initialize as zero
    u_r.resize(u_lon.size());
    #pragma omp parallel default(none) private(II) shared(u_r)
    { 
        #pragma omp for collapse(1) schedule(static)
        for (II = 0; II < u_r.size(); II++) {
            u_r.at(II) = 0.;
        }
    }

    if (constants::COMP_BC_TRANSFERS) {
        // If desired, read in rho and p
        read_var_from_file(rho, "rho", "input.nc");
        read_var_from_file(p,   "p",   "input.nc");
    }

    //const int Ntime  = time.size();
    //const int Ndepth = depth.size();
    const int Nlon   = longitude.size();
    const int Nlat   = latitude.size();

    if (not(constants::CARTESIAN)) {
        // Convert coordinate to radians
        int ii;
        #pragma omp parallel default(none) private(ii) shared(longitude)
        { 
            #pragma omp for collapse(1) schedule(static)
            for (ii = 0; ii < Nlon; ii++) {
                longitude.at(ii) = longitude.at(ii) * M_PI / 180;
            }
        }
        #pragma omp parallel default(none) private(ii) shared(latitude)
        {
            #pragma omp for collapse(1) schedule(static)
            for (ii = 0; ii < Nlat; ii++) {
                latitude.at(ii) = latitude.at(ii) * M_PI / 180;
            }
        }
    }

    // Compute the area of each 'cell'
    //   which will be necessary for integration
    #if DEBUG >= 1
    if (wRank == 0) { fprintf(stdout, "Computing the cell areas.\n\n"); }
    #endif

    std::vector<double> areas(Nlon * Nlat);
    compute_areas(areas, longitude, latitude);

    // Now pass the arrays along to the filtering routines
    const double pre_filter_time = MPI_Wtime();
    filtering(u_r, u_lon, u_lat,
              rho, p,
              filter_scales,
              areas, 
              time, depth,
              longitude, latitude,
              mask, myCounts, myStarts);
    const double post_filter_time = MPI_Wtime();

    // Done!
    #if DEBUG >= 0
    const double delta_clock = MPI_Wtick();
    if (wRank == 0) {
        fprintf(stdout, "\n\n");
        fprintf(stdout, "Process completed.\n");
        fprintf(stdout, "\n");
        fprintf(stdout, "Start-up time:  %g\n", pre_filter_time - start_time);
        fprintf(stdout, "Filtering time: %g\n", post_filter_time - pre_filter_time);
        fprintf(stdout, "   (clock resolution = %g)\n", delta_clock);
    }
    #endif

    fprintf(stdout, "Processor %d / %d waiting to finalize.\n", wRank + 1, wSize);
    MPI_Finalize();
    return 0;
}
