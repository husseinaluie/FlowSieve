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

#include "netcdf_io.hpp"
#include "functions.hpp"
#include "constants.hpp"

int main(int argc, char *argv[]) {
    
    if (constants::PERIODIC_Y) {
        // The non-uniform lat routine cannot handle
        //   periodic lat (y) grids
        assert(constants::UNIFORM_LAT_GRID);
    }

    static_assert(constants::CARTESIAN);
    static_assert(constants::PERIODIC_X);
    static_assert(constants::PERIODIC_Y);
    static_assert(constants::UNIFORM_LAT_GRID);

    // Enable all floating point exceptions but FE_INEXACT
    //feenableexcept(FE_ALL_EXCEPT & ~FE_INEXACT);

    // Specify the number of OpenMP threads
    //   and initialize the MPI world
    int thread_safety_provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &thread_safety_provided);
    MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI::ERRORS_THROW_EXCEPTIONS);
    //MPI_Status status;

    int wRank=-1, wSize=-1;
    MPI_Comm_rank( MPI_COMM_WORLD, &wRank );
    MPI_Comm_size( MPI_COMM_WORLD, &wSize );

    // For the time being, hard-code the filter scales
    //   include scales as a comma-separated list
    //   scales are given in metres
    // A zero scale will cause everything to nan out
    const double Lx = 25e3;
    std::vector<double> filter_scales { 10., 0.1*Lx, 0.25*Lx };

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

    // Print processor assignments
    MPI_Barrier(MPI_COMM_WORLD);
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
    MPI_Barrier(MPI_COMM_WORLD);
    #endif

    std::vector<double> longitude;
    std::vector<double> latitude;
    std::vector<double> time;
    std::vector<double> depth;

    std::vector<double> u, v, h;

    std::vector<double> mask;
    std::vector<int> myCounts, myStarts;

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
    read_var_from_file(u, "u", "input.nc", &mask, &myCounts, &myStarts);
    read_var_from_file(v, "v", "input.nc", &mask, &myCounts, &myStarts);
    read_var_from_file(h, "h", "input.nc", &mask, &myCounts, &myStarts);

    //const int Ntime  = time.size();
    //const int Ndepth = depth.size();
    const int Nlon   = longitude.size();
    const int Nlat   = latitude.size();

    // Compute the area of each 'cell'
    //   which will be necessary for integration
    #if DEBUG >= 1
    if (wRank == 0) { fprintf(stdout, "Computing the cell areas.\n\n"); }
    #endif

    std::vector<double> areas(Nlon * Nlat);
    compute_areas(areas, longitude, latitude);

    // Now pass the arrays along to the filtering routines
    filtering_sw(u, v, h,
              filter_scales, areas, 
              time, depth, longitude, latitude,
              mask, myCounts, myStarts);

    // Done!
    fprintf(stdout, "Processor %d / %d waiting to finalize.\n", wRank + 1, wSize);
    MPI_Finalize();
    return 0;
}
