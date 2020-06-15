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

int main(int argc, char *argv[]) {
    
    // PERIODIC_Y implies UNIFORM_LAT_GRID
    static_assert( (constants::UNIFORM_LAT_GRID) or (not(constants::PERIODIC_Y)),
            "PERIODIC_Y requires UNIFORM_LAT_GRID.\n"
            "Please update constants.hpp accordingly.\n");

    // NO_FULL_OUTPUTS implies APPLY_POSTPROCESS
    static_assert( (constants::APPLY_POSTPROCESS) or (not(constants::NO_FULL_OUTPUTS)),
            "If NO_FULL_OUTPUTS is true, then APPLY_POSTPROCESS must also be true, "
            "otherwise no outputs will be produced.\n"
            "Please update constants.hpp accordingly.");

    // NO_FULL_OUTPUTS implies MINIMAL_OUTPUT
    static_assert( (constants::MINIMAL_OUTPUT) or (not(constants::NO_FULL_OUTPUTS)),
            "NO_FULL_OUTPUTS implies MINIMAL_OUTPUT. "
            "You must either change NO_FULL_OUTPUTS to false, "
            "or MINIMAL_OUTPUT to true.\n" 
            "Please update constants.hpp accordingly.");

    // Enable all floating point exceptions but FE_INEXACT
    //feenableexcept(FE_ALL_EXCEPT & ~FE_INEXACT);

    // Specify the number of OpenMP threads
    //   and initialize the MPI world
    int thread_safety_provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &thread_safety_provided);
    //MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI::ERRORS_THROW_EXCEPTIONS);
    const double start_time = MPI_Wtime();

    int wRank=-1, wSize=-1;
    MPI_Comm_rank( MPI_COMM_WORLD, &wRank );
    MPI_Comm_size( MPI_COMM_WORLD, &wSize );

    // For the time being, hard-code the filter scales
    //   include scales as a comma-separated list
    //   scales are given in metres
    // A zero scale will cause everything to nan out
    std::vector<double> filter_scales { 
        1e3,
        1e4, 1.58e4, 2.51e4, 3.98e4, 6.31e4,
        1e5, 1.58e5, 2.51e5, 3.98e5, 6.31e5,
        1e6, //1.58e6, 2.51e6, 3.98e6, 6.31e6,
        //1e7//, 1.58e7, 2.51e7, 3.98e7, 6.31e7
    };

    // Default input file name
    char input_fname[50], flag_version[50], flag_input[50];
    snprintf(input_fname,  50, "input.nc");
    snprintf(flag_version, 50, "--version");
    snprintf(flag_input,   50, "--input_file");

    // Parse command-line flags
    int ii = 1;
    //for(int ii = 1; ii < argc; ++ii) {  
    while(ii < argc) {  
        if (wRank == 0) {
            fprintf(stdout, "Argument %d : %s\n", ii, argv[ii]);
        }

        // check if the flag is 'version'
        if (strcmp(argv[ii], flag_version) == 0) {
            if (wRank == 0) {
                print_compile_info(&filter_scales);
            }
            return 0;

            // check if the flag is 'input_file' and, if it is
            //   the next input is then used as the filename
            //   of the input
        } else if (strcmp(argv[ii], flag_input) == 0) {
            snprintf(input_fname, 50, argv[ii+1]);
            ++ii;

            // Otherwise, the flag is unrecognized
        } else {
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
    std::vector<double> F_potential, F_toroidal, u_tor;
    std::vector<double> mask;
    std::vector<int> myCounts, myStarts;

    // Read in source data / get size information
    #if DEBUG >= 1
    if (wRank == 0) { fprintf(stdout, "Reading in source data.\n\n"); }
    #endif

    // Read in the grid coordinates
    read_var_from_file(longitude, "longitude", input_fname);
    read_var_from_file(latitude,  "latitude",  input_fname);
    read_var_from_file(time,      "time",      input_fname);
    read_var_from_file(depth,     "depth",     input_fname);

    // Read in the toroidal and potential fields
    read_var_from_file(F_potential, "F_potential", 
            input_fname, NULL, &myCounts, &myStarts);
    read_var_from_file(F_toroidal,  "F_toroidal", 
            input_fname, NULL, &myCounts, &myStarts);

    // read in velocity to get the mask
    read_var_from_file(u_tor, "u_tor", 
            input_fname, &mask, &myCounts, &myStarts);

    #if DEBUG >= 1
    const int Ntime  = time.size();
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
            for (ii = 0; ii < Nlon; ii++) {
                longitude.at(ii) = longitude.at(ii) * D2R;
            }

            #pragma omp for collapse(1) schedule(static)
            for (ii = 0; ii < Nlat; ii++) {
                latitude.at(ii) = latitude.at(ii) * D2R;
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
    filtering_helmholtz(
            F_potential, F_toroidal,
            filter_scales, areas, 
            time, depth, longitude, latitude,
            mask, myCounts, myStarts);
    const double post_filter_time = MPI_Wtime();

    // Done!
    #if DEBUG >= 0
    const double delta_clock = MPI_Wtick();
    if (wRank == 0) {
        fprintf(stdout, "\n\n");
        fprintf(stdout, "Process completed.\n");
        fprintf(stdout, "\n");
        fprintf(stdout, "Start-up time  = %.13g\n", pre_filter_time - start_time);
        fprintf(stdout, "Filtering time = %.13g\n", post_filter_time - pre_filter_time);
        fprintf(stdout, "   (clock resolution = %.13g)\n", delta_clock);
    }
    #endif

    fprintf(stdout, "Processor %d / %d waiting to finalize.\n", wRank + 1, wSize);
    MPI_Finalize();
    return 0;
}
