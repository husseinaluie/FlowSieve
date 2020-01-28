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
        10e3
    };

    // Parse command-line flags
    char buffer [50];
    if (wRank == 0) {
        for(int ii = 1; ii < argc; ++ii) {  
            fprintf(stdout, "Argument %d : %s\n", ii, argv[ii]);

            // check if we passed the '--version falg'
            snprintf(buffer, 50, "--version");
            if (strcmp(argv[ii], buffer) == 0) {
                print_compile_info(filter_scales);
                return 0;

            // otherwise, flag unrecognized
            } else {
                fprintf(stderr, "Flag %s not recognized.\n", argv[ii]);
                return -1;
            }
        }
    }

    char input_fname [50];
    snprintf(input_fname, 50, "input.nc");

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

    // Load in the variables that we want to filter
    std::vector<const std::vector<double>* > fields_to_filter;
    std::vector<std::string> var_names;

    std::vector<double> eta;
    read_var_from_file(eta, "eta", input_fname, &mask, &myCounts, &myStarts);
    fields_to_filter.push_back(&eta);
    var_names.push_back("eta");


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
    filter_fields(
            fields_to_filter, var_names,    
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
