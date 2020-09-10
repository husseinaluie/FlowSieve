#include <fenv.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <algorithm>
#include <math.h>
#include <vector>
#include <mpi.h>
#include <omp.h>

#include "../netcdf_io.hpp"
#include "../functions.hpp"
#include "../constants.hpp"
#include "../fft_based.hpp"

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
    std::vector<double> filter_scales {
        2, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100
    };

    //
    //// Parse command-line arguments
    //
    InputParser input(argc, argv);
    if(input.cmdOptionExists("--version")){
        if (wRank == 0) { print_compile_info(NULL); } 
        return 0;
    }

    // first argument is the flag, second argument is default value (for when flag is not present)
    const std::string &input_fname       = input.getCmdOption("--input_file",  "input.nc");

    const std::string &time_dim_name      = input.getCmdOption("--time",        "time");
    const std::string &depth_dim_name     = input.getCmdOption("--depth",       "depth");
    const std::string &latitude_dim_name  = input.getCmdOption("--latitude",    "latitude");
    const std::string &longitude_dim_name = input.getCmdOption("--longitude",   "longitude");

    const std::string &zonal_vel_name    = input.getCmdOption("--zonal_vel",   "uo");
    const std::string &merid_vel_name    = input.getCmdOption("--merid_vel",   "vo");

    // Set OpenMP thread number
    const int max_threads = omp_get_max_threads();
    omp_set_num_threads( max_threads );

    // Print some header info, depending on debug level
    print_header_info();

    std::vector<double> longitude, latitude, time, depth;
    std::vector<double> u_r, u_lon, u_lat, rho, p;
    std::vector<bool> mask;
    std::vector<int> myCounts, myStarts;
    size_t II;

    // Read in source data / get size information
    #if DEBUG >= 1
    if (wRank == 0) { fprintf(stdout, "Reading in source data.\n\n"); }
    #endif

    // Read in the grid coordinates
    read_var_from_file(time,      time_dim_name,      input_fname);
    read_var_from_file(depth,     depth_dim_name,     input_fname);
    read_var_from_file(latitude,  latitude_dim_name,  input_fname);
    read_var_from_file(longitude, longitude_dim_name, input_fname);
     
    convert_coordinates(longitude, latitude);

    // Read in the velocity fields
    read_var_from_file(u_lon, zonal_vel_name, input_fname, &mask, &myCounts, &myStarts);
    read_var_from_file(u_lat, merid_vel_name, input_fname, &mask, &myCounts, &myStarts);

    // No u_r in inputs, so initialize as zero
    u_r.resize(u_lon.size());
    #pragma omp parallel default(none) private(II) shared(u_r)
    { 
        #pragma omp for collapse(1) schedule(static)
        for (II = 0; II < u_r.size(); II++) {
            u_r.at(II) = 0.;
        }
    }

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

    // Compute the area of each 'cell'
    //   which will be necessary for integration
    #if DEBUG >= 1
    if (wRank == 0) { fprintf(stdout, "Computing the cell areas.\n\n"); }
    #endif

    std::vector<double> areas(Nlon * Nlat);
    compute_areas(areas, longitude, latitude);

    // Now pass the arrays along to the filtering routines
    const double pre_filter_time = MPI_Wtime();
    filtering_fftw(
            u_r, u_lon, u_lat,
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
        fprintf(stdout, "Start-up time  = %.13g\n", pre_filter_time - start_time);
        fprintf(stdout, "Filtering time = %.13g\n", post_filter_time - pre_filter_time);
        fprintf(stdout, "   (clock resolution = %.13g)\n", delta_clock);
    }
    #endif

    fprintf(stdout, "Processor %d / %d waiting to finalize.\n", wRank + 1, wSize);
    MPI_Finalize();
    return 0;
}
