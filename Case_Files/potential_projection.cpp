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
#include "../preprocess.hpp"

int main(int argc, char *argv[]) {
    
    // PERIODIC_Y implies UNIFORM_LAT_GRID
    static_assert( (constants::UNIFORM_LAT_GRID) or (not(constants::PERIODIC_Y)),
            "PERIODIC_Y requires UNIFORM_LAT_GRID.\n"
            "Please update constants.hpp accordingly.\n");
    static_assert( not(constants::CARTESIAN),
            "Potential projection not set to handle Cartesian coordinates.\n"
            );

    // Enable all floating point exceptions but FE_INEXACT
    //feenableexcept(FE_ALL_EXCEPT & ~FE_INEXACT);

    // Specify the number of OpenMP threads
    //   and initialize the MPI world
    int thread_safety_provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &thread_safety_provided);
    //MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI::ERRORS_THROW_EXCEPTIONS);

    int wRank=-1, wSize=-1;
    MPI_Comm_rank( MPI_COMM_WORLD, &wRank );
    MPI_Comm_size( MPI_COMM_WORLD, &wSize );

    //
    //// Parse command-line arguments
    //
    InputParser input(argc, argv);
    if(input.cmdOptionExists("--version")){
        if (wRank == 0) { print_compile_info(NULL); } 
        return 0;
    }

    // first argument is the flag, second argument is default value (for when flag is not present)
    const std::string &input_fname  = input.getCmdOption("--input_file",   "input.nc");
    const std::string &output_fname = input.getCmdOption("--output_file",  "potential_projection.nc");
    const std::string &seed_fname   = input.getCmdOption("--seed_file",    "seed.nc");

    const std::string &time_dim_name      = input.getCmdOption("--time",        "time");
    const std::string &depth_dim_name     = input.getCmdOption("--depth",       "depth");
    const std::string &latitude_dim_name  = input.getCmdOption("--latitude",    "latitude");
    const std::string &longitude_dim_name = input.getCmdOption("--longitude",   "longitude");

    const std::string &zonal_vel_name    = input.getCmdOption("--zonal_vel",   "uo");
    const std::string &merid_vel_name    = input.getCmdOption("--merid_vel",   "vo");

    const std::string &tolerance_string = input.getCmdOption("--tolerance", "5e-3");
    const double tolerance = stod(tolerance_string);  

    const std::string &iteration_string = input.getCmdOption("--max_iterations", "100000");
    const int max_iterations = stoi(iteration_string);  

    // Print processor assignments
    const int max_threads = omp_get_max_threads();
    omp_set_num_threads( max_threads );

    // Print some header info, depending on debug level
    print_header_info();

    std::vector<double> longitude, latitude, time, depth;
    std::vector<double> u_lon, u_lat, seed;
    std::vector<bool> mask;
    std::vector<int> myCounts, myStarts;

    // Read in source data / get size information
    #if DEBUG >= 1
    if (wRank == 0) { fprintf(stdout, "Reading in source data.\n\n"); }
    #endif

    // Read in the grid coordinates
    read_var_from_file(longitude, longitude_dim_name, input_fname);
    read_var_from_file(latitude,  latitude_dim_name,  input_fname);
    read_var_from_file(time,      time_dim_name,      input_fname);
    read_var_from_file(depth,     depth_dim_name,     input_fname);
     
    convert_coordinates(longitude, latitude);

    // Read in the velocity fields
    read_var_from_file(u_lon, zonal_vel_name,  input_fname, &mask, &myCounts, &myStarts);
    read_var_from_file(u_lat, merid_vel_name,  input_fname, &mask, &myCounts, &myStarts);

    // Read in the seed
    double seed_count;
    read_attr_from_file(seed_count, "seed_count", seed_fname);
    const bool single_seed = (seed_count < 2);
    if (wRank == 0) { fprintf(stdout, " single_seed = %s\n", single_seed ? "true" : "false"); }
    read_var_from_file(seed, "seed",  seed_fname, NULL, NULL, NULL, not(single_seed));

    const int Nlon   = longitude.size();
    const int Nlat   = latitude.size();

    // Compute the area of each 'cell'
    //   which will be necessary for integration
    #if DEBUG >= 1
    if (wRank == 0) { fprintf(stdout, "Computing the cell areas.\n\n"); }
    #endif
    std::vector<double> areas(Nlon * Nlat);
    compute_areas(areas, longitude, latitude);

    Apply_Potential_Projection(
            output_fname,
            u_lon, u_lat, time, depth, latitude, longitude,
            areas, mask, myCounts, myStarts, seed, single_seed,
            tolerance, max_iterations
            );


    // Done!
    #if DEBUG >= 0
    if (wRank == 0) {
        fprintf(stdout, "\n\n");
        fprintf(stdout, "Process completed.\n");
        fprintf(stdout, "\n");
    }
    #endif

    fprintf(stdout, "Processor %d / %d waiting to finalize.\n", wRank + 1, wSize);
    MPI_Finalize();
    return 0;
}
