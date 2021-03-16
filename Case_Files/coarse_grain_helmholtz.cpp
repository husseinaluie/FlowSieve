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
        100e3, 250e3, 400e3

        //1.e4 , 1.29e4, 1.67e4, 2.15e4, 2.78e4, 3.59e4,

        /*
                                                      4.64e4, 5.99e4, 7.74e4,
        1.e5, 1.29e5, 1.67e5, 2.15e5, 2.78e5, 3.59e5, 4.64e5, 5.99e5, 7.74e5,
        1.e6, 1.29e6, 1.67e6, 2.15e6, 
        */
        
        /*
        1.e4 , 1.29e4, 1.67e4, 2.15e4, 2.78e4, 3.59e4, 4.64e4, 5.99e4, 7.74e4,
        1.e5 , 1.29e5, 1.67e5, 2.15e5, 2.78e5, 3.59e5, 4.64e5, 5.99e5, 7.74e5,
        1.e6 , 1.29e6, 1.67e6, 2.15e6, 2.78e6, 3.59e6, 4.64e6, 5.99e6, 7.74e6,
        1.e7
        */

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
    const std::string &tor_input_fname = input.getCmdOption("--toroidal_input_file",   "toroidal_projection.nc");
    const std::string &pot_input_fname = input.getCmdOption("--potential_input_file",  "potential_projection.nc");
    const std::string &vel_input_fname = input.getCmdOption("--velocity_input_file",   "toroidal_projection.nc");

    const std::string &time_dim_name      = input.getCmdOption("--time",        "time");
    const std::string &depth_dim_name     = input.getCmdOption("--depth",       "depth");
    const std::string &latitude_dim_name  = input.getCmdOption("--latitude",    "latitude");
    const std::string &longitude_dim_name = input.getCmdOption("--longitude",   "longitude");

    const std::string &Nprocs_in_time_string  = input.getCmdOption("--Nprocs_in_time",  "1");
    const std::string &Nprocs_in_depth_string = input.getCmdOption("--Nprocs_in_depth", "1");
    const int Nprocs_in_time_input  = stoi(Nprocs_in_time_string);
    const int Nprocs_in_depth_input = stoi(Nprocs_in_depth_string);

    const std::string &tor_field_var_name = input.getCmdOption("--tor_field",   "F");
    const std::string &pot_field_var_name = input.getCmdOption("--pot_field",   "F");
    const std::string &vel_field_var_name = input.getCmdOption("--vel_field",   "u_lat");

    // Print processor assignments
    const int max_threads = omp_get_max_threads();
    omp_set_num_threads( max_threads );

    // Print some header info, depending on debug level
    print_header_info();

    std::vector<double> longitude, latitude, time, depth;
    std::vector<double> F_potential, F_toroidal, u_tor;
    std::vector<bool> mask;
    std::vector<int> myCounts, myStarts;

    // Read in source data / get size information
    #if DEBUG >= 1
    if (wRank == 0) { fprintf(stdout, "Reading in source data.\n\n"); }
    #endif

    // Read in the grid coordinates
    //   implicitely assume coordinates are the same between input files
    read_var_from_file(longitude, longitude_dim_name, tor_input_fname);
    read_var_from_file(latitude,  latitude_dim_name,  tor_input_fname);
    read_var_from_file(time,      time_dim_name,      tor_input_fname);
    read_var_from_file(depth,     depth_dim_name,     tor_input_fname);

    const int Ntime  = time.size();
    const int Ndepth = depth.size();
    const int Nlon   = longitude.size();
    const int Nlat   = latitude.size();

    // Apply some cleaning to the processor allotments if necessary. 
    const int Nprocs_in_time  = ( Ntime  == 1 ) ? 1 : 
                                ( Ndepth == 1 ) ? wSize : 
                                                  Nprocs_in_time_input;
    const int Nprocs_in_depth = ( Ndepth == 1 ) ? 1 : 
                                ( Ntime  == 1 ) ? wSize : 
                                                  Nprocs_in_depth_input;
    #if DEBUG >= 0
    if (Nprocs_in_time != Nprocs_in_time_input) { 
        if (wRank == 0) { fprintf(stdout, " WARNING!! Changing number of processors in time to %'d from %'d\n", Nprocs_in_time, Nprocs_in_time_input); }
    }
    if (Nprocs_in_depth != Nprocs_in_depth_input) { 
        if (wRank == 0) { fprintf(stdout, " WARNING!! Changing number of processors in depth to %'d from %'d\n", Nprocs_in_depth, Nprocs_in_depth_input); }
    }
    if (wRank == 0) { fprintf(stdout, " Nproc(time, depth) = (%'d, %'d)\n", Nprocs_in_time, Nprocs_in_depth); }
    #endif
    assert( Nprocs_in_time * Nprocs_in_depth == wSize );
     
    convert_coordinates(longitude, latitude);

    // Read in the toroidal and potential fields
    read_var_from_file(F_potential, pot_field_var_name, pot_input_fname, NULL, &myCounts, &myStarts, Nprocs_in_time, Nprocs_in_depth);
    read_var_from_file(F_toroidal,  tor_field_var_name, tor_input_fname, NULL, &myCounts, &myStarts, Nprocs_in_time, Nprocs_in_depth);

    // read in velocity to get the mask
    read_var_from_file(u_tor, vel_field_var_name, vel_input_fname, &mask, &myCounts, &myStarts, Nprocs_in_time, Nprocs_in_depth);

    // Mask out the pole, if necessary (i.e. set lat = 90 to land)
    mask_out_pole(latitude, mask, myCounts[0], myCounts[1], myCounts[2], myCounts[3]);

    // Compute the area of each 'cell'
    //   which will be necessary for integration
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
