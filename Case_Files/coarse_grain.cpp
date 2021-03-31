#include <fenv.h>
#include <stdio.h>
#include <string>
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
        // Michele 2d-turbulence
        //2, 5, 10, 15, 20, 30, 45, 65, 90, 125, 180

        // 2d fill-in
        //7, 12, 17, 24, 37, 54, 72, 81, 256, 500

        // 2d full-output
        //10, 20, 45, 65

        // AVISO full output
        //50e3, 100e3, 215e3, 460e3

        //10e3, 100e3, 250e3, 500e3, 1000e3

        // Artifical Dataset
        //0.1, 0.2, 0.3, 0.4, 0.5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12

        // AVISO paper
        1.e4, 1.29e4, 1.67e4, 2.15e4, 2.78e4, 3.59e4, 4.64e4, 5.99e4, 7.74e4,
        1.e5, 1.29e5, 1.67e5, 2.15e5, 2.78e5, 3.59e5, 4.64e5, 5.99e5, 7.74e5,
        1.e6, 1.29e6, 1.67e6, 2.15e6//, 2.78e6, 3.59e6

        // Lo-res (NEMO part of AVISO paper)
        /*
        1e4, 1.58e4, 2.51e4, 3.98e4, 6.31e4,
        1e5, 1.58e5, 2.51e5, 3.98e5, 6.31e5,
        1e6, 
        */

        // Double resolution of above (i.e. fill in gaps)
        /*
        1.26e4,  2.00e4  3.16e4,  5.01e4,  7.94e4,
        1.26e5,  2.00e5  3.16e5,  5.01e5,  7.94e5,
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
    const std::string &input_fname       = input.getCmdOption("--input_file",  "input.nc");

    const std::string   &time_dim_name      = input.getCmdOption("--time",        "time"),
                        &depth_dim_name     = input.getCmdOption("--depth",       "depth"),
                        &latitude_dim_name  = input.getCmdOption("--latitude",    "latitude"),
                        &longitude_dim_name = input.getCmdOption("--longitude",   "longitude");

    const std::string &latlon_in_degrees  = input.getCmdOption("--is_degrees",   "true");

    const std::string   &Nprocs_in_time_string  = input.getCmdOption("--Nprocs_in_time",  "1"),
                        &Nprocs_in_depth_string = input.getCmdOption("--Nprocs_in_depth", "1");
    const int   Nprocs_in_time_input  = stoi(Nprocs_in_time_string),
                Nprocs_in_depth_input = stoi(Nprocs_in_depth_string);

    const std::string   &zonal_vel_name    = input.getCmdOption("--zonal_vel",   "uo"),
                        &merid_vel_name    = input.getCmdOption("--merid_vel",   "vo"),
                        &density_var_name  = input.getCmdOption("--density",     "rho"),
                        &pressure_var_name = input.getCmdOption("--pressure",    "p");

    const std::string   &region_defs_fname    = input.getCmdOption("--region_definitions_file",    "region_definitions.nc"),
                        &region_defs_dim_name = input.getCmdOption("--region_definitions_dim",     "region"),
                        &region_defs_var_name = input.getCmdOption("--region_definitions_var",     "region_definition");

    // Set OpenMP thread number
    const int max_threads = omp_get_max_threads();
    omp_set_num_threads( max_threads );

    // Print some header info, depending on debug level
    print_header_info();

    // Initialize dataset class instance
    dataset source_data;

    // Read in source data / get size information
    #if DEBUG >= 1
    if (wRank == 0) { fprintf(stdout, "Reading in source data.\n\n"); }
    #endif

    // Read in the grid coordinates
    source_data.load_time(      time_dim_name,      input_fname );
    source_data.load_depth(     depth_dim_name,     input_fname );
    source_data.load_latitude(  latitude_dim_name,  input_fname );
    source_data.load_longitude( longitude_dim_name, input_fname );

    // Apply some cleaning to the processor allotments if necessary. 
    source_data.check_processor_divisions( Nprocs_in_time_input, Nprocs_in_depth_input );
     
    // Convert to radians, if appropriate
    if ( latlon_in_degrees == "true" ) {
        convert_coordinates( source_data.longitude, source_data.latitude );
    }

    // Compute the area of each 'cell' which will be necessary for integration
    source_data.compute_cell_areas();

    // Read in the velocity fields
    source_data.load_variable( "u_lon", zonal_vel_name, input_fname, true, true );
    source_data.load_variable( "u_lat", merid_vel_name, input_fname, true, true );

    // Get the MPI-local dimension sizes
    source_data.Ntime  = source_data.myCounts[0];
    source_data.Ndepth = source_data.myCounts[1];

    // No u_r in inputs, so initialize as zero
    source_data.variables.insert( std::pair< std::string, std::vector<double> >
                                           ( "u_r",       std::vector<double>(0, source_data.variables.at("u_lon").size()) ) 
                                );

    if (constants::COMP_BC_TRANSFERS) {
        // If desired, read in rho and p
        source_data.load_variable( "rho", density_var_name, input_fname, false, false );
        source_data.load_variable( "p", pressure_var_name, input_fname, false, false );
    }

    // Mask out the pole, if necessary (i.e. set lat = 90 to land)
    mask_out_pole( source_data.latitude, source_data.mask, source_data.Ntime, source_data.Ndepth, source_data.Nlat, source_data.Nlon );

    // Read in the region definitions and compute region areas
    source_data.load_region_definitions( region_defs_fname, region_defs_dim_name, region_defs_var_name );

    // Now pass the arrays along to the filtering routines
    const double pre_filter_time = MPI_Wtime();
    filtering( source_data, filter_scales );
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
