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
    
    // Cannot extend to poles AND be Cartesian
    static_assert( not( (constants::EXTEND_DOMAIN_TO_POLES) and (constants::CARTESIAN) ),
            "Cartesian implies that there are no poles, so cannot extend to poles."
            "Please update constants.hpp accordingly.");

    // Enable all floating point exceptions but FE_INEXACT
    //feenableexcept( FE_ALL_EXCEPT & ~FE_INEXACT & ~FE_UNDERFLOW );
    //fprintf( stdout, " %d : %d \n", FE_ALL_EXCEPT, FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW | FE_INEXACT | FE_UNDERFLOW );
    feenableexcept( FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW );

    // Specify the number of OpenMP threads
    //   and initialize the MPI world
    int thread_safety_provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &thread_safety_provided);
    //MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI::ERRORS_THROW_EXCEPTIONS);
    const double start_time = MPI_Wtime();

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
    const bool asked_help = input.cmdOptionExists("--help");
    if (asked_help) {
        fprintf( stdout, "\033[1;4mThe command-line input arguments [and default values] are:\033[0m\n" );
    }

    // first argument is the flag, second argument is default value (for when flag is not present)
    const std::string   &Helm_input_fname  = input.getCmdOption("--Helmholtz_input_file",       
                                                                "Helmholtz_projection.nc",      
                                                                asked_help,
                                                                "netCDF file containing streamfunction and potential function."),
                        &adjacency_fname  = input.getCmdOption("--adjacency_file",     
                                                               "adjacency.nc",  
                                                               asked_help,
                                                               "Filename for the adjacency data (from the LLC_build_adjacency routine)"),
                        &vel_input_fname   = input.getCmdOption("--velocity_input_file",        
                                                                "vels.nc",                      
                                                                asked_help,
                                                                "netCDF file containing velocity field (only used to get land mask)");

    const std::string   &time_dim_name      = input.getCmdOption("--time",        
                                                                 "time",       
                                                                 asked_help,
                                                                 "Name of 'time' dimension in netCDF input file."),
                        &depth_dim_name     = input.getCmdOption("--depth",       
                                                                 "depth",      
                                                                 asked_help,
                                                                 "Name of 'depth' dimension in netCDF input file."),
                        &latitude_dim_name  = input.getCmdOption("--latitude",    
                                                                 "latitude",   
                                                                 asked_help,
                                                                 "Name of 'latitude' dimension in netCDF input file."),
                        &longitude_dim_name = input.getCmdOption("--longitude",   
                                                                 "longitude",  
                                                                 asked_help,
                                                                 "Name of 'longitude' dimension in netCDF input file.");

    const std::string &latlon_in_degrees  = input.getCmdOption("--is_degrees",   
                                                               "true", 
                                                               asked_help,
                                                               "Boolean (true/false) indicating if the grid is in degrees (true) or radians (false).");

    const std::string   &Nprocs_in_time_string  = input.getCmdOption("--Nprocs_in_time",  
                                                                     "1", 
                                                                     asked_help,
                                                                     "The number of MPI divisions in time. Optimally divides Ntime evenly.\nIf Ndepth = 1, Nprocs_in_time is automatically determined."),
                        &Nprocs_in_depth_string = input.getCmdOption("--Nprocs_in_depth", 
                                                                     "1", 
                                                                     asked_help,
                                                                     "The number of MPI divisions in depth. Optimally divides Ndepth evenly.\nIf Ntime = 1, Nprocs_in_depth is automatically determined.");
    const int   Nprocs_in_time_input  = stoi(Nprocs_in_time_string),
                Nprocs_in_depth_input = stoi(Nprocs_in_depth_string);

    const std::string   &tor_field_var_name     = input.getCmdOption("--tor_field",     "Psi",   asked_help, "Name of toroidal field (streamfunction) in input file."),
                        &pot_field_var_name     = input.getCmdOption("--pot_field",     "Phi",   asked_help, "Name of potential field (potential function) in input file."),
                        &dArea_field_var_name   = input.getCmdOption("--dArea_field",   "dA",    asked_help, "Name of cell areas in input file."),
                        &vel_field_var_name     = input.getCmdOption("--vel_field",     "u_lat", asked_help, "Name of a velocity field in input file (used to get land information).");

    const std::string   &region_defs_fname    = input.getCmdOption("--region_definitions_file",    
                                                                   "region_definitions.nc", 
                                                                   asked_help,
                                                                   "netCDF file containing user-specified region definitions."),
                        &region_defs_dim_name = input.getCmdOption("--region_definitions_dim",     
                                                                   "region",                
                                                                   asked_help,
                                                                   "Name of the region dimension in the regions file."),
                        &region_defs_var_name = input.getCmdOption("--region_definitions_var",     
                                                                   "region_definition",     
                                                                   asked_help,
                                                                   "Name of the variable in the regions file that provides the region definitions.");

    // Also read in the filter scales from the commandline
    //   e.g. --filter_scales "10.e3 150.76e3 1000e3" (units are in metres)
    std::vector<double> filter_scales;
    input.getFilterScales( filter_scales, "--filter_scales", asked_help );

    if (asked_help) { return 0; }

    // Print processor assignments
    const int max_threads = omp_get_max_threads();
    omp_set_num_threads( max_threads );

    // Print some header info, depending on debug level
    print_header_info();

    // Initialize dataset class instance
    dataset source_data, mask_data;

    // Read in source data / get size information
    #if DEBUG >= 1
    if (wRank == 0) { fprintf(stdout, "Reading in source data.\n\n"); }
    #endif

    // Read in the grid coordinates
    //   implicitely assume coordinates are the same between input files
    source_data.load_time(      time_dim_name,      Helm_input_fname );
    source_data.load_depth(     depth_dim_name,     Helm_input_fname );
    read_LLC_latlon_from_file( source_data.latitude,  latitude_dim_name,  Helm_input_fname );
    read_LLC_latlon_from_file( source_data.longitude, longitude_dim_name, Helm_input_fname );

    source_data.Nlat = 1;
    source_data.Nlon = 1;

    // Convert to radians, if appropriate
    if ( latlon_in_degrees == "true" ) {
        convert_coordinates( source_data.longitude, source_data.latitude );
    }

    // Load the adjacency matrix and other adjacency-adjacent arrays
    source_data.load_adjacency( adjacency_fname );

    // Cell areas are trickier, so they will be passed in as an input.
    // ... but, for right now, just treat them as 1. everywhere
    //source_data.areas.resize( source_data.latitude.size(), 1. );
    //source_data.load_variable( "dA", "dA", vel_input_fname, false, false);
    //source_data.areas = source_data.variables["dA"];
    read_LLC_latlon_from_file( source_data.areas, dArea_field_var_name, vel_input_fname );

    // Read in the toroidal and potential fields
    source_data.load_variable( "F_potential", pot_field_var_name, Helm_input_fname, false, true );
    source_data.load_variable( "F_toroidal",  tor_field_var_name, Helm_input_fname, false, true );

    // Get the MPI-local dimension sizes
    source_data.Ntime  = source_data.myCounts[0];
    source_data.Ndepth = source_data.myCounts[1];

    // Get mask : read in velocity to get the mask
    source_data.load_variable( "sample_velocity", vel_field_var_name, vel_input_fname, true, true);

    // Now pass the data along to the filtering routines
    const double pre_filter_time = MPI_Wtime();
    LLC_filtering_helmholtz( source_data, filter_scales );
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

    #if DEBUG >= 1
    fprintf(stdout, "Processor %d / %d waiting to finalize.\n", wRank + 1, wSize);
    #endif
    MPI_Finalize();
    return 0;
}
