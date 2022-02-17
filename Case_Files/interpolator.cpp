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
            "Toroidal projection now set to handle Cartesian coordinates.\n"
            );

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
    const std::string   &input_fname    = input.getCmdOption("--input_file",   "input.nc"),
                        &output_fname   = input.getCmdOption("--output_file",  "output.nc");

    const std::string   &time_dim_name      = input.getCmdOption("--time",        "time"),
                        &depth_dim_name     = input.getCmdOption("--depth",       "depth"),
                        &latitude_dim_name  = input.getCmdOption("--latitude",    "latitude"),
                        &longitude_dim_name = input.getCmdOption("--longitude",   "longitude");

    const std::string &latlon_in_degrees  = input.getCmdOption("--is_degrees",   "true");

    const std::string   &Nprocs_in_time_string  = input.getCmdOption("--Nprocs_in_time",  "1"),
                        &Nprocs_in_depth_string = input.getCmdOption("--Nprocs_in_depth", "1"),
                        &Nlayers_string         = input.getCmdOption("--Nlayers", "7");
    const int   Nprocs_in_time_input  = stoi(Nprocs_in_time_string),
                Nprocs_in_depth_input = stoi(Nprocs_in_depth_string),
                num_interp_layers     = stoi(Nlayers_string);

    std::vector< std::string > vars_to_refine, vars_in_output;
    input.getListofStrings( vars_to_refine, "--input_variables" );
    input.getListofStrings( vars_in_output, "--output_variables" );
    const int Nvars = vars_to_refine.size();

    // Print processor assignments
    const int max_threads = omp_get_max_threads();
    omp_set_num_threads( max_threads );

    // Print some header info, depending on debug level
    print_header_info();

    // Initialize dataset class instance
    dataset orig_data;

    // Read in source data / get size information
    #if DEBUG >= 1
    if (wRank == 0) { fprintf(stdout, "Reading in source data.\n\n"); }
    #endif

    // Read in the grid coordinates
    orig_data.load_time(      time_dim_name,      input_fname );
    orig_data.load_depth(     depth_dim_name,     input_fname );
    orig_data.load_latitude(  latitude_dim_name,  input_fname );
    orig_data.load_longitude( longitude_dim_name, input_fname );

    // Apply some cleaning to the processor allotments if necessary. 
    orig_data.check_processor_divisions( Nprocs_in_time_input, Nprocs_in_depth_input );
    orig_data.Nprocs_in_time  = orig_data.Nprocs_in_time;
    orig_data.Nprocs_in_depth = orig_data.Nprocs_in_depth;
     
    // Convert to radians, if appropriate
    if ( latlon_in_degrees == "true" ) {
        convert_coordinates( orig_data.longitude,   orig_data.latitude );
    }

    // Read in field to get size info
    orig_data.load_variable( "to_interp", vars_to_refine.at(0), input_fname, true, true );

    const int   full_Ntime  = orig_data.full_Ntime,
                Ntime       = orig_data.myCounts[0],
                Ndepth      = orig_data.myCounts[1],
                Nlat        = orig_data.Nlat,
                Nlon        = orig_data.Nlon;

    size_t starts[4] = { orig_data.myStarts.at(0), orig_data.myStarts.at(1), 0,    0    };
    size_t counts[4] = { orig_data.myCounts.at(0), orig_data.myCounts.at(1), Nlat, Nlon };

    // Compute the area of each 'cell' which will be necessary for creating the output file
    if (wRank == 0) { fprintf( stdout, "Computing cell areas.\n" ); }
    orig_data.compute_cell_areas();

    // Initialize file and write out coarsened fields
    if (wRank == 0) { fprintf( stdout, "Preparing output file\n" ); }
    initialize_output_file( orig_data, vars_in_output, output_fname.c_str() );

    // Now coarsen the velocity fields
    const size_t Npts = Ntime * Ndepth * Nlat * Nlon;
    std::vector<double> interped_field(Npts);

    // Next, the coarse velocities
    for ( int Ivar = 0; Ivar < Nvars; Ivar++ ) {

        orig_data.load_variable( "to_interp", vars_to_refine.at(Ivar), input_fname, true, true );

        interpolate_over_land_from_coast( interped_field, orig_data.variables["to_interp"], num_interp_layers,
                orig_data.time, orig_data.depth, orig_data.latitude, orig_data.longitude, orig_data.mask, orig_data.myCounts);

        write_field_to_output( interped_field, vars_in_output.at(Ivar), starts, counts, output_fname );
    }

    #if DEBUG >= 1
    fprintf(stdout, "Processor %d / %d waiting to finalize.\n", wRank + 1, wSize);
    #endif
    MPI_Finalize();
    return 0;
}
