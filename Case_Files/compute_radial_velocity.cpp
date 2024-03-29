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
                                                                "netCDF file where Helmholtz scalars are stored."),
                        &vel_input_fname   = input.getCmdOption("--velocity_input_file",        
                                                                "vels.nc",                      
                                                                asked_help,
                                                                "netCDF file where velocities are stored (to provide land information)."),
                        &output_fname      = input.getCmdOption("--output_file",                
                                                                "radial_vel.nc",                
                                                                asked_help,
                                                                "Filename for where the radial velocity should be stored.");

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
                                                                     "The number of MPI divisions in time. Optimally divides Ntime evenly.\nIf Ndepth = 1, Nprocs_in_time is automatically determined.");
    const int   Nprocs_in_time_input  = stoi(Nprocs_in_time_string),
                Nprocs_in_depth_input = 1;

    const std::string   &tor_field_var_name     = input.getCmdOption("--tor_field",     
                                                                     "Psi",          
                                                                     asked_help,
                                                                     "Name of the streamfunction (Psi) in the input file."),
                        &pot_field_var_name     = input.getCmdOption("--pot_field",     
                                                                     "Phi",          
                                                                     asked_help,
                                                                     "Name of the potential function (Phi) in the input file."),
                        &vel_field_var_name     = input.getCmdOption("--vel_field",     
                                                                     "u_lat",        
                                                                     asked_help,
                                                                     "Name of a velocity field in the velocity input file.");

    if (asked_help) { return 0; }

    // Print processor assignments
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
    //   implicitely assume coordinates are the same between input files
    source_data.load_time(      time_dim_name,      Helm_input_fname );
    source_data.load_depth(     depth_dim_name,     Helm_input_fname );
    source_data.load_latitude(  latitude_dim_name,  Helm_input_fname );
    source_data.load_longitude( longitude_dim_name, Helm_input_fname );

    // Apply some cleaning to the processor allotments if necessary. 
    source_data.check_processor_divisions( Nprocs_in_time_input, Nprocs_in_depth_input );
     
    // Convert to radians, if appropriate
    if ( latlon_in_degrees == "true" ) {
        convert_coordinates( source_data.longitude, source_data.latitude );
    }

    // Compute the area of each 'cell' which will be necessary for integration
    source_data.compute_cell_areas();

    // Read in the toroidal and potential fields
    source_data.load_variable( "F_potential", pot_field_var_name, Helm_input_fname, false, true );
    source_data.load_variable( "F_toroidal",  tor_field_var_name, Helm_input_fname, false, true );

    // Get the MPI-local dimension sizes
    source_data.Ntime  = source_data.myCounts[0];
    source_data.Ndepth = source_data.myCounts[1];

    source_data.compute_cell_areas();

    // Mask out the pole, if necessary (i.e. set lat = 90 to land)
    mask_out_pole( source_data.latitude, source_data.mask, source_data.Ntime, source_data.Ndepth, source_data.Nlat, source_data.Nlon );

    const std::vector<double> &F_potential = source_data.variables.at("F_potential");
    std::vector<double> null_vector(0);
    const int   Ntime   = source_data.Ntime,    // this is the MPI-local Ntime, not the full Ntime
                Ndepth  = source_data.Ndepth,   // this is the MPI-local Ndepth, not the full Ndepth
                Nlat    = source_data.Nlat,
                Nlon    = source_data.Nlon;
    const std::vector<double>   &time       = source_data.time,
                                &depth      = source_data.depth,
                                &latitude   = source_data.latitude,
                                &longitude  = source_data.longitude,
                                &dAreas     = source_data.areas;
    const unsigned int num_pts = Ntime * Ndepth * Nlat * Nlon;
    std::vector<double> negative_divergence( num_pts, 0. ), 
                        u_r( num_pts, 0. ), 
                        u_lon_pot( num_pts, 0. ),
                        u_lat_pot( num_pts, 0. );
    const std::vector<bool> mask( num_pts, true );
    source_data.mask = mask;

    // We only need the potential velocity, since div(tor vel) = 0 by definition
    #if DEBUG >= 2
    if (wRank == 0) { fprintf(stdout, "\nGet potential velocity.\n"); }
    #endif
    potential_vel_from_F( u_lon_pot, u_lat_pot, F_potential, longitude, latitude, Ntime, Ndepth, Nlat, Nlon, mask );

    #if DEBUG >= 2
    if (wRank == 0) { fprintf(stdout, "\nComputing the divergence\n"); }
    fflush(stdout);
    #endif
    compute_vorticity(
            null_vector, null_vector, null_vector, negative_divergence, null_vector,
            source_data, u_r, u_lon_pot, u_lat_pot );
    int Itime, Idepth, Ilat, Ilon;
    for ( size_t index = 0; index < negative_divergence.size(); ++index) {
        negative_divergence.at(index) = -1. * negative_divergence.at(index);
    }

    #if DEBUG >= 2
    if (wRank == 0) { fprintf(stdout, "\nStarting depth integral.\n"); }
    #endif
    depth_integrate( u_r, negative_divergence, source_data );

    // And write to a file
    std::vector<std::string> vars_to_write;
    vars_to_write.push_back("u_r");
    initialize_output_file( source_data, vars_to_write, output_fname.c_str(), -1);
    const std::vector<int>  &myCounts = source_data.myCounts,
                            &myStarts = source_data.myStarts;
    const int ndims = 4;
    size_t starts[ndims] = { size_t(myStarts.at(0)), size_t(myStarts.at(1)), size_t(myStarts.at(2)), size_t(myStarts.at(3)) };
    size_t counts[ndims] = { size_t(Ntime),          size_t(Ndepth),         size_t(Nlat),           size_t(Nlon)           };
    write_field_to_output(u_r, "u_r", starts, counts, output_fname, NULL);

    // Done!
    #if DEBUG >= 0
    const double delta_clock = MPI_Wtick();
    if (wRank == 0) {
        fprintf(stdout, "\n\n");
        fprintf(stdout, "Process completed.\n");
        fprintf(stdout, "\n");
    }
    #endif

    #if DEBUG >= 1
    fprintf(stdout, "Processor %d / %d waiting to finalize.\n", wRank + 1, wSize);
    #endif
    MPI_Finalize();
    return 0;
}
