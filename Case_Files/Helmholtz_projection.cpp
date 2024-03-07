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

    // Currently cannot be Cartesian
    static_assert( not(constants::CARTESIAN),
            "Toroidal projection not set to handle Cartesian coordinates.\n"
            );

    // Enable all floating point exceptions but FE_INEXACT
    feenableexcept( FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW );

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
    const bool asked_help = input.cmdOptionExists("--help");
    if (asked_help) {
        fprintf( stdout, "\033[1;4mThe command-line input arguments [and default values] are:\033[0m\n" );
    }

    // first argument is the flag, second argument is default value (for when flag is not present)
    const std::string   &input_fname      = input.getCmdOption("--input_file",      
                                                               "input.nc",                 
                                                               asked_help,
                                                               "netCDF file containing vector field to Helmholtz decompose"),
                        &output_fname     = input.getCmdOption("--output_file",     
                                                               "projection_Helmholtz.nc",  
                                                               asked_help,
                                                               "Filename for the output (netCDF file of Helmholtz scalars)"),
                        &seed_fname       = input.getCmdOption("--seed_file",       
                                                               "zero",                  
                                                               asked_help,
                                                               "netCDF file containing initial guesses for Helmholtz scalars.\nUse 'zero' (the default) for no seed.");

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

    const std::string   &zonal_vel_name    = input.getCmdOption("--zonal_vel",   
                                                                "uo",          
                                                                asked_help,
                                                                "Name of zonal (eastward) velocity in netCDF input file"),
                        &merid_vel_name    = input.getCmdOption("--merid_vel",   
                                                                "vo",          
                                                                asked_help,
                                                                "Name of meridional (northward) velocity in netCDF input file"),
                        &tor_seed_name     = input.getCmdOption("--tor_seed",    
                                                                "Psi_seed",    
                                                                asked_help,
                                                                "Name of streamfunction (Psi) seed variable in seed file (if applicable)"),
                        &pot_seed_name     = input.getCmdOption("--pot_seed",    
                                                                "Phi_seed",    
                                                                asked_help,
                                                                "Name of potential function (Phi) seed variable in seed file (if applicable)");

    const std::string &tolerance_string = input.getCmdOption("--tolerance", 
                                                             "1e-120", 
                                                             asked_help,
                                                             "Termination tolerance. Note that this is based on matrix norms, and not a great measure of convergence.\nIn general, set to be very small (default) and check the convergence metrics\nthat are provided in the output file.");
    const double tolerance = stod(tolerance_string);  

    const std::string &iteration_string = input.getCmdOption("--max_iterations", 
                                                             "100000", 
                                                             asked_help,
                                                             "Maximum number of iterations for the solver before terminating. Can use exponential notation (e.g. '5e3')");
    const int max_iterations = stod(iteration_string);  

    const std::string &Tikhov_Lap_string = input.getCmdOption("--Tikhov_Laplace", "1.", asked_help);
    const double Tikhov_Laplace = stod(Tikhov_Lap_string);  

    const std::string &use_mask_string = input.getCmdOption("--use_mask", 
                                                            "false", 
                                                            asked_help,
                                                            "Boolean (true/false) indicating if land masking should be accounted for in the projection.\nThis is generally not advised, especially if coarse-graining will use 'filter over land' anyways.");
    const bool use_mask = string_to_bool(use_mask_string);

    const std::string &use_area_weight_string = input.getCmdOption("--use_area_weight", 
                                                                   "true", 
                                                                   asked_help,
                                                                   "Boolean (true/false) indicating if the least-squares problem should be weighted by area, so that larger cells have more priority.\nSetting to true is generally advised because of the poles.");
    const bool use_area_weight = string_to_bool(use_area_weight_string);

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
    source_data.load_time(      time_dim_name,      input_fname );
    source_data.load_depth(     depth_dim_name,     input_fname );
    source_data.load_latitude(  latitude_dim_name,  input_fname );
    source_data.load_longitude( longitude_dim_name, input_fname );

    // Apply some cleaning to the processor allotments if necessary. 
    source_data.check_processor_divisions( Nprocs_in_time_input, Nprocs_in_depth_input );
     
    // Convert to radians, if appropriate
    if ( (latlon_in_degrees == "true") and (not(constants::CARTESIAN)) ) {
        convert_coordinates( source_data.longitude, source_data.latitude );
    }

    // Read in the velocity fields
    source_data.load_variable( "u_lon", zonal_vel_name, input_fname, true, true );
    source_data.load_variable( "u_lat", merid_vel_name, input_fname, true, true );

    // Get the MPI-local dimension sizes
    source_data.Ntime  = source_data.myCounts[0];
    source_data.Ndepth = source_data.myCounts[1];

    //
    //// If necessary, extend the domain to reach the poles
    //
    if ( constants::EXTEND_DOMAIN_TO_POLES ) {
        #if DEBUG >= 0
        if (wRank == 0) { fprintf( stdout, "Extending the domain to the poles\n" ); }
        #endif

        // Extend the latitude grid to reach the poles and update source_data with the new info.
        std::vector<double> extended_latitude;
        int orig_lat_start_in_extend;
        #if DEBUG >= 1
        if (wRank == 0) { fprintf( stdout, "    Extending latitude to poles\n" ); }
        #endif
        extend_latitude_to_poles( source_data.latitude, extended_latitude, orig_lat_start_in_extend );

        // Extend out the mask
        #if DEBUG >= 1
        if (wRank == 0) { fprintf( stdout, "    Extending mask to poles\n" ); }
        #endif
        extend_mask_to_poles( source_data.mask, source_data, extended_latitude, orig_lat_start_in_extend );

        // Extend out all of the variable fields
        for(const auto& var_data : source_data.variables) {
            #if DEBUG >= 1
            if (wRank == 0) { fprintf( stdout, "    Extending variable %s to poles\n", var_data.first.c_str() ); }
            #endif
            extend_field_to_poles( source_data.variables[var_data.first], source_data, extended_latitude, orig_lat_start_in_extend );
        }

        // Update source_data to use the extended latitude
        source_data.latitude = extended_latitude;
        source_data.Nlat = source_data.latitude.size();
        source_data.myCounts[2] = source_data.Nlat;

    }

    // Compute the area of each 'cell' which will be necessary for integration
    source_data.compute_cell_areas();

    // Mask out the pole, if necessary (i.e. set lat = 90 to land)
    mask_out_pole( source_data.latitude, source_data.mask, source_data.Ntime, source_data.Ndepth, source_data.Nlat, source_data.Nlon );

    // Read in the seed
    // If extending to poles, then assume that the seed is already on the extended grid
    // since otherwise extending with zeros (or some constant) could be messy
    // the refine seed code includes the grid extensions
    double seed_count;
    bool single_seed;
    std::vector<double> Psi_seed, Phi_seed;
    if (seed_fname == "zero") {
        seed_count = 1.;
        single_seed = (seed_count == 1);
        Psi_seed.resize( source_data.Nlat * source_data.Nlon, 0.);
        Phi_seed.resize( source_data.Nlat * source_data.Nlon, 0.);
    } else {
        read_attr_from_file(seed_count, "seed_count", seed_fname);
        const int Nprocs_in_time  = source_data.Nprocs_in_time,
                  Nprocs_in_depth = source_data.Nprocs_in_depth;
        single_seed = (seed_count == 1);
        read_var_from_file( Psi_seed, tor_seed_name, seed_fname, NULL, NULL, NULL, Nprocs_in_time, Nprocs_in_depth, not(single_seed) );
        read_var_from_file( Phi_seed, pot_seed_name, seed_fname, NULL, NULL, NULL, Nprocs_in_time, Nprocs_in_depth, not(single_seed) );
    }

    // Apply to projection routine
    Apply_Helmholtz_Projection( output_fname, source_data, Psi_seed, Phi_seed, single_seed, 
            tolerance, max_iterations, use_area_weight, use_mask, Tikhov_Laplace );

    // Done!
    #if DEBUG >= 0
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
