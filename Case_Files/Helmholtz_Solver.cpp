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
#include <locale>
#include <iostream>

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
    //feenableexcept(FE_ALL_EXCEPT & ~FE_INEXACT);
    feenableexcept( FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW );

    // Specify the number of OpenMP threads
    //   and initialize the MPI world
    int thread_safety_provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &thread_safety_provided);
    //MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI::ERRORS_THROW_EXCEPTIONS);

    int wRank=-1, wSize=-1;
    MPI_Comm_rank( MPI_COMM_WORLD, &wRank );
    MPI_Comm_size( MPI_COMM_WORLD, &wSize );

    // Set the locale to the system/user specs
    //  [ this modifies how numbers are printed ]
    // on startup, the global locale is the "C" locale. we'll replace the C++ global 
    // locale and the "C" locale with the user-preferred locale for future wide character output
    //std::cout << L"User-preferred locale setting is " << std::locale("").name().c_str() << L'\n';
    std::locale::global(std::locale(""));
    std::cout.imbue(std::locale());

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
                        &adjacency_fname  = input.getCmdOption("--adjacency_file",     
                                                               "adjacency.nc",  
                                                               asked_help,
                                                               "Filename for the adjacency data (from the LLC_build_adjacency routine)"),
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

    const std::string &use_vels_string  = input.getCmdOption("--use_vels",   
                                                               "true", 
                                                               asked_help,
                                                               "Boolean (true/false) indicating if velocities should be included in solver (true) or just vorticity/divergence (false).");
    const bool use_vel = ( use_vels_string == "true" );

    const std::string &use_vort_div_string  = input.getCmdOption("--use_vort_div",   
                                                               "true", 
                                                               asked_help,
                                                               "Boolean (true/false) indicating if vorticity and divergence should be included in solver (true) or just velocity (false).");
    const bool use_vort_div = ( use_vort_div_string == "true" );

    const std::string &collapse_land_string  = input.getCmdOption("--collapse_land",   
                                                               "false", 
                                                               asked_help,
                                                               "Boolean (true/false) indicating if contiguous land should be collapse to enforce identically zero velocity.");
    const bool collapse_land = ( collapse_land_string == "true" );

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
                        &dArea_field_var_name   = input.getCmdOption("--dArea_field",   "dA",    asked_help, "Name of cell areas in input file."),
                        &tor_seed_name     = input.getCmdOption("--tor_seed",    
                                                                "Psi_seed",    
                                                                asked_help,
                                                                "Name of streamfunction (Psi) seed variable in seed file (if applicable)"),
                        &pot_seed_name     = input.getCmdOption("--pot_seed",    
                                                                "Phi_seed",    
                                                                asked_help,
                                                                "Name of potential function (Phi) seed variable in seed file (if applicable)");

    const std::string &CFL_string = input.getCmdOption("--CFL", 
                                                       "1e-6", 
                                                       asked_help,
                                                       "CFL time-step factor.");
    const double CFL = stod(CFL_string);  

    const std::string &hyper_visc_string = input.getCmdOption("--hyper_viscosity", 
                                                       "0", 
                                                       asked_help,
                                                       "Hyperviscosity coefficient [0,1].");
    const double hyper_visc = stod(hyper_visc_string);  

    const std::string &tolerance_string = input.getCmdOption("--tolerance", 
                                                             "1e-10", 
                                                             asked_help,
                                                             "Termination tolerance in terms of relative two norm.");
    const double tolerance = stod(tolerance_string);  

    const std::string &iteration_string = input.getCmdOption("--max_iterations", 
                                                             "1e5", 
                                                             asked_help,
                                                             "Maximum number of iterations for the solver before terminating. Can use exponential notation (e.g. '5e3')");
    const int max_iterations = stod(iteration_string);  

    const std::string &iter_cycle_string = input.getCmdOption("--iterations_per_cycle", 
                                                             "1e3", 
                                                             asked_help,
                                                             "Number of iterations between error calculations. Can use exponential notation (e.g. '5e3')");
    const int iters_per_cycle = stod(iter_cycle_string);  

    const std::string &num_refine_string = input.getCmdOption("--num_refinements", 
                                                             "-1", 
                                                             asked_help,
                                                             "Number of grid refinements (halvings) to use in solving");
    int num_refinements = stod(num_refine_string);  

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
    source_data.load_time(  time_dim_name,  input_fname );
    source_data.load_depth( depth_dim_name, input_fname );
    if ( constants::GRID_TYPE == constants::GridType::MeshGrid ) {
        source_data.load_latitude(  latitude_dim_name,  input_fname );
        source_data.load_longitude( longitude_dim_name, input_fname );
    } else {
        read_LLC_latlon_from_file( source_data.latitude,  latitude_dim_name,  input_fname );
        read_LLC_latlon_from_file( source_data.longitude, longitude_dim_name, input_fname );
    }
     
    // Convert to radians, if appropriate
    if ( (latlon_in_degrees == "true") and (not(constants::CARTESIAN)) ) {
        convert_coordinates( source_data.longitude, source_data.latitude );
    }

    // Apply some cleaning to the processor allotments if necessary. 
    if ( not(constants::GRID_TYPE == constants::GridType::MeshGrid ) ) {
        source_data.Nlat = 1;
        source_data.Nlon = 1;
    }
    source_data.check_processor_divisions( Nprocs_in_time_input, Nprocs_in_depth_input );

    // Read in the velocity fields
    source_data.load_variable( "u_lon", zonal_vel_name, input_fname, true, true );
    source_data.load_variable( "u_lat", merid_vel_name, input_fname, true, true );

    // Get the MPI-local dimension sizes
    source_data.Ntime  = source_data.myCounts[0];
    source_data.Ndepth = source_data.myCounts[1];

    // If we're using FILTER_OVER_LAND, then the mask has been wiped out. Load in a mask that still includes land references
    //      so that we have both. Will be used to get 'water-only' region areas.
    if (constants::FILTER_OVER_LAND) { 
        read_mask_from_file( source_data.reference_mask, zonal_vel_name, input_fname,
               source_data.Nprocs_in_time, source_data.Nprocs_in_depth );
    }

    if ( constants::GRID_TYPE == constants::GridType::MeshGrid ) {
        source_data.compute_cell_areas();
    } else {
        // Compute the area of each 'cell' which will be necessary for integration
        // Cell areas are trickier, so they will be passed in as an input.
        read_LLC_latlon_from_file( source_data.areas, dArea_field_var_name, adjacency_fname );

        // Build the adjacency matrix and other adjacency-adjacent arrays
        //  down the road, just load in a pre-built one, but for right now
        //  this is easier.
        // Need to have area loaded already!
        source_data.load_adjacency( adjacency_fname );
    }


    if ( ( constants::EXTEND_DOMAIN_TO_POLES ) and 
         ( constants::GRID_TYPE == constants::GridType::MeshGrid ) 
    ) {

        // Extend the latitude grid to reach the poles and update source_data with the new info.
        std::vector<double> extended_latitude;
        int orig_lat_start_in_extend;
        #if DEBUG >= 2
        if (wRank == 0) { fprintf( stdout, "    Extending latitude to poles\n" ); }
        #endif
        extend_latitude_to_poles( source_data.latitude, extended_latitude, orig_lat_start_in_extend );

        // Extend out the mask
        #if DEBUG >= 2
        if (wRank == 0) { fprintf( stdout, "    Extending mask to poles\n" ); }
        #endif
        extend_mask_to_poles( source_data.mask, source_data, extended_latitude, orig_lat_start_in_extend );


        extend_field_to_poles( source_data.variables.at("u_lon"), source_data, extended_latitude, orig_lat_start_in_extend );
        extend_field_to_poles( source_data.variables.at("u_lat"), source_data, extended_latitude, orig_lat_start_in_extend );

        // Extend out all of the region definitions
        source_data.latitude = extended_latitude;
        source_data.Nlat = extended_latitude.size();
        source_data.myCounts[2] = source_data.Nlat;
        source_data.compute_cell_areas();
    }
    source_data.initialize_index_functions();

    // If requested, compute (and store) sin(lat) and cos(lat)
    if ( constants::PRECOMPUTE_SIN_COS_LATITUDE ) {
        source_data.compute_sin_cos_lat();
    }

    // Number of lat/lon points
    size_t Npts;
    if ( constants::GRID_TYPE == constants::GridType::MeshGrid ) {
        Npts = source_data.areas.size();
    } else {
        Npts = source_data.latitude.size();
    }

    fprintf( stdout, "Npts (%'zu) len(uo) = %'zu, len(vo) = %'zu, len(dA) = %'zu, len(mask) = %'zu\n",
          Npts, 
          source_data.variables.at("u_lon").size(),
          source_data.variables.at("u_lat").size(),
          source_data.areas.size(),
          source_data.mask.size()
          );

    // Read in the seed
    // If extending to poles, then assume that the seed is already on the extended grid
    // since otherwise extending with zeros (or some constant) could be messy
    // the refine seed code includes the grid extensions
    double seed_count;
    bool single_seed;
    std::vector<double> Psi_seed, Phi_seed;
    bool seed_provided = false;
    if (seed_fname == "zero") {
        seed_count = 1.;
        single_seed = (seed_count == 1);
        Psi_seed.resize( Npts, 0.);
        Phi_seed.resize( Npts, 0.);

        seed_provided = false;
        fprintf( stdout, "No seed provided, will initialize with zeros.\n" );
    } else {
        //read_attr_from_file(seed_count, "seed_count", seed_fname);
        //const int Nprocs_in_time  = source_data.Nprocs_in_time,
        //          Nprocs_in_depth = source_data.Nprocs_in_depth;
        //single_seed = (seed_count == 1);
        //read_var_from_file( Psi_seed, tor_seed_name, seed_fname, NULL, NULL, NULL, Nprocs_in_time, Nprocs_in_depth, not(single_seed) );
        //read_var_from_file( Phi_seed, pot_seed_name, seed_fname, NULL, NULL, NULL, Nprocs_in_time, Nprocs_in_depth, not(single_seed) );

        read_var_from_file( Psi_seed, tor_seed_name, seed_fname, NULL, NULL, NULL );
        read_var_from_file( Phi_seed, pot_seed_name, seed_fname, NULL, NULL, NULL );

        seed_provided = true;
        fprintf( stdout, "Using provided seed.\n" );
    }

    // Auto-determine number of refinements if not set
    //  aims for ~1-degree at coarsest
    if ( num_refinements < 0 ) {
        int counter = 0;
        if ( constants::GRID_TYPE == constants::GridType::MeshGrid ) {
            while ( (180*360) * pow( 4, 2*counter ) < Npts ) { counter++; };
        } else {
            while ( 4000 * pow( 4, 2*counter ) < Npts ) { counter++; };
            num_refinements = counter+1;
        }
        #if DEBUG >= 0
        fprintf( stdout, "Setting the number of refinement levels to %d.\n", num_refinements );
        #endif
    }

    // Apply to projection routine
    if ( constants::GRID_TYPE == constants::GridType::MeshGrid ) {
        Helmholtz_Solver_MeshGrid( "output.nc", source_data, tolerance, max_iterations, iters_per_cycle,
                num_refinements, seed_provided, Psi_seed, Phi_seed );
    } else {
        Helmholtz_Solver_NonMeshGrid( "output.nc", source_data, tolerance, max_iterations, iters_per_cycle, use_area_weight,
                use_mask, use_vel, use_vort_div, collapse_land, num_refinements, CFL, 
                seed_provided, Psi_seed, Phi_seed );
    }


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
