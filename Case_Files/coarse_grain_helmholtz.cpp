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
                        &vel_input_fname   = input.getCmdOption("--velocity_input_file",        
                                                                "vels.nc",                      
                                                                asked_help,
                                                                "netCDF file containing velocity field (only used to get land mask)"),
                        &u_r_input_fname   = input.getCmdOption("--radial_velocity_input_file", 
                                                                "NONE",                         
                                                                asked_help,
                                                                "netCDF file containing radial/vertical velocity, if applicable. \nNONE if no radial velocity."),
                        &wind_input_fname  = constants::COMP_WIND_FORCE ?
                                                input.getCmdOption("--wind_tau_input_file",        
                                                                   "wind_tau_projection.nc",       
                                                                   asked_help,
                                                                   "netCDF file containing surface wind stress Helmholtz Psi and Phi")
                                                : "",
                        &quad_input_fname  = constants::COMP_PI_HELMHOLTZ ? 
                                                input.getCmdOption("--uiuj_Helmholtz_input_file",  
                                                                   "helmholtz_projection_uiuj.nc", 
                                                                   asked_help,
                                                                   "netCDF file containing Helmholtz components of u*u tensor")
                                                : "";

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

    const std::string &do_radial_derivs   = input.getCmdOption("--use_depth_derivs", 
                                                               "false", 
                                                               asked_help,
                                                               "Boolean (true/false) indicating if vertical derivatives should be used.");

    const std::string &depth_is_elevation = input.getCmdOption("--depth_is_elevation", 
                                                               "false", 
                                                               asked_help,
                                                               "Boolean (true/false) indicating if the depth grid is actually elevation \n(i.e. points down [true] or up [false])\nOnly impacts vertical derivatives / integrals.");

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

    const std::string   &tor_field_var_name     = input.getCmdOption("--tor_field", "Psi",   asked_help, "Name of toroidal field (streamfunction) in input file."),
                        &pot_field_var_name     = input.getCmdOption("--pot_field", "Phi",   asked_help, "Name of potential field (potential function) in input file."),
                        &vel_field_var_name     = input.getCmdOption("--vel_field", "u_lat", asked_help, "Name of a velocity field in input file (used to get land information)."),
                        &u_r_field_var_name     = input.getCmdOption("--u_r_field", "u_r",   asked_help, "Name of vertical/radial velocity field in input file (if used)."),
                        &wind_tau_Psi_var_name  = constants::COMP_WIND_FORCE ? input.getCmdOption("--wind_tau_Psi",  "wind_tau_Psi", asked_help) : "",
                        &wind_tau_Phi_var_name  = constants::COMP_WIND_FORCE ? input.getCmdOption("--wind_tau_Phi",  "wind_tau_Phi", asked_help) : "",
                        &uiuj_F_r_var_name      = constants::COMP_PI_HELMHOLTZ ? input.getCmdOption("--uiuj_F_r",      "uiuj_F_r",     asked_help) : "",
                        &uiuj_F_Phi_var_name    = constants::COMP_PI_HELMHOLTZ ? input.getCmdOption("--uiuj_F_Phi",    "uiuj_F_Phi",   asked_help) : "",
                        &uiuj_F_Psi_var_name    = constants::COMP_PI_HELMHOLTZ ? input.getCmdOption("--uiuj_F_Psi",    "uiuj_F_Psi",   asked_help) : "";

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
                                                                   "Name of the variable in the regions file that provides the region definitions."),
                        &coarse_map_grid_fname = input.getCmdOption("--coarse_map_grid",     
                                                                   "none",     
                                                                   asked_help,
                                                                   "netCDF file containing user-specified lat/lon grid for coarsened maps." );

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
    source_data.load_latitude(  latitude_dim_name,  Helm_input_fname );
    source_data.load_longitude( longitude_dim_name, Helm_input_fname );

    // Apply some cleaning to the processor allotments if necessary. 
    source_data.check_processor_divisions( Nprocs_in_time_input, Nprocs_in_depth_input );

    // Load in the coarsened grid, if applicable
    if ( not( coarse_map_grid_fname == "none" ) ) {
        source_data.prepare_for_coarsened_grids( coarse_map_grid_fname );
    }
     
    // Convert to radians, if appropriate
    if ( latlon_in_degrees == "true" ) {
        convert_coordinates( source_data.longitude, source_data.latitude );
    }

    // Turn on depth derivatives, if appropriate
    if ( do_radial_derivs == "true" ) {
        source_data.use_depth_derivatives = true;
        if ( constants::DiffOrd != 2 ) {
            if (wRank == 0) {
                fprintf( stderr, "WARNING: Depth grid is assumed to be non-uniform.\n" );
                fprintf( stderr, "         Non-uniform can currently only use 2nd order\n" );
                fprintf( stderr, "         derivatives. Compile-time setting was %dth order\n", constants::DiffOrd );
                fprintf( stderr, "         Will adjust to only use 2nd order for depth.\n" );
            }
        }
    }

    // Check if the depth grid is actually an elevation grid
    if ( depth_is_elevation == "true" ) {
        source_data.depth_is_elevation = true;
    }

    // Compute the area of each 'cell' which will be necessary for integration
    source_data.compute_cell_areas();

    // Read in the toroidal and potential fields
    source_data.load_variable( "F_potential", pot_field_var_name, Helm_input_fname, false, true );
    source_data.load_variable( "F_toroidal",  tor_field_var_name, Helm_input_fname, false, true );

    if ( u_r_input_fname == "NONE" ) {
        // If no u_r provided, just assume it's zero
        source_data.variables["u_r"] = std::vector<double>( source_data.variables["F_potential"].size(), 0. );
        source_data.compute_radial_vel = false;
    } else {
        source_data.load_variable( "u_r",  u_r_field_var_name, u_r_input_fname, false, true );
        source_data.compute_radial_vel = true;
    }

    // Read in the Helmholtz fields for uiuj
    if ( constants::COMP_PI_HELMHOLTZ ) {
        source_data.load_variable( "uiuj_F_r",      uiuj_F_r_var_name,      quad_input_fname, false, true );
        source_data.load_variable( "uiuj_F_Phi",    uiuj_F_Phi_var_name,    quad_input_fname, false, true );
        source_data.load_variable( "uiuj_F_Psi",    uiuj_F_Psi_var_name,    quad_input_fname, false, true );
    }

    if ( constants::COMP_WIND_FORCE ) {
        source_data.load_variable( "wind_tau_Psi", wind_tau_Psi_var_name, wind_input_fname, false, true );
        source_data.load_variable( "wind_tau_Phi", wind_tau_Phi_var_name, wind_input_fname, false, true );
    }

    // Get the MPI-local dimension sizes
    source_data.Ntime  = source_data.myCounts[0];
    source_data.Ndepth = source_data.myCounts[1];

    // Get mask : read in velocity to get the mask, and extend to the poles if needed
    source_data.load_variable( "sample_velocity", vel_field_var_name, vel_input_fname, true, true);

    // If we're using FILTER_OVER_LAND, then the mask has been wiped out. Load in a mask that still includes land references
    //      so that we have both. Will be used to get 'water-only' region areas.
    if (constants::FILTER_OVER_LAND) { 
        read_mask_from_file( source_data.reference_mask, vel_field_var_name, vel_input_fname,
               source_data.Nprocs_in_time, source_data.Nprocs_in_depth );
    }

    if ( constants::EXTEND_DOMAIN_TO_POLES ) {

        // Need to get the latitude grid for the mask
        mask_data.load_latitude(  latitude_dim_name,  vel_input_fname );
        mask_data.load_longitude( longitude_dim_name, vel_input_fname );
        mask_data.Ntime  = source_data.Ntime;
        mask_data.Ndepth = source_data.Ndepth;

        mask_data.mask = source_data.mask;
        mask_data.reference_mask = source_data.reference_mask;
     
        // Convert to radians, if appropriate
        if ( latlon_in_degrees == "true" ) {
            convert_coordinates( mask_data.longitude, mask_data.latitude );
        }

        // Extend the latitude grid to reach the poles and update source_data with the new info.
        std::vector<double> extended_latitude;
        int orig_lat_start_in_extend;
        #if DEBUG >= 2
        if (wRank == 0) { fprintf( stdout, "    Extending latitude to poles\n" ); }
        #endif
        extend_latitude_to_poles( mask_data.latitude, extended_latitude, orig_lat_start_in_extend );

        // Extend out the mask
        #if DEBUG >= 2
        if (wRank == 0) { fprintf( stdout, "    Extending mask to poles\n" ); }
        #endif
        extend_mask_to_poles( source_data.mask, mask_data, extended_latitude, orig_lat_start_in_extend );
        if (constants::FILTER_OVER_LAND) { 
            extend_mask_to_poles( source_data.reference_mask, mask_data, extended_latitude, orig_lat_start_in_extend, false );
        }

        // Extend out all of the region definitions
        mask_data.compute_cell_areas();
        mask_data.compute_region_areas();

        // Read in the region definitions and compute region areas
        if ( check_file_existence( region_defs_fname ) ) {
            // If the file exists, then read in from that
            mask_data.load_region_definitions( region_defs_fname, region_defs_dim_name, region_defs_var_name );
        } else {
            // Otherwise, just make a single region which is the entire domain
            mask_data.region_names.push_back("full_domain");
            mask_data.regions.insert( std::pair< std::string, std::vector<bool> >( 
                        "full_domain", std::vector<bool>( mask_data.Nlat * mask_data.Nlon, true) ) 
                    );
        }

        for(const auto& reg_data : mask_data.regions) {
            #if DEBUG >= 2
            if (wRank == 0) { fprintf( stdout, "    Extending region %s to poles\n", reg_data.first.c_str() ); }
            #endif
            extend_mask_to_poles( mask_data.regions[reg_data.first], mask_data, extended_latitude, orig_lat_start_in_extend, false );

            // After it's extended, copy it over into source data
            source_data.region_names.push_back(reg_data.first);
            source_data.regions.insert( std::pair< std::string, std::vector<bool> >( reg_data.first, reg_data.second ) );
        }

        source_data.region_names.push_back("AUTO_ALL");
        source_data.regions.insert( std::pair< std::string, std::vector<bool> >( 
                    "AUTO_ALL", std::vector<bool>( source_data.Nlat * source_data.Nlon, true) ) 
                );

    } else {

        // Read in the region definitions and compute region areas
        if ( check_file_existence( region_defs_fname ) ) {
            // If the file exists, then read in from that
            source_data.load_region_definitions( region_defs_fname, region_defs_dim_name, region_defs_var_name );
        } else {
            // Otherwise, just make a single region which is the entire domain
            source_data.region_names.push_back("full_domain");
            source_data.regions.insert( std::pair< std::string, std::vector<bool> >( 
                        "full_domain", std::vector<bool>( source_data.Nlat * source_data.Nlon, true) ) 
                    );
        }
    }
    source_data.compute_cell_areas();
    source_data.compute_region_areas();

    // Mask out the pole, if necessary (i.e. set lat = 90 to land)
    mask_out_pole( source_data.latitude, source_data.mask, source_data.Ntime, source_data.Ndepth, source_data.Nlat, source_data.Nlon );

    // If we need it, go ahead an merge the mask across processors now
    if (source_data.use_depth_derivatives and ( source_data.Nprocs_in_depth > 1 )) {
        source_data.gather_mask_across_depth( source_data.mask, source_data.mask_DEPTH );
    }

    // Now pass the data along to the filtering routines
    const double pre_filter_time = MPI_Wtime();
    filtering_helmholtz( source_data, filter_scales );
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
