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

    // first argument is the flag, second argument is default value (for when flag is not present)
    //const std::string   &tor_input_fname   = input.getCmdOption("--toroidal_input_file",        "toroidal_projection.nc"),
    //                    &pot_input_fname   = input.getCmdOption("--potential_input_file",       "potential_projection.nc"),
    const std::string   &Helm_input_fname  = input.getCmdOption("--Helmholtz_input_file",       "Helmholtz_projection.nc"),
                        &vel_input_fname   = input.getCmdOption("--velocity_input_file",        "toroidal_projection.nc"),
                        &quad_input_fname  = input.getCmdOption("--uiuj_Helmholtz_input_file",  "helmholtz_projection_uiuj.nc");

    const std::string   &time_dim_name      = input.getCmdOption("--time",        "time"),
                        &depth_dim_name     = input.getCmdOption("--depth",       "depth"),
                        &latitude_dim_name  = input.getCmdOption("--latitude",    "latitude"),
                        &longitude_dim_name = input.getCmdOption("--longitude",   "longitude");

    const std::string &latlon_in_degrees  = input.getCmdOption("--is_degrees",   "true");

    const std::string   &Nprocs_in_time_string  = input.getCmdOption("--Nprocs_in_time",  "1"),
                        &Nprocs_in_depth_string = input.getCmdOption("--Nprocs_in_depth", "1");
    const int   Nprocs_in_time_input  = stoi(Nprocs_in_time_string),
                Nprocs_in_depth_input = stoi(Nprocs_in_depth_string);

    const std::string   &tor_field_var_name     = input.getCmdOption("--tor_field",     "Psi"),
                        &pot_field_var_name     = input.getCmdOption("--pot_field",     "Phi"),
                        &vel_field_var_name     = input.getCmdOption("--vel_field",     "u_lat"),
                        &uiuj_F_r_var_name      = input.getCmdOption("--uiuj_F_r",      "uiuj_F_r"),
                        &uiuj_F_Phi_var_name    = input.getCmdOption("--uiuj_F_Phi",    "uiuj_F_Phi"),
                        &uiuj_F_Psi_var_name    = input.getCmdOption("--uiuj_F_Psi",    "uiuj_F_Psi");

    const std::string   &region_defs_fname    = input.getCmdOption("--region_definitions_file",    "region_definitions.nc"),
                        &region_defs_dim_name = input.getCmdOption("--region_definitions_dim",     "region"),
                        &region_defs_var_name = input.getCmdOption("--region_definitions_var",     "region_definition");

    // Also read in the filter scales from the commandline
    //   e.g. --filter_scales "10.e3 150.76e3 1000e3" (units are in metres)
    std::vector<double> filter_scales;
    input.getFilterScales( filter_scales, "--filter_scales" );

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
     
    // Convert to radians, if appropriate
    if ( latlon_in_degrees == "true" ) {
        convert_coordinates( source_data.longitude, source_data.latitude );
    }

    // Compute the area of each 'cell' which will be necessary for integration
    source_data.compute_cell_areas();

    // Read in the toroidal and potential fields
    source_data.load_variable( "F_potential", pot_field_var_name, Helm_input_fname, false, true );
    source_data.load_variable( "F_toroidal",  tor_field_var_name, Helm_input_fname, false, true );

    // Read in the Helmholtz fields for uiuj
    if ( constants::COMP_PI_HELMHOLTZ ) {
        source_data.load_variable( "uiuj_F_r",      uiuj_F_r_var_name,      quad_input_fname, false, true );
        source_data.load_variable( "uiuj_F_Phi",    uiuj_F_Phi_var_name,    quad_input_fname, false, true );
        source_data.load_variable( "uiuj_F_Psi",    uiuj_F_Psi_var_name,    quad_input_fname, false, true );
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
