#include <fenv.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <algorithm>
#include <utility>
#include <math.h>
#include <vector>
#include <mpi.h>
#include <omp.h>
#include <cassert>

#include "../netcdf_io.hpp"
#include "../functions.hpp"
#include "../constants.hpp"
#include "../postprocess.hpp"
#include "../preprocess.hpp"

#include "../ALGLIB/integration.h"

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

    // Cannot have OkuboWeiss and postprocessing turned on
    static_assert( not( (constants::DO_OKUBOWEISS_ANALYSIS) and (constants::APPLY_POSTPROCESS) ),
           "" );

    // Enable all floating point exceptions but FE_INEXACT
    //feenableexcept( FE_ALL_EXCEPT & ~FE_INEXACT & ~FE_UNDERFLOW );
    //fprintf( stdout, " %d : %d \n", FE_ALL_EXCEPT, FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW | FE_INEXACT | FE_UNDERFLOW );
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
        fprintf( stdout, "The command-line input arguments [and default values] are:\n" );
    }

    // first argument is the flag, second argument is default value (for when flag is not present)
    const std::string   &input_fname   = input.getCmdOption("--input_file",     "input.nc", asked_help),
                        &vel_input_fname   = input.getCmdOption("--velocity_input_file",        
                                                                "vels.nc",                      
                                                                asked_help,
                                                                "netCDF file containing velocity field (only used to get land mask)");

    const std::string   &time_dim_name      = input.getCmdOption("--time",        "time",      asked_help),
                        &depth_dim_name     = input.getCmdOption("--depth",       "depth",     asked_help),
                        &latitude_dim_name  = input.getCmdOption("--latitude",    "latitude",  asked_help),
                        &longitude_dim_name = input.getCmdOption("--longitude",   "longitude", asked_help);

    const std::string &latlon_in_degrees  = input.getCmdOption("--is_degrees",   "true", asked_help);

    const std::string   &Nprocs_in_time_string  = input.getCmdOption("--Nprocs_in_time",  "1", asked_help),
                        &Nprocs_in_depth_string = input.getCmdOption("--Nprocs_in_depth", "1", asked_help);
    const int   Nprocs_in_time_input  = stoi(Nprocs_in_time_string),
                Nprocs_in_depth_input = stoi(Nprocs_in_depth_string);

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

    const std::string   &tor_field_var_name     = input.getCmdOption("--tor_field", "Psi",   asked_help, "Name of toroidal field (streamfunction) in input file."),
                        &pot_field_var_name     = input.getCmdOption("--pot_field", "Phi",   asked_help, "Name of potential field (potential function) in input file."),
                        &vel_field_var_name     = input.getCmdOption("--vel_field", "u_lat", asked_help, "Name of a velocity field in input file (used to get land information).");

    // Also read in the filter scale from the commandline
    //  this script only accepts a *single* filter scale
    const std::string   &filter_scale_string  = input.getCmdOption("--filter_scale",  "1e4", asked_help);
    const double   filter_scale  = stod(filter_scale_string);

    const std::string   &num_integration_steps_string  = input.getCmdOption("--integration_steps",  "2", asked_help);
    const int num_integration_steps = stoi(num_integration_steps_string);

    if (asked_help) { return 0; }

    // Print processor assignments
    const int max_threads = omp_get_max_threads();
    omp_set_num_threads( max_threads );

    // Print some header info, depending on debug level
    print_header_info();

    Timing_Records timing_records;

    // Initialize dataset class instance
    dataset source_data, mask_data;

    // Read in source data / get size information
    #if DEBUG >= 1
    if (wRank == 0) { fprintf(stdout, "Reading in source data.\n\n"); }
    #endif

    // Read in the grid coordinates
    //   implicitely assume coordinates are the same between input files
    source_data.load_time(      time_dim_name,      input_fname );
    source_data.load_depth(     depth_dim_name,     input_fname );
    source_data.load_latitude(  latitude_dim_name,  input_fname );
    source_data.load_longitude( longitude_dim_name, input_fname );

    const std::vector<double>   &latitude   = source_data.latitude,
                                &longitude  = source_data.longitude;

    const bool one_snapshot = (     ( (time_dim_name  == "DNE") or (time_dim_name  == "DOES_NOT_EXIST") )
                                and ( (depth_dim_name == "DNE") or (depth_dim_name == "DOES_NOT_EXIST") )
                              );

    // Apply some cleaning to the processor allotments if necessary. 
    source_data.check_processor_divisions( Nprocs_in_time_input, Nprocs_in_depth_input );

    // Load in the coarsened grid, if applicable
    if ( not( coarse_map_grid_fname == "none" ) ) {
        source_data.prepare_for_coarsened_grids( coarse_map_grid_fname );
    }
     
    // Convert to radians, if appropriate
    if ( latlon_in_degrees == "true" ) { convert_coordinates( source_data.longitude, source_data.latitude ); }

    // Compute the area of each 'cell' which will be necessary for integration
    #if DEBUG >= 2
    if (wRank == 0) { fprintf( stdout, "Computing cell areas.\n" ); }
    #endif
    source_data.compute_cell_areas();

    // Read in the toroidal and potential fields
    source_data.load_variable( "F_potential", pot_field_var_name, input_fname, false, true );
    source_data.load_variable( "F_toroidal",  tor_field_var_name, input_fname, false, true );

    // Get the MPI-local dimension sizes
    source_data.Ntime  = one_snapshot ? 1 : source_data.myCounts[0];
    source_data.Ndepth = one_snapshot ? 1 : source_data.myCounts[1];
    const size_t Npts = source_data.Ntime * source_data.Ndepth * source_data.Nlat * source_data.Nlon;
    const std::vector<int>  &myStarts = source_data.myStarts;

    //
    int LAT_lb, LAT_ub, Itime, Idepth, Ilat, Ilon;
    const int   Ntime   = source_data.Ntime,
                Ndepth  = source_data.Ndepth,
                Nlat    = source_data.Nlat,
                Nlon    = source_data.Nlon;
    size_t Ivar, index;

    // Get mask : read in velocity to get the mask, and extend to the poles if needed
    source_data.load_variable( "sample_velocity", vel_field_var_name, vel_input_fname, true, true);
    const std::vector<bool> &mask = source_data.mask;

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
        #if DEBUG >= 2
        if (wRank == 0) { fprintf( stdout, "    Extending region definitions\n" ); }
        #endif
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


    //
    std::vector< double > 
            strain_energy(Npts, 0.), KE_strain(Npts,0.),
            cyclonic_energy(Npts, 0.), KE_cyclonic(Npts,0.),
            anticyclonic_energy(Npts, 0.), KE_anticyclonic(Npts,0.),
            prev_strain(Npts, 0.), current_strain(Npts,0.),
            prev_cyclonic(Npts, 0.), current_cyclonic(Npts,0.),
            prev_anticyclonic(Npts, 0.), current_anticyclonic(Npts,0.),
            u_lon_tor(Npts,0), u_lon_pot(Npts, 0), u_lon_tot(Npts, 0),
            u_lat_tor(Npts,0), u_lat_pot(Npts, 0), u_lat_tot(Npts, 0);
    std::vector<double> null_vector(0);

    //
    //// Set up filtering vectors
    //
    #if DEBUG >= 1
    if (wRank == 0) { fprintf( stdout, "Setting up filtering values.\n" ); fflush(stdout); }
    #endif
    std::vector<double > filter_values_doubles, local_kernel(Nlat * Nlon, 0.);
    std::vector<double*> filter_values_ptrs, null_ptrs_vector;
    std::vector<double> coarse_Psi( Npts, 0. ),
                        coarse_Phi( Npts, 0. );
    std::vector<const std::vector<double>*> filter_fields;
    filter_fields.push_back( &source_data.variables.at("F_potential") );
    filter_fields.push_back( &source_data.variables.at("F_toroidal") );
    double dl_kern, dll_kern;
    
    #if DEBUG >= 0
    if (wRank == 0) { fprintf( stdout, "Will filter from 0 km to %'g km, with %'d steps\n",
                        filter_scale / 1e3, num_integration_steps ); }
    #endif

    //
    //// Get the Gauss-Legendre quadrature weights and nodes
    //
    alglib::real_1d_array weights_array, nodes_array;
    std::vector<double> quad_weights(num_integration_steps, 0),
                        quad_nodes(num_integration_steps, 0);
    weights_array.attach_to_ptr( num_integration_steps, &quad_weights[0] );
    nodes_array.attach_to_ptr( num_integration_steps, &quad_nodes[0] );

    alglib::ae_int_t info;
    //alglib::gqgenerategausslegendre( num_integration_steps, info, nodes_array, weights_array );
    alglib::gqgenerategaussjacobi( num_integration_steps, 0, 1, info, nodes_array, weights_array );
    for (size_t II = 0; II < num_integration_steps; ++II ) {
        // Adjust the weights and nodes to [0,ell] instead of [-1,1]
        nodes_array[II]   = (filter_scale/2.) * (1 + nodes_array[II]);
        weights_array[II] = (filter_scale/2.) * weights_array[II];
        //if (wRank == 0) { fprintf( stdout, "%g -> %g\n", nodes_array[II], weights_array[II] ); }
    }
    

    //
    //// Apply filtering
    //
    double prev_scale = 0., scale_delta, scale_l2_d2;
    for (size_t ell_ind = 0; ell_ind < num_integration_steps; ell_ind++) {

        //scale_delta = filter_scale * ((double)ell_ind) / (num_integration_steps - 1.0);
        scale_delta = nodes_array[ell_ind];

        #if DEBUG >= 0
        if (wRank == 0) { fprintf( stdout, "First filter scale %'g km.\n", scale_delta / 1e3 ); }
        #endif

        filter_fields.clear();
        filter_fields.push_back( &source_data.variables.at("F_potential") );
        filter_fields.push_back( &source_data.variables.at("F_toroidal") );

        #pragma omp parallel \
        default(none) \
        shared( source_data, filter_fields, coarse_Psi, coarse_Phi, null_vector ) \
        private( filter_values_doubles, filter_values_ptrs, null_ptrs_vector, \
                 Itime, Idepth, Ilat, Ilon, Ivar, index, \
                 LAT_lb, LAT_ub, dl_kern, dll_kern ) \
        firstprivate( local_kernel, Nlon, Nlat, Ndepth, Ntime, scale_delta )
        {

            filter_values_doubles.clear();
            filter_values_doubles.resize( 2 );

            null_ptrs_vector.clear();
            filter_values_ptrs.clear();
            filter_values_ptrs.resize( 2 );
            for ( Ivar = 0; Ivar < 2; Ivar++ ) { filter_values_ptrs.at(Ivar) = &(filter_values_doubles.at(Ivar)); }

            #pragma omp for collapse(1) schedule(dynamic)
            for (Ilat = 0; Ilat < Nlat; Ilat++) {
                get_lat_bounds(LAT_lb, LAT_ub, source_data.latitude,  Ilat, scale_delta); 

                // If our longitude grid is uniform, and spans the full periodic domain,
                // then we can just compute it once and translate it at each lon index
                if ( (constants::PERIODIC_X) and (constants::UNIFORM_LON_GRID) and (constants::FULL_LON_SPAN) ) {
                    std::fill(local_kernel.begin(), local_kernel.end(), 0);
                    compute_local_kernel( local_kernel, null_vector, null_vector, scale_delta, source_data, 
                            Ilat, 0, LAT_lb, LAT_ub );
                }

                for (Ilon = 0; Ilon < Nlon; Ilon++) {

                    if ( not( (constants::PERIODIC_X) and (constants::UNIFORM_LON_GRID) and (constants::FULL_LON_SPAN) ) ) {
                        // If we couldn't precompute the kernel earlier, then do it now
                        std::fill(local_kernel.begin(), local_kernel.end(), 0);
                        compute_local_kernel( local_kernel, null_vector, null_vector, scale_delta, source_data, 
                                Ilat, Ilon, LAT_lb, LAT_ub );
                    }

                    for (Itime = 0; Itime < Ntime; Itime++) {
                        for (Idepth = 0; Idepth < Ndepth; Idepth++) {

                            // Convert our four-index to a one-index
                            index = Index(Itime, Idepth, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon);

                            if ( not(constants::FILTER_OVER_LAND) and not(source_data.mask.at(index)) ) {
                                coarse_Phi.at(index) = constants::fill_value;
                                coarse_Psi.at(index) = constants::fill_value;
                            } else{
                                // Apply the filter at the point
                                apply_filter_at_point(  filter_values_ptrs, null_ptrs_vector, null_ptrs_vector, 
                                                        dl_kern, dll_kern,
                                                        filter_fields, source_data, Itime, Idepth, Ilat, Ilon, 
                                                        LAT_lb, LAT_ub, scale_delta, std::vector<bool>(2,false), 
                                                        local_kernel, null_vector, null_vector );

                                coarse_Phi.at(index) = filter_values_doubles[0];
                                coarse_Psi.at(index) = filter_values_doubles[1];
                            }
                        }
                    }
                }
            }
        }

        // Get pot and tor velocities
        toroidal_vel_from_F( u_lon_tor, u_lat_tor, coarse_Psi, longitude, latitude, Ntime, Ndepth, Nlat, Nlon, mask);
        potential_vel_from_F(u_lon_pot, u_lat_pot, coarse_Phi, longitude, latitude, Ntime, Ndepth, Nlat, Nlon, mask);

        #pragma omp parallel \
        default( none ) \
        shared( mask, u_lon_tor, u_lat_tor, u_lon_pot, u_lat_pot, u_lon_tot, u_lat_tot ) \
        private( index )
        {
            #pragma omp for collapse(1) schedule(static)
            for (index = 0; index < u_lon_tor.size(); ++index) {
                if ( mask.at(index) ) {
                    u_lon_tot.at(index) = u_lon_tor.at(index) + u_lon_pot.at(index);
                    u_lat_tot.at(index) = u_lat_tor.at(index) + u_lat_pot.at(index);
                }
            }
        }

        // Now that we have filtered Psi and Phi, we need to get the scalars |S|^2 and |Omega|^2
        std::vector<double> zero_array( Npts, 0 );
        compute_vorticity(
                null_vector, null_vector, null_vector, null_vector, null_vector, 
                KE_cyclonic, KE_anticyclonic, KE_strain,
                source_data, zero_array, u_lon_tot, u_lat_tot);
    
        // And now we need to filter those scalar fields, using scale
        scale_l2_d2 = sqrt( pow(filter_scale,2) - pow(scale_delta,2) );

        #if DEBUG >= 0
        if (wRank == 0) { fprintf( stdout, "secondary filtering at scale %'g km.\n", scale_l2_d2 / 1e3 ); }
        #endif

        filter_fields.clear();
        filter_fields.push_back( &KE_strain );
        filter_fields.push_back( &KE_cyclonic );
        filter_fields.push_back( &KE_anticyclonic );

        #pragma omp parallel \
        default(none) \
        shared( source_data, filter_fields, stdout, current_strain, KE_strain, \
                current_cyclonic, KE_cyclonic, current_anticyclonic, KE_anticyclonic, \
                null_vector \
              ) \
        private( filter_values_doubles, filter_values_ptrs, null_ptrs_vector, \
                 Itime, Idepth, Ilat, Ilon, Ivar, index, \
                 LAT_lb, LAT_ub, dl_kern, dll_kern ) \
        firstprivate( local_kernel, Nlon, Nlat, Ndepth, Ntime, scale_l2_d2 )
        {

            filter_values_doubles.clear();
            filter_values_doubles.resize( 3 );

            null_ptrs_vector.clear();
            filter_values_ptrs.clear();
            filter_values_ptrs.resize( 3 );
            for ( Ivar = 0; Ivar < 3; Ivar++ ) { filter_values_ptrs.at(Ivar) = &(filter_values_doubles.at(Ivar)); }

            #pragma omp for collapse(1) schedule(dynamic)
            for (Ilat = 0; Ilat < Nlat; Ilat++) {
                get_lat_bounds(LAT_lb, LAT_ub, source_data.latitude,  Ilat, scale_l2_d2); 

                // If our longitude grid is uniform, and spans the full periodic domain,
                // then we can just compute it once and translate it at each lon index
                if ( (constants::PERIODIC_X) and (constants::UNIFORM_LON_GRID) and (constants::FULL_LON_SPAN) ) {
                    std::fill(local_kernel.begin(), local_kernel.end(), 0);
                    compute_local_kernel( local_kernel, null_vector, null_vector, scale_l2_d2, source_data, 
                                            Ilat, 0, LAT_lb, LAT_ub );
                }

                for (Ilon = 0; Ilon < Nlon; Ilon++) {

                    if ( not( (constants::PERIODIC_X) and (constants::UNIFORM_LON_GRID) and (constants::FULL_LON_SPAN) ) ) {
                        // If we couldn't precompute the kernel earlier, then do it now
                        std::fill(local_kernel.begin(), local_kernel.end(), 0);
                        compute_local_kernel( local_kernel, null_vector, null_vector, scale_l2_d2, source_data, 
                                                Ilat, Ilon, LAT_lb, LAT_ub );
                    }

                    for (Itime = 0; Itime < Ntime; Itime++) {
                        for (Idepth = 0; Idepth < Ndepth; Idepth++) {

                            // Convert our four-index to a one-index
                            index = Index(Itime, Idepth, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon);

                            if ( not(constants::FILTER_OVER_LAND) and not(source_data.mask.at(index)) ) {
                                current_strain.at(index) = constants::fill_value;
                                current_cyclonic.at(index) = constants::fill_value;
                                current_anticyclonic.at(index) = constants::fill_value;
                            } else{
                                // Apply the filter at the point
                                apply_filter_at_point(  filter_values_ptrs, null_ptrs_vector, null_ptrs_vector, 
                                                        dl_kern, dll_kern,
                                                        filter_fields, source_data, Itime, Idepth, Ilat, Ilon, 
                                                        LAT_lb, LAT_ub, scale_l2_d2, std::vector<bool>(3,false), 
                                                        local_kernel, null_vector, null_vector );

                                current_strain.at(index) = filter_values_doubles[0];
                                current_cyclonic.at(index) = filter_values_doubles[1];
                                current_anticyclonic.at(index) = filter_values_doubles[2];
                            }
                        }
                    }
                }
            }
        }


        // Now that we have computed the double-filtered fields, accumulate the
        // integral, using trapezoid rule
        //if ( ell_ind > 0 ) {
        // Gauss-Legendre
        //double current_weight = 2 * scale_delta * weights_array[ell_ind];

        // Gauss-Jacobi
        double current_weight = 2 * (filter_scale / 2.) * weights_array[ell_ind];

        // Trapezoidal
        //double current_weight = 2 * ( scale_delta - prev_scale ); 
        if ( true ) {
            #pragma omp parallel \
            default( none ) \
            shared( mask, strain_energy, current_strain, prev_strain, \
                    cyclonic_energy, current_cyclonic, prev_cyclonic, \
                    anticyclonic_energy, current_anticyclonic, prev_anticyclonic ) \
            private( index ) \
            firstprivate( current_weight, scale_delta, prev_scale )
            {
                #pragma omp for collapse(1) schedule(static)
                for (index = 0; index < strain_energy.size(); ++index) {
                    if ( mask.at(index) ) {
                        strain_energy.at(index)       += current_weight * current_strain.at(index);
                        cyclonic_energy.at(index)     += current_weight * current_cyclonic.at(index);
                        anticyclonic_energy.at(index) += current_weight * current_anticyclonic.at(index);
                        /*
                        strain_energy.at(index)       += current_weight * 
                                                         0.5 * (  scale_delta * current_strain.at(index) 
                                                                + prev_scale * prev_strain.at(index) );
                        cyclonic_energy.at(index)     += current_weight * 
                                                         0.5 * (  scale_delta * current_cyclonic.at(index) 
                                                                + prev_scale * prev_cyclonic.at(index) );
                        anticyclonic_energy.at(index) += current_weight * 
                                                         0.5 * (  scale_delta * current_anticyclonic.at(index) 
                                                                + prev_scale * prev_anticyclonic.at(index) );
                        */
                    }
                }
            }
        }
        /*
        std::swap( current_strain, prev_strain );
        std::swap( current_cyclonic, prev_cyclonic );
        std::swap( current_anticyclonic, prev_anticyclonic );
        */
        prev_scale = scale_delta;


    } // end ell loop

    
    // Create output file
    std::vector<std::string> output_names;
    output_names.push_back( "strain_KE" );
    output_names.push_back( "cyclonic_KE" );
    output_names.push_back( "anticyclonic_KE" );

    if (not(constants::NO_FULL_OUTPUTS)) {
        char fname [50];
        snprintf(fname, 50, "filter_%.6gkm.nc", filter_scale / 1e3);
        initialize_output_file( source_data, output_names, fname, filter_scale );

        const int ndims = 4;
        size_t starts[ndims];
        if (one_snapshot) {
            starts[0] = 0;
            starts[1] = 0;
            starts[2] = size_t(myStarts.at(0));
            starts[3] = size_t(myStarts.at(1));
        } else {
       starts[0] = size_t(myStarts.at(0));
       starts[1] = size_t(myStarts.at(1));
       starts[2] = size_t(myStarts.at(2));
       starts[3] = size_t(myStarts.at(3));
        }
        size_t counts[ndims] = { size_t(Ntime),          size_t(Ndepth),         size_t(Nlat),           size_t(Nlon)           };

        write_field_to_output( strain_energy, "strain_KE", starts, counts, fname, &source_data.mask );
        write_field_to_output( cyclonic_energy, "cyclonic_KE", starts, counts, fname, &source_data.mask );
        write_field_to_output( anticyclonic_energy, "anticyclonic_KE", starts, counts, fname, &source_data.mask );
    }


    //
    //// on-line postprocessing, if desired
    //
    if (constants::APPLY_POSTPROCESS) {

        // Set up postprocessing fields
        std::vector<const std::vector<double>*> postprocess_fields;
        std::vector<std::string> postprocess_names;

        postprocess_fields.push_back( &strain_energy );
        postprocess_fields.push_back( &cyclonic_energy );
        postprocess_fields.push_back( &anticyclonic_energy );
        postprocess_names.push_back( "strain_KE" );
        postprocess_names.push_back( "cyclonic_KE" );
        postprocess_names.push_back( "anticyclonic_KE" );

        MPI_Barrier(MPI_COMM_WORLD);

        #if DEBUG >= 1
        if (wRank == 0) { fprintf(stdout, "Beginning post-process routines\n"); }
        fflush(stdout);
        #endif

        std::vector<double> null_vector(0);
        Apply_Postprocess_Routines(
                source_data, postprocess_fields, postprocess_names, null_vector,
                filter_scale, timing_records);

        #if DEBUG >= 1
        if (wRank == 0) { fprintf(stdout, "Finished post-process routines\n"); }
        fflush(stdout);
        #endif
    }



    // DONE!
    #if DEBUG >= 1
    fprintf(stdout, "Processor %d / %d waiting to finalize.\n", wRank + 1, wSize);
    #endif
    MPI_Finalize();
    return 0;

}
