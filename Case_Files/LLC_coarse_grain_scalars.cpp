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
#include <deque>

#include "../netcdf_io.hpp"
#include "../functions.hpp"
#include "../constants.hpp"
#include "../postprocess.hpp"

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
    const std::string   &input_fname      = input.getCmdOption("--input_file",      
                                                               "input.nc",                 
                                                               asked_help,
                                                               "netCDF file containing vector field to Helmholtz decompose"),
                        &adjacency_fname  = input.getCmdOption("--adjacency_file",     
                                                               "adjacency.nc",  
                                                               asked_help,
                                                               "Filename for the adjacency data (from the LLC_build_adjacency routine)"),
                        &dArea_fname      = input.getCmdOption("--area_file",     
                                                               "input.nc",  
                                                               asked_help,
                                                               "Filename for the cell areas");
    const std::string   &dArea_field_var_name   = input.getCmdOption("--dArea_field",   "dA",    asked_help, "Name of cell areas in input file.");

    const std::string   &time_dim_name      = input.getCmdOption("--time",        "time",      asked_help),
                        &depth_dim_name     = input.getCmdOption("--depth",       "depth",     asked_help),
                        &latitude_dim_name  = input.getCmdOption("--latitude",    "latitude",  asked_help),
                        &longitude_dim_name = input.getCmdOption("--longitude",   "longitude", asked_help);

    const std::string &latlon_in_degrees  = input.getCmdOption("--is_degrees",   "true", asked_help);

    const std::string &compute_PEKE_conv  = input.getCmdOption("--do_PEKE_conversion", "false", asked_help);

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
                                                                   "Name of the variable in the regions file that provides the region definitions.");

    // Also read in the fields to be filtered from commandline
    //   e.g. --variables "rho salinity p" (names must match with input netcdf file)
    std::vector< std::string > vars_to_filter;
    input.getListofStrings( vars_to_filter, "--variables", asked_help );
    const size_t Nvars = vars_to_filter.size();

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

    Timing_Records timing_records;

    // Initialize dataset class instance
    dataset source_data;

    // Read in source data / get size information
    #if DEBUG >= 1
    if (wRank == 0) { fprintf(stdout, "Reading in source data.\n\n"); }
    #endif

    // Read in the grid coordinates
    //   implicitely assume coordinates are the same between input files
    source_data.load_time(      time_dim_name,      input_fname );
    source_data.load_depth(     depth_dim_name,     input_fname );
    read_LLC_latlon_from_file( source_data.latitude,  latitude_dim_name,  input_fname );
    read_LLC_latlon_from_file( source_data.longitude, longitude_dim_name, input_fname );

    const bool one_snapshot = (     ( (time_dim_name  == "DNE") or (time_dim_name  == "DOES_NOT_EXIST") )
                                and ( (depth_dim_name == "DNE") or (depth_dim_name == "DOES_NOT_EXIST") )
                              );

    // Apply some cleaning to the processor allotments if necessary. 
    source_data.check_processor_divisions( Nprocs_in_time_input, Nprocs_in_depth_input );
     
    // Convert to radians, if appropriate
    if ( ( latlon_in_degrees == "true" ) and not( constants::CARTESIAN ) ) { 
        convert_coordinates( source_data.longitude, source_data.latitude ); 
    }

    // Build the adjacency matrix and other adjacency-adjacent arrays
    //  down the road, just load in a pre-built one, but for right now
    //  this is easier.
    source_data.load_adjacency( adjacency_fname );

    // Compute the area of each 'cell' which will be necessary for integration
    #if DEBUG >= 2
    if (wRank == 0) { fprintf( stdout, "Computing cell areas.\n" ); }
    #endif
    read_LLC_latlon_from_file( source_data.areas, dArea_field_var_name, dArea_fname );

    // Read in the scalar fields to filter
    #if DEBUG >= 1
    if (wRank == 0) { fprintf( stdout, "Reading in original fields.\n" ); }
    #endif
    for (size_t field_ind = 0; field_ind < vars_to_filter.size(); field_ind++) {
        source_data.load_variable( vars_to_filter.at(field_ind), vars_to_filter.at(field_ind), input_fname, true, true, !one_snapshot );
    }

    // Get the MPI-local dimension sizes
    source_data.Ntime  = one_snapshot ? 1 : source_data.myCounts[0];
    source_data.Ndepth = one_snapshot ? 1 : source_data.myCounts[1];

    const std::vector<int>  &myStarts = source_data.myStarts;

    //
    int LAT_lb, LAT_ub, Itime, Idepth, Ilat, Ilon;
    const int   Ntime   = source_data.Ntime,
                Ndepth  = source_data.Ndepth,
                Nlatlon = source_data.latitude.size();
    const size_t num_pts = Ntime * Ndepth * ( (size_t) Nlatlon );
    size_t Ivar, index;

    //
    #if DEBUG >= 1
    if (wRank == 0) { fprintf( stdout, "Setting up %'zu coarse fields.\n", Nvars ); fflush(stdout); }
    #endif
    std::vector< std::vector<double> > coarse_fields(Nvars);
    for (size_t field_ind = 0; field_ind < Nvars; field_ind++) {
        coarse_fields.at(field_ind).resize( num_pts );
    }

    //
    //// Set up filtering vectors
    //
    #if DEBUG >= 1
    if (wRank == 0) { fprintf( stdout, "Setting up filtering values.\n" ); fflush(stdout); }
    #endif
    std::vector<double > filter_values_doubles;
    std::vector<const std::vector<double>*> filter_fields;
    for (size_t field_ind = 0; field_ind < vars_to_filter.size(); field_ind++) {
        filter_fields.push_back( &source_data.variables.at(vars_to_filter.at(field_ind)) );
    }
    

    // Set up postprocessing fields
    std::vector<const std::vector<double>*> postprocess_fields;
    std::vector<std::string> postprocess_names;
    for ( Ivar = 0; Ivar < Nvars; Ivar++ ) {
        postprocess_fields.push_back( &coarse_fields.at(Ivar) );
        postprocess_names.push_back( vars_to_filter.at(Ivar) );
    }

    // Filter-loop variables
    size_t target_index, search_index, neighbour_index;
    const size_t num_neighbours = source_data.num_neighbours;
    double target_lat, target_lon, kern_val, dArea, local_val, kernel_normalization,
           local_lat, local_lon, local_dist;
    
    std::deque<size_t> points_to_test;
    std::vector<bool>   was_rejected(num_pts, false), 
                        was_accepted(num_pts, false), 
                        planned_for_testing(num_pts, false);
    // Want to use bitset directly, but does not accept runtime-specified size
    // so will use std::vector<bool>, which is a bitset under the hood
    // Note that each element is 1 bit, not byte (because of the 
    // specialization of vector<bool>), so that, while e.g. LLC4320 requires
    // ~1.8GB to store a single field at double precision, these bitsets 
    // only require <30MB. So, having each thread / processor keep their own
    // copy shouldn't be too onerous memory-wise

    //
    //// Apply filtering
    //
    for (size_t ell_ind = 0; ell_ind < filter_scales.size(); ell_ind++) {

        double scale = filter_scales.at(ell_ind);

        #if DEBUG >= 0
        if (wRank == 0) { fprintf( stdout, "Filter scale %'g km.\n", scale / 1e3 ); }
        #endif

        #pragma omp parallel \
        default(none) \
        shared( source_data, filter_fields, coarse_fields, scale, filter_scales, stdout ) \
        private( filter_values_doubles, Itime, Idepth, Ivar, \
                 target_index, search_index, neighbour_index, kern_val, dArea, local_val, \
                 was_rejected, was_accepted, planned_for_testing, points_to_test, \
                 target_lat, target_lon, kernel_normalization, \
                 local_lat, local_lon, local_dist ) \
        firstprivate( Nlatlon, Ndepth, Ntime, Nvars, ell_ind )
        {

            filter_values_doubles.clear();
            filter_values_doubles.resize( Nvars );

            was_rejected.resize( num_pts );
            was_accepted.resize( num_pts );
            planned_for_testing.resize( num_pts );

            #pragma omp for collapse(1) schedule(dynamic)
            for ( target_index = 0; target_index < Nlatlon; target_index++) {

                std::fill(was_rejected.begin(), was_rejected.end(), false);
                std::fill(was_accepted.begin(), was_accepted.end(), false);
                std::fill(planned_for_testing.begin(), planned_for_testing.end(), false);

                target_lat = source_data.latitude.at(  target_index );
                target_lon = source_data.longitude.at( target_index );


                kernel_normalization = 0.0;
                dArea = source_data.areas.at(target_index);
                kernel_normalization += 1. * dArea;
                for (Ivar = 0; Ivar < Nvars; Ivar++) {
                    local_val = filter_fields.at(Ivar)->at(target_index);
                    filter_values_doubles.at(Ivar) = dArea * local_val;
                    // not +=, here we re-set the values to the local value
                    // before looping to accumulate over the kernel
                }

                // Next, seed the 'points to test' with the adjacent points
                points_to_test.clear();
                for (neighbour_index = 0; neighbour_index < num_neighbours; neighbour_index++ ) {
                    points_to_test.push_back( source_data.adjacency_indices.at(target_index)[neighbour_index] );
                    planned_for_testing[ source_data.adjacency_indices.at(target_index)[neighbour_index] ] = true;
                }

                // So long as we still have points to test, keep testing!
                // This implements a breadth-first search through the adjacency matrix to build
                //      filtering kernel. If a point is within distance, add it to the kernel,
                //      and then test it's neighbours. If those are in, test their neighbours,
                //      and so on. If a point is too far away, we do not test it's neighbours.
                // This then assumes that for any point X within distance L of Y, that there
                //      is an adjacency path from Y to X strictly using points within distance
                //      L of Y.
                // Along the way, test for double inclusing to make sure that we don't test 
                //      points repeatedly etc.
                while ( points_to_test.size() > 0 ) {

                    // Pull out the most-recently-added point, and remove it from the 'to test' list,
                    // since we're testing it now.
                    // Since we're pulling out the most-recently-added, that effectively makes this a
                    // depth-first search to build the kernel.
                    search_index = points_to_test.front();
                    points_to_test.pop_front();
                    planned_for_testing[ search_index ] = false;

                    local_lat = source_data.latitude.at(  search_index );
                    local_lon = source_data.longitude.at( search_index );
                    local_dist = distance( target_lon, target_lat, local_lon, local_lat );

                    kern_val = kernel( local_dist, filter_scales.at(ell_ind) );

                    if ( ( kern_val > 1e-10 ) or ( local_dist <= 0.5 * filter_scales.at(ell_ind) )  ) {
                        was_accepted[search_index] = true;

                        // Accumulate the coarse values
                        dArea = source_data.areas.at(search_index);
                        kernel_normalization += kern_val * dArea;
                        for (Ivar = 0; Ivar < Nvars; Ivar++) {
                            local_val = filter_fields.at(Ivar)->at(search_index);
                            filter_values_doubles.at(Ivar) += kern_val * dArea * local_val;
                        }

                        // Since this point is in the kernel
                        // add its neighbours to the 'to search' list
                        for ( int ii = 0; ii < num_neighbours; ii++ ) {
                            neighbour_index = source_data.adjacency_indices.at(search_index)[ii];

                            // but first check if that neighbour has been rejected already
                            if ( was_rejected[neighbour_index] ) { continue; }

                            // then check if that neighbour is already accepted
                            if ( was_accepted[neighbour_index] ) { continue; }

                            // then check if that neighbour is already on the search list
                            if ( planned_for_testing[neighbour_index] ) { continue; }

                            // if it's not on any of those list already
                            // then add it to the 'points to test'
                            points_to_test.push_back( neighbour_index );
                            planned_for_testing[neighbour_index] = true;

                        }
                    } else {
                        // Otherwise, record this point as rejected, and move on
                        was_rejected[search_index] = true;
                    }

                } // loop through to make the kernel

                // Store the filtered values in the appropriate arrays
                for ( Ivar = 0; Ivar < Nvars; Ivar++ ) {
                    coarse_fields.at(Ivar).at(target_index) = filter_values_doubles.at(Ivar) / kernel_normalization;
                }
            }
        }

        //
        //// Create output file
        //
        if (not(constants::NO_FULL_OUTPUTS)) {
            char fname [50];
            snprintf(fname, 50, "filter_%.6gkm.nc", scale / 1e3);
            initialize_output_file( source_data, vars_to_filter, fname, scale );

            const int ndims = 3;
            size_t starts[ndims];
            if (one_snapshot) {
                starts[0] = 0;
                starts[1] = 0;
                starts[2] = size_t(myStarts.at(0));
            } else {
                starts[0] = size_t(myStarts.at(0));
                starts[1] = size_t(myStarts.at(1));
                starts[2] = size_t(myStarts.at(2));
            }
            size_t counts[ndims] = { size_t(Ntime), size_t(Ndepth), num_pts };

            for ( Ivar = 0; Ivar < Nvars; Ivar++ ) {
                write_field_to_output( coarse_fields.at(Ivar), vars_to_filter.at(Ivar), starts, counts, fname, &source_data.mask );
            }
        }


        //
        //// on-line postprocessing, if desired
        //

        /*
        if (constants::APPLY_POSTPROCESS) {
            MPI_Barrier(MPI_COMM_WORLD);

            #if DEBUG >= 1
            if (wRank == 0) { fprintf(stdout, "Beginning post-process routines\n"); }
            fflush(stdout);
            #endif

            std::vector<double> null_vector(0);
            Apply_Postprocess_Routines(
                    source_data, postprocess_fields, postprocess_names, null_vector,
                    scale, timing_records);

            #if DEBUG >= 1
            if (wRank == 0) { fprintf(stdout, "Finished post-process routines\n"); }
            fflush(stdout);
            #endif
        }
        */

    } // end ell loop


    // DONE!
    #if DEBUG >= 1
    fprintf(stdout, "Processor %d / %d waiting to finalize.\n", wRank + 1, wSize);
    #endif
    MPI_Finalize();
    return 0;

}
