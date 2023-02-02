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
    const std::string   &input_fname   = input.getCmdOption("--input_file",     "input.nc", asked_help);

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
    source_data.load_latitude(  latitude_dim_name,  input_fname );
    source_data.load_longitude( longitude_dim_name, input_fname );

    const bool one_snapshot = (     ( (time_dim_name  == "DNE") or (time_dim_name  == "DOES_NOT_EXIST") )
                                and ( (depth_dim_name == "DNE") or (depth_dim_name == "DOES_NOT_EXIST") )
                              );

    // Apply some cleaning to the processor allotments if necessary. 
    source_data.check_processor_divisions( Nprocs_in_time_input, Nprocs_in_depth_input );
     
    // Convert to radians, if appropriate
    if ( latlon_in_degrees == "true" ) { convert_coordinates( source_data.longitude, source_data.latitude ); }

    // Compute the area of each 'cell' which will be necessary for integration
    #if DEBUG >= 2
    if (wRank == 0) { fprintf( stdout, "Computing cell areas.\n" ); }
    #endif
    source_data.compute_cell_areas();

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
    const size_t Npts = source_data.Ntime * source_data.Ndepth * source_data.Nlat * source_data.Nlon;

    const std::vector<int>  &myStarts = source_data.myStarts;

    //
    int LAT_lb, LAT_ub, Itime, Idepth, Ilat, Ilon;
    const int   Ntime   = source_data.Ntime,
                Ndepth  = source_data.Ndepth,
                Nlat    = source_data.Nlat,
                Nlon    = source_data.Nlon;
    size_t Ivar, index;

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
        source_data.compute_region_areas();
    }

    //
    #if DEBUG >= 1
    if (wRank == 0) { fprintf( stdout, "Setting up coarse fields.\n" ); fflush(stdout); }
    #endif
    std::vector< std::vector<double> > coarse_fields(Nvars);
    for (size_t field_ind = 0; field_ind < vars_to_filter.size(); field_ind++) {
        coarse_fields.at(field_ind).resize( Npts );
    }

    //
    //// Set up filtering vectors
    //
    #if DEBUG >= 1
    if (wRank == 0) { fprintf( stdout, "Setting up filtering values.\n" ); fflush(stdout); }
    #endif
    std::vector<double > filter_values_doubles, local_kernel(Nlat * Nlon, 0.);
    std::vector<double*> filter_values_ptrs;
    std::vector<const std::vector<double>*> filter_fields;
    for (size_t field_ind = 0; field_ind < vars_to_filter.size(); field_ind++) {
        filter_fields.push_back( &source_data.variables.at(vars_to_filter.at(field_ind)) );
    }
    

    // Set up postprocessing fields
    std::vector<const std::vector<double>*> postprocess_fields;
    std::vector<std::string> postprocess_names;
    int rho_ind = -1, wo_ind = -1;
    for ( Ivar = 0; Ivar < Nvars; Ivar++ ) {
        postprocess_fields.push_back( &coarse_fields.at(Ivar) );
        postprocess_names.push_back( vars_to_filter.at(Ivar) );
        if ( compute_PEKE_conv == "true" ) {
            if (vars_to_filter.at(Ivar) == "rho") { rho_ind = Ivar; }
            if (vars_to_filter.at(Ivar) == "wo" ) { wo_ind = Ivar; }
        }
    }
    if ( compute_PEKE_conv == "true" ) {
        postprocess_names.push_back( "barrho_barwo" );
        std::vector<double> barrho_barwo( filter_fields[0]->size(), 0. );
        postprocess_fields.push_back( &barrho_barwo );
        assert( rho_ind >= 0 );
        assert( wo_ind >= 0 );
    }

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
        shared( source_data, filter_fields, coarse_fields, scale, stdout, \
              barrho_barwo, rho_ind, wo_ind ) \
        private( filter_values_doubles, filter_values_ptrs, \
                 Itime, Idepth, Ilat, Ilon, Ivar, index, \
                 LAT_lb, LAT_ub ) \
        firstprivate( local_kernel, Nlon, Nlat, Ndepth, Ntime, Nvars )
        {

            filter_values_doubles.clear();
            filter_values_doubles.resize( Nvars );

            filter_values_ptrs.clear();
            filter_values_ptrs.resize( Nvars );
            for ( Ivar = 0; Ivar < Nvars; Ivar++ ) { filter_values_ptrs.at(Ivar) = &(filter_values_doubles.at(Ivar)); }

            #pragma omp for collapse(1) schedule(dynamic)
            for (Ilat = 0; Ilat < Nlat; Ilat++) {
                get_lat_bounds(LAT_lb, LAT_ub, source_data.latitude,  Ilat, scale); 

                // If our longitude grid is uniform, and spans the full periodic domain,
                // then we can just compute it once and translate it at each lon index
                if ( (constants::PERIODIC_X) and (constants::UNIFORM_LON_GRID) and (constants::FULL_LON_SPAN) ) {
                    std::fill(local_kernel.begin(), local_kernel.end(), 0);
                    compute_local_kernel( local_kernel, scale, source_data, Ilat, 0, LAT_lb, LAT_ub );
                }

                for (Ilon = 0; Ilon < Nlon; Ilon++) {

                    if ( not( (constants::PERIODIC_X) and (constants::UNIFORM_LON_GRID) and (constants::FULL_LON_SPAN) ) ) {
                        // If we couldn't precompute the kernel earlier, then do it now
                        std::fill(local_kernel.begin(), local_kernel.end(), 0);
                        compute_local_kernel( local_kernel, scale, source_data, Ilat, Ilon, LAT_lb, LAT_ub );
                    }

                    for (Itime = 0; Itime < Ntime; Itime++) {
                        for (Idepth = 0; Idepth < Ndepth; Idepth++) {

                            // Convert our four-index to a one-index
                            index = Index(Itime, Idepth, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon);

                            // Apply the filter at the point
                            apply_filter_at_point(  filter_values_ptrs, filter_fields, source_data, Itime, Idepth, Ilat, Ilon, 
                                                    LAT_lb, LAT_ub, scale, std::vector<bool>(Nvars,false), local_kernel );

                            // Store the filtered values in the appropriate arrays
                            for ( Ivar = 0; Ivar < Nvars; Ivar++ ) {
                                coarse_fields.at(Ivar).at(index) = filter_values_doubles.at(Ivar);
                            }

                            barrho_barwo.at(index) = coarse_fields.at(rho_ind).at(index) * coarse_fields.at(wo_ind).at(index);
                        }
                    }
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

            for ( Ivar = 0; Ivar < Nvars; Ivar++ ) {
                write_field_to_output( coarse_fields.at(Ivar), vars_to_filter.at(Ivar), starts, counts, fname, &source_data.mask );
            }
        }


        //
        //// on-line postprocessing, if desired
        //

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

    } // end ell loop


    // DONE!
    #if DEBUG >= 1
    fprintf(stdout, "Processor %d / %d waiting to finalize.\n", wRank + 1, wSize);
    #endif
    MPI_Finalize();
    return 0;

}
