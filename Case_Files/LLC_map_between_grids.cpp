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
#include <random>

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
    const std::string   &source_adjacency_fname = input.getCmdOption("--source_adjacency_file",     
                                                               "source_adjacency.nc",  
                                                               asked_help,
                                                               "Filename for the source adjacency data (from the LLC_build_adjacency routine)"),
                        &target_adjacency_fname = input.getCmdOption("--target_adjacency_file",        
                                                                "target_adjacency.nc",                      
                                                                asked_help,
                                                               "Filename for the destination adjacency data (from the LLC_build_adjacency routine)"),
                        &source_fields_fname    = input.getCmdOption("--fields_file",
                                                                "fields.nc",
                                                                asked_help,
                                                                "Filename for the data fields to be mapped from the source grid to the target grid."),
                        &output_fname           = input.getCmdOption("--output_file",        
                                                                "grid_map.nc",                      
                                                                asked_help,
                                                               "Filename for the output");

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

    std::vector< std::string > vars_to_refine, vars_in_output;
    input.getListofStrings( vars_to_refine, "--input_variables", asked_help,
            "List of variable names (space-separated) that you want to refine / up-sample.\nNames must correspond to the name of the variable in the input file.\ne.g. 'rho u v w'" );
    input.getListofStrings( vars_in_output, "--output_variables", asked_help,
            "List of names (space-separated) that you want the variables to be called in the output file.\nNote that these must be in the same order!\ne.g. 'rho ulon ulat ur'");
    const int Nvars = vars_to_refine.size();


    if (asked_help) { return 0; }

    // Print processor assignments
    const int max_threads = omp_get_max_threads();
    omp_set_num_threads( max_threads );

    // Print some header info, depending on debug level
    print_header_info();

    // Initialize dataset class instance
    dataset source_data, target_data;

    // Read in source data / get size information
    #if DEBUG >= 1
    if (wRank == 0) { fprintf(stdout, "Reading in source data.\n\n"); }
    #endif

    // Read in the grid coordinates
    source_data.load_time(      time_dim_name,      source_fields_fname );
    source_data.load_depth(     depth_dim_name,     source_fields_fname );
    read_LLC_latlon_from_file( source_data.latitude,  latitude_dim_name,  source_adjacency_fname );
    read_LLC_latlon_from_file( source_data.longitude, longitude_dim_name, source_adjacency_fname );
    source_data.Nlat = 1;
    source_data.Nlon = 1;

    const unsigned long long int Nlatlon_source = source_data.latitude.size();

    target_data.load_time(      time_dim_name,      source_fields_fname );
    target_data.load_depth(     depth_dim_name,     source_fields_fname );
    read_LLC_latlon_from_file( target_data.latitude,  latitude_dim_name,  target_adjacency_fname );
    read_LLC_latlon_from_file( target_data.longitude, longitude_dim_name, target_adjacency_fname );
    target_data.Nlat = 1;
    target_data.Nlon = 1;

    const unsigned long long int Nlatlon_target = target_data.latitude.size();

    // Convert to radians, if appropriate
    if ( latlon_in_degrees == "true" ) {
        convert_coordinates( source_data.longitude, source_data.latitude );
        convert_coordinates( target_data.longitude, target_data.latitude );
    }

    // Load the adjacency matrix and other adjacency-adjacent arrays
    source_data.load_adjacency( source_adjacency_fname );
    target_data.load_adjacency( target_adjacency_fname );

    // Make a random number generator
    std::random_device random_device; // obtain a random number from hardware
    std::mt19937_64 random_generator(random_device()); // seed the generator
    std::uniform_int_distribution<> random_index(0, Nlatlon_source-1); // define the range

    //
    //// Now do the actual work
    //
    std::vector<unsigned long long int> nearest_source_to_target(Nlatlon_target);
    unsigned long long int source_index = 0, target_index = 0, best_index, II, best_II,
                  perc_incr = (unsigned long long int) (0.05 * (double) Nlatlon_target),
                  next_perc = perc_incr;
    double best_dist, local_dist;
    int thread_id;
    const double    typical_spacing_target = sqrt( 4 * M_PI / Nlatlon_target ) * constants::R_earth,
                    typical_spacing_source = sqrt( 4 * M_PI / Nlatlon_source ) * constants::R_earth;
    const double    typical_spacing = std::fmax( typical_spacing_target, typical_spacing_source );
    if (wRank == 0) {
        fprintf( stdout, "Typical spacing is: %g km.\n", typical_spacing / 1e3 );
        fflush(stdout);
    }

    #pragma omp parallel \
    default(none) \
    shared( source_data, target_data, nearest_source_to_target, stdout ) \
    private( best_index, target_index, source_index, local_dist, best_dist, best_II, II, \
             thread_id, next_perc ) \
    firstprivate( Nlatlon_target, Nlatlon_source, random_index, random_generator, perc_incr, typical_spacing )
    {

        thread_id = omp_get_thread_num();  // thread ID
        best_index = 0;
        //best_index = random_index(random_generator);
        next_perc = perc_incr;
        #pragma omp for collapse(1) schedule(guided)
        for ( target_index = 0; target_index < Nlatlon_target; target_index++ ) {

            /*
            if (thread_id == 0) {
                while ( target_index >= next_perc ) {
                    fprintf( stdout, "." );
                    next_perc += perc_incr;
                }
            }
            */

            best_dist = distance(
                    source_data.longitude.at(best_index),
                    source_data.latitude.at(best_index),
                    target_data.longitude.at(target_index),
                    target_data.latitude.at(target_index)
                    );

            /*
            for ( source_index = 1; source_index < Nlatlon_source; source_index++ ) {
                local_dist = distance(
                        source_data.longitude.at(source_index),
                        source_data.latitude.at(source_index),
                        target_data.longitude.at(target_index),
                        target_data.latitude.at(target_index)
                        );
                if (local_dist < best_dist) {
                    best_dist  = local_dist;
                    best_index = source_index;
                }
            }
            */

            while (true) {
                best_II = source_data.num_neighbours;
                for ( II = 0; II < source_data.num_neighbours; II++ ) {

                    source_index = source_data.adjacency_indices.at(best_index).at(II);

                    local_dist = distance(
                            source_data.longitude.at(source_index),
                            source_data.latitude.at(source_index),
                            target_data.longitude.at(target_index),
                            target_data.latitude.at(target_index)
                            );
                    if (local_dist < best_dist) {
                        best_dist = local_dist;
                        best_II = II;
                    }
                }

                if ( best_II < source_data.num_neighbours ) {
                    // If we found a closer point, move to it, and repeat
                    best_index = source_data.adjacency_indices.at(best_index).at(best_II);
                } else if ( best_dist > 10 * typical_spacing ) {
                    // If we didn't find a closer point in the neighbours, but the current
                    // point is still too far away, then just jump to a new random starting point
                    best_index = random_index(random_generator);
                } else {
                    // Otherwise, we've found the closest point. So break out of the loop now.
                    break;
                }
            }

            nearest_source_to_target.at(target_index) = best_index;
        }
    }
    fprintf( stdout, "\n" );



    //
    //// Now interpolate between the two grids
    ///     since we have the spline data already,
    ///     we'll just use that
    //

    initialize_output_file( target_data, vars_in_output, output_fname.c_str() );

    size_t starts[3] = {0,0,0};
    size_t counts[3] = {1,1,Nlatlon_target};

    std::vector<double> var_fine(Nlatlon_target);
    double II_val, mx, my, del_lon, del_lat, tmp_val;
    unsigned long long int II_index, JJ;
    for ( int Ivar = 0; Ivar < Nvars; Ivar++ ) {

        fprintf( stdout, "\n" );
        next_perc = perc_incr;

        source_data.load_variable( "coarse_field", vars_to_refine.at(Ivar), source_fields_fname, true, true );

        std::fill( var_fine.begin(), var_fine.end(), 0. );

        #pragma omp parallel \
        default(none) \
        shared( nearest_source_to_target, source_data, target_data, var_fine ) \
        private( target_index, source_index, II, mx, my, II_index, II_val, del_lon, del_lat, \
                 thread_id, next_perc ) \
        firstprivate( Nlatlon_target, perc_incr, stdout )
        {

            thread_id = omp_get_thread_num();  // thread ID
            #pragma omp for collapse(1) schedule(static)
            for ( target_index = 0; target_index < Nlatlon_target; target_index++ ) {

                /*
                if (thread_id == 0) {
                    while ( target_index >= next_perc ) {
                        fprintf( stdout, "." );
                        next_perc += perc_incr;
                    }
                }
                */

                source_index = nearest_source_to_target.at(target_index);

                // just nearest-neighbour
                //var_fine.at(target_index) = source_data.variables.at("coarse_field").at(source_index);

                mx = 0;
                my = 0;

                for ( II = 0; II < source_data.num_neighbours + 1; II++ ) {
                    II_index = source_data.adjacency_indices.at(source_index).at(II);

                    II_val = source_data.variables.at("coarse_field").at(II_index);

                    mx += source_data.adjacency_ddlon_weights.at(source_index).at(II) * II_val;
                    my += source_data.adjacency_ddlat_weights.at(source_index).at(II) * II_val;
                }

                // Right now, only uses the linear parts of the spline
                del_lon = target_data.longitude.at(target_index) - source_data.longitude.at(source_index);
                del_lat = target_data.latitude.at(target_index) - source_data.latitude.at(source_index);
                if (del_lon >   M_PI) { del_lon = 2 * M_PI - del_lon; }
                if (del_lon < - M_PI) { del_lon = 2 * M_PI + del_lon; }

                var_fine.at(target_index) = source_data.variables.at("coarse_field").at(source_index)  +  mx * del_lon  +  my * del_lat;

            }
        }

        write_field_to_output( var_fine, vars_in_output.at(Ivar), starts, counts, output_fname, NULL );
    }


    // Done!
    #if DEBUG >= 0
    const double delta_clock = MPI_Wtick();
    if (wRank == 0) {
        fprintf(stdout, "\n\n");
        fprintf(stdout, "Process completed.\n");
    }
    #endif

    #if DEBUG >= 1
    fprintf(stdout, "Processor %d / %d waiting to finalize.\n", wRank + 1, wSize);
    #endif
    MPI_Finalize();
    return 0;
}
