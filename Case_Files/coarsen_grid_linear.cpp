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
    const bool asked_help = input.cmdOptionExists("--help");
    if (asked_help) {
        fprintf( stdout, "\033[1;4mThe command-line input arguments [and default values] are:\033[0m\n" );
    }

    // first argument is the flag, second argument is default value (for when flag is not present)
    const std::string   &coarse_fname   = input.getCmdOption("--coarse_file",   
                                                             "coarse.nc",        
                                                             asked_help,
                                                             "netCDF file containing the grid onto which you want to downsample / coarsen."),
                        &fine_fname     = input.getCmdOption("--fine_file",     
                                                             "fine.nc",          
                                                             asked_help,
                                                             "netCDf file containing the variables that you want to downsample / coarsen."),
                        &output_fname   = input.getCmdOption("--output_file",   
                                                             "coarse_vel.nc",    
                                                             asked_help,
                                                             "Filename for where the downsampled variables will be stored.");

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

    std::vector< std::string > vars_to_refine, vars_in_output;
    input.getListofStrings( vars_to_refine, "--input_variables",  asked_help, 
            "List of variable names (space-separated) that you want to down-sample.\nNames must correspond to the name of the variable in the input file.\ne.g. 'rho u v w'" );
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
    dataset coarse_data, fine_data;

    // Read in source data / get size information
    #if DEBUG >= 1
    if (wRank == 0) { fprintf(stdout, "Reading in source data.\n\n"); }
    #endif

    // Read in the grid coordinates
    //coarse_data.load_time(      time_dim_name,      coarse_fname );
    coarse_data.load_depth(     depth_dim_name,     coarse_fname );
    coarse_data.load_latitude(  latitude_dim_name,  coarse_fname );
    coarse_data.load_longitude( longitude_dim_name, coarse_fname );

    fine_data.load_time(      time_dim_name,      fine_fname );
    coarse_data.time.resize( fine_data.time.size() );
    coarse_data.time = fine_data.time;
    coarse_data.full_Ntime = fine_data.full_Ntime;

    fine_data.load_depth(     depth_dim_name,     fine_fname );
    fine_data.load_latitude(  latitude_dim_name,  fine_fname );
    fine_data.load_longitude( longitude_dim_name, fine_fname );

    const bool  COARSE_LAT_GRID_INCREASING = ( coarse_data.latitude[1]  > coarse_data.latitude[0]  ) ? true : false,
                COARSE_LON_GRID_INCREASING = ( coarse_data.longitude[1] > coarse_data.longitude[0] ) ? true : false;

    const bool  FINE_LAT_GRID_INCREASING = ( fine_data.latitude[1]  > fine_data.latitude[0]  ) ? true : false,
                FINE_LON_GRID_INCREASING = ( fine_data.longitude[1] > fine_data.longitude[0] ) ? true : false;

    #if DEBUG >= 0
    if (wRank == 0) {
        fprintf( stdout, "Coarse latitude grid is %s and longitude grid is %s\n",
                COARSE_LAT_GRID_INCREASING ? "increasing" : "decreasing",
                COARSE_LON_GRID_INCREASING ? "increasing" : "decreasing"
               );
        fprintf( stdout, "Fine latitude grid is %s and longitude grid is %s\n",
                FINE_LAT_GRID_INCREASING ? "increasing" : "decreasing",
                FINE_LON_GRID_INCREASING ? "increasing" : "decreasing"
               );
    }
    #endif

    // Apply some cleaning to the processor allotments if necessary. 
    coarse_data.check_processor_divisions( Nprocs_in_time_input, Nprocs_in_depth_input );
    fine_data.Nprocs_in_time = coarse_data.Nprocs_in_time;
    fine_data.Nprocs_in_depth = coarse_data.Nprocs_in_depth;
     
    // Convert to radians, if appropriate
    if ( latlon_in_degrees == "true" ) {
        convert_coordinates( coarse_data.longitude, coarse_data.latitude );
        convert_coordinates( fine_data.longitude,   fine_data.latitude );
    }

    // Read in velocity fields
    fine_data.load_variable( "fine_field", vars_to_refine.at(0), fine_fname, true, true );

    const int   full_Ntime  = fine_data.full_Ntime,
                Ntime       = fine_data.myCounts[0],
                Ndepth      = fine_data.myCounts[1],
                Nlat_coarse = coarse_data.Nlat,
                Nlon_coarse = coarse_data.Nlon,
                Nlat_fine   = fine_data.Nlat,
                Nlon_fine   = fine_data.Nlon;

    size_t starts[4] = { (size_t) fine_data.myStarts.at(0), 
                         (size_t) fine_data.myStarts.at(1), 
                         0,                
                         0                
                       };
    size_t counts[4] = { (size_t) fine_data.myCounts.at(0), 
                         (size_t) fine_data.myCounts.at(1), 
                         (size_t) coarse_data.Nlat, 
                         (size_t) coarse_data.Nlon };

    #if DEBUG >= 1
    if (wRank == 0) {
        fprintf( stdout, " (%zu,%zu,%zu,%zu) -> (%zu,%zu,%zu,%zu)\n", starts[0], starts[1], starts[2], starts[3],
                                                                      counts[0], counts[1], counts[2], counts[3]);
    }
    #endif

    // Compute the area of each 'cell' which will be necessary for creating the output file
    if (wRank == 0) { fprintf( stdout, "Computing cell areas.\n" ); }
    coarse_data.compute_cell_areas();

    // Initialize file and write out coarsened fields
    if (wRank == 0) { fprintf( stdout, "Preparing output file\n" ); }
    vars_in_output.push_back("mask");
    initialize_output_file( coarse_data, vars_in_output, output_fname.c_str() );

    #if DEBUG >= 1
    if (wRank == 0) {
        fprintf( stdout, " c(%d,%d,%d,%d) -> f(%d,%d,%d,%d)\n", full_Ntime, Ndepth, Nlat_fine,   Nlon_fine,
                                                                full_Ntime, Ndepth, Nlat_coarse, Nlon_coarse );
    }
    #endif

    // Now coarsen the velocity fields
    const size_t    Npts_coarse = Ntime * Ndepth * ((size_t) Nlat_coarse) * ((size_t) Nlon_coarse);
    const size_t    Npts_fine   = Ntime * Ndepth * ((size_t) Nlat_fine)   * ((size_t) Nlon_fine);

    std::vector<double> var_coarse(Npts_coarse);
    std::vector<bool> mask_coarse(Npts_coarse, false);

    // Next, the coarse velocities
    int cnt, land_cnt, Itime, Idepth, LEFT, RIGHT, BOT, TOP, Ilat_fine, Ilon_fine, Ilat_coarse, Ilon_coarse;
    double target_lat, target_lon, interp_val;
    size_t II_fine, II_coarse, coarse_mask_count;

    for ( int Ivar = 0; Ivar < Nvars; Ivar++ ) {

        fine_data.load_variable( "fine_field", vars_to_refine.at(Ivar), fine_fname, true, true );
    
        std::fill( var_coarse.begin(), var_coarse.end(), 0. );
        std::fill( mask_coarse.begin(), mask_coarse.end(), false );

        coarse_mask_count = 0;

        #pragma omp parallel \
        default(none) \
        shared( coarse_data, fine_data, var_coarse, mask_coarse, vars_to_refine, stdout, stderr ) \
        private( target_lat, target_lon, Itime, Idepth, II_coarse, Ilat_coarse, Ilon_coarse, \
                 RIGHT, LEFT, BOT, TOP, II_fine, Ilat_fine, Ilon_fine, cnt, interp_val, land_cnt ) \
        firstprivate( Npts_coarse, Nlon_coarse, Nlat_coarse, Nlon_fine, Nlat_fine, Ntime, Ndepth, \
                      COARSE_LAT_GRID_INCREASING, COARSE_LON_GRID_INCREASING, FINE_LAT_GRID_INCREASING, FINE_LON_GRID_INCREASING ) \
        reduction( + : coarse_mask_count )
        {
            #pragma omp for collapse(1) schedule(guided)
            for (II_coarse = 0; II_coarse < Npts_coarse; ++II_coarse) {

                Index1to4( II_coarse, Itime, Idepth, Ilat_coarse, Ilon_coarse, Ntime, Ndepth, Nlat_coarse, Nlon_coarse );


                //
                //// Get bounding box indices in the fine grid.
                //

                // bottom
                if ( COARSE_LAT_GRID_INCREASING ) {
                    if ( Ilat_coarse == 0 ) {
                        target_lat = coarse_data.latitude.at(Ilat_coarse);
                    } else {
                        target_lat = 0.5 * ( coarse_data.latitude.at(Ilat_coarse) + coarse_data.latitude.at(Ilat_coarse - 1) );
                    }
                } else {
                    if ( Ilat_coarse < Nlat_coarse - 1 ) {
                        target_lat = 0.5 * ( coarse_data.latitude.at(Ilat_coarse) + coarse_data.latitude.at(Ilat_coarse + 1) );
                    } else {
                        target_lat = coarse_data.latitude.at(Ilat_coarse);
                    }
                }
                if ( FINE_LAT_GRID_INCREASING ) {
                    BOT = std::lower_bound( fine_data.latitude.begin(), fine_data.latitude.end(), target_lat ) 
                          - fine_data.latitude.begin();
                } else {
                    //BOT = std::upper_bound( fine_data.latitude.begin(), fine_data.latitude.end(), target_lat ) 
                    BOT = std::upper_bound( fine_data.latitude.rbegin(), fine_data.latitude.rend(), target_lat ) 
                          - fine_data.latitude.rbegin();
                    BOT = (Nlat_fine - 1) - BOT;
                }
                BOT = (BOT < 0) ? 0 : (BOT >= Nlat_fine) ? Nlat_fine - 1 : BOT;

                // top
                if ( COARSE_LAT_GRID_INCREASING ) {
                    if ( Ilat_coarse < Nlat_coarse - 1 ) {
                        target_lat = 0.5 * ( coarse_data.latitude.at(Ilat_coarse) + coarse_data.latitude.at(Ilat_coarse + 1) );
                    } else {
                        target_lat = coarse_data.latitude.at(Ilat_coarse);
                    }
                } else {
                    if ( Ilat_coarse == 0 ) {
                        target_lat = coarse_data.latitude.at(Ilat_coarse);
                    } else {
                        target_lat = 0.5 * ( coarse_data.latitude.at(Ilat_coarse) + coarse_data.latitude.at(Ilat_coarse - 1) );
                    }
                }
                if ( FINE_LAT_GRID_INCREASING ) {
                    TOP =  std::lower_bound( fine_data.latitude.begin(), fine_data.latitude.end(), target_lat ) 
                            - fine_data.latitude.begin();
                } else {
                    //TOP =  std::upper_bound( fine_data.latitude.begin(), fine_data.latitude.end(), target_lat ) 
                    TOP =  std::upper_bound( fine_data.latitude.rbegin(), fine_data.latitude.rend(), target_lat ) 
                            - fine_data.latitude.rbegin();
                    TOP = (Nlat_fine - 1) - TOP;
                }
                TOP = (TOP < 0) ? 0 : (TOP >= Nlat_fine) ? Nlat_fine - 1 : TOP;

                // left
                if ( COARSE_LON_GRID_INCREASING ) {
                    if ( Ilon_coarse == 0 ) {
                        target_lon = coarse_data.longitude.at(Ilon_coarse);
                    } else {
                        target_lon = 0.5 * ( coarse_data.longitude.at(Ilon_coarse) + coarse_data.longitude.at(Ilon_coarse - 1) );
                    }
                } else {
                    if ( Ilon_coarse < Nlon_coarse - 1 ) {
                        target_lon = 0.5 * ( coarse_data.longitude.at(Ilon_coarse) + coarse_data.longitude.at(Ilon_coarse + 1) );
                    } else {
                        target_lon = coarse_data.longitude.at(Ilon_coarse);
                    }
                }
                if ( FINE_LON_GRID_INCREASING ) {
                    LEFT = std::lower_bound( fine_data.longitude.begin(), fine_data.longitude.end(), target_lon ) 
                                - fine_data.longitude.begin();
                } else {
                    //LEFT = std::upper_bound( fine_data.longitude.begin(), fine_data.longitude.end(), target_lon ) 
                    LEFT = std::upper_bound( fine_data.longitude.rbegin(), fine_data.longitude.rend(), target_lon ) 
                                - fine_data.longitude.rbegin();
                    LEFT = (Nlon_fine - 1) - LEFT;
                }
                LEFT = (LEFT < 0) ? 0 : (LEFT >= Nlon_fine) ? Nlon_fine - 1 : LEFT;

                //fprintf( stdout, " %g, %d", target_lon, LEFT );

                // right
                if ( COARSE_LON_GRID_INCREASING ) {
                    if ( Ilon_coarse < Nlon_coarse - 1 ) {
                        target_lon = 0.5 * ( coarse_data.longitude.at(Ilon_coarse) + coarse_data.longitude.at(Ilon_coarse + 1) );
                    } else {
                        target_lon = coarse_data.longitude.at(Ilon_coarse);
                    }
                } else {
                    if ( Ilon_coarse == 0 ) {
                        target_lon = coarse_data.longitude.at(Ilon_coarse);
                    } else {
                        target_lon = 0.5 * ( coarse_data.longitude.at(Ilon_coarse) + coarse_data.longitude.at(Ilon_coarse - 1) );
                    }
                }
                if ( FINE_LON_GRID_INCREASING ) {
                    RIGHT = std::lower_bound( fine_data.longitude.begin(), fine_data.longitude.end(), target_lon ) 
                                - fine_data.longitude.begin();
                } else {
                    //RIGHT = std::upper_bound( fine_data.longitude.begin(), fine_data.longitude.end(), target_lon ) 
                    RIGHT = std::upper_bound( fine_data.longitude.rbegin(), fine_data.longitude.rend(), target_lon ) 
                                - fine_data.longitude.rbegin();
                    RIGHT = (Nlon_fine - 1) - RIGHT;
                }
                RIGHT = (RIGHT < 0) ? 0 : (RIGHT >= Nlon_fine) ? Nlon_fine - 1 : RIGHT;

                //fprintf( stdout, " : %g, %d\n", target_lon, RIGHT );
                //fprintf( stdout, " %'zu : [ %d, %d, %d, %d ] \n", II_coarse, LEFT, BOT, TOP, RIGHT );

                // Loop through the fine grid and build up the coarsened value.
                interp_val = 0.;
                land_cnt = 0.;
                cnt = 0;
                for (Ilat_fine = std::min(BOT,TOP); Ilat_fine <= std::max(BOT,TOP); Ilat_fine++) {
                    for (Ilon_fine = std::min(LEFT,RIGHT); Ilon_fine <= std::max(LEFT,RIGHT); Ilon_fine++) {
                        II_fine = Index( Itime, Idepth, Ilat_fine, Ilon_fine, Ntime, Ndepth, Nlat_fine, Nlon_fine );
                        if ( fine_data.mask.at( II_fine ) ) {
                            interp_val += fine_data.variables.at( "fine_field" ).at( II_fine );
                            cnt++;
                        } else {
                            land_cnt++;
                        }
                    }
                }

                if ( (cnt == 0) and (land_cnt == 0) ) {
                    fprintf( stderr, "No points in fine-grid correspond to coarse grid (%d, %d)\n", Ilat_coarse, Ilon_coarse );
                }

                if ( (cnt == 0) or (cnt < land_cnt) ) {
                    mask_coarse.at(II_coarse) = false;
                    var_coarse.at(II_coarse) = constants::FILTER_OVER_LAND ? 0 : constants::fill_value;
                } else {
                    mask_coarse.at(II_coarse) = true;
                    var_coarse.at(II_coarse) = interp_val / cnt;
                }
                if (not(mask_coarse.at(II_coarse))) { coarse_mask_count++; }

                #if DEBUG >= 3
                fprintf( stdout, " %'zu : [ %d, %d, %d, %d ] : %g, %d, %d : %g %s \n", 
                        II_coarse, LEFT, BOT, TOP, RIGHT, interp_val, cnt, land_cnt, var_coarse.at(II_coarse), mask_coarse.at(II_coarse) ? "Water" : "Land" );
                #endif
            }
        }
        #if DEBUG >= 1
        if (wRank == 0) { fprintf( stdout, "Coarse var has %'zu masked values.\n", coarse_mask_count ); }
        coarse_mask_count = 0;
        for (II_coarse = 0; II_coarse < Npts_coarse; ++II_coarse) {
            if ( var_coarse.at(II_coarse) == 0. ) { coarse_mask_count++; }
        }
        if (wRank == 0) { fprintf( stdout, "Coarse var has %'zu zero values.\n", coarse_mask_count ); }
        #endif
        write_field_to_output( var_coarse, vars_in_output.at(Ivar), starts, counts, output_fname, &mask_coarse );
    }

    coarse_mask_count = 0;
    for (II_coarse = 0; II_coarse < Npts_coarse; ++II_coarse) {
        var_coarse.at(II_coarse) = mask_coarse.at(II_coarse) ? 1. : 0.;
        if ( mask_coarse.at(II_coarse) ) { coarse_mask_count++; }
    }
    write_field_to_output( var_coarse, "mask", starts, counts, output_fname, NULL );
    if (wRank == 0) { fprintf( stdout, "Coarse mask has %'zu water values.\n", coarse_mask_count ); }

    if (wRank == 0) {fprintf( stdout, "Storing seed count to file\n" ); }
    add_attr_to_file( "seed_count", full_Ntime, output_fname.c_str() );

    #if DEBUG >= 1
    fprintf(stdout, "Processor %d / %d waiting to finalize.\n", wRank + 1, wSize);
    #endif
    MPI_Finalize();
    return 0;
}
