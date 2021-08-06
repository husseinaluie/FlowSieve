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

    // first argument is the flag, second argument is default value (for when flag is not present)
    const std::string   &coarse_fname   = input.getCmdOption("--coarse_file",   "coarse.nc"),
                        &fine_fname     = input.getCmdOption("--fine_file",     "fine.nc"),
                        &output_fname   = input.getCmdOption("--output_file",   "coarse_vel.nc");

    const std::string   &time_dim_name      = input.getCmdOption("--time",        "time"),
                        &depth_dim_name     = input.getCmdOption("--depth",       "depth"),
                        &latitude_dim_name  = input.getCmdOption("--latitude",    "latitude"),
                        &longitude_dim_name = input.getCmdOption("--longitude",   "longitude");

    const std::string &latlon_in_degrees  = input.getCmdOption("--is_degrees",   "true");

    const std::string   &Nprocs_in_time_string  = input.getCmdOption("--Nprocs_in_time",  "1"),
                        &Nprocs_in_depth_string = input.getCmdOption("--Nprocs_in_depth", "1");
    const int   Nprocs_in_time_input  = stoi(Nprocs_in_time_string),
                Nprocs_in_depth_input = stoi(Nprocs_in_depth_string);

    //const std::string  &var_name_coarse = input.getCmdOption("--var_in_coarse", "F");
    //const std::string  &var_name_output = input.getCmdOption("--var_in_output", "seed");

    std::vector< std::string > vars_to_refine, vars_in_output;
    input.getListofStrings( vars_to_refine, "--input_variables" );
    input.getListofStrings( vars_in_output, "--output_variables" );
    const int Nvars = vars_to_refine.size();

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
    coarse_data.load_time(      time_dim_name,      coarse_fname );
    coarse_data.load_depth(     depth_dim_name,     coarse_fname );
    coarse_data.load_latitude(  latitude_dim_name,  coarse_fname );
    coarse_data.load_longitude( longitude_dim_name, coarse_fname );

    fine_data.load_time(      time_dim_name,      fine_fname );
    fine_data.load_depth(     depth_dim_name,     fine_fname );
    fine_data.load_latitude(  latitude_dim_name,  fine_fname );
    fine_data.load_longitude( longitude_dim_name, fine_fname );

    // Apply some cleaning to the processor allotments if necessary. 
    coarse_data.check_processor_divisions( Nprocs_in_time_input, Nprocs_in_depth_input );
     
    // Convert to radians, if appropriate
    if ( latlon_in_degrees == "true" ) {
        convert_coordinates( coarse_data.longitude, coarse_data.latitude );
        convert_coordinates( fine_data.longitude,   fine_data.latitude );
    }

    // Read in velocity fields
    coarse_data.load_variable( "coarse_field", vars_to_refine.at(0), coarse_fname, true, true );

    const int   full_Ntime  = coarse_data.full_Ntime,
                Ntime       = coarse_data.myCounts[0],
                Ndepth      = coarse_data.myCounts[1],
                Nlat_coarse = coarse_data.Nlat,
                Nlon_coarse = coarse_data.Nlon,
                Nlat_fine   = fine_data.Nlat,
                Nlon_fine   = fine_data.Nlon;

    size_t starts[4] = { coarse_data.myStarts.at(0), coarse_data.myStarts.at(1), 0,              0              };
    size_t counts[4] = { coarse_data.myCounts.at(0), coarse_data.myCounts.at(1), fine_data.Nlat, fine_data.Nlon };

    // Compute the area of each 'cell' which will be necessary for creating the output file
    if (wRank == 0) { fprintf( stdout, "Computing cell areas.\n" ); }
    fine_data.compute_cell_areas();

    // Initialize file and write out coarsened fields
    if (wRank == 0) { fprintf( stdout, "Preparing output file\n" ); }
    initialize_output_file( fine_data, vars_in_output, output_fname.c_str() );

    #if DEBUG >= 1
    if (wRank == 0) {
        fprintf( stdout, " c(%d,%d,%d,%d) -> f(%d,%d,%d,%d)\n", Ntime, Ndepth, Nlat_coarse, Nlon_coarse,
                                                                Ntime, Ndepth, Nlat_fine,   Nlon_fine );
    }
    #endif

    // Now coarsen the velocity fields
    const size_t    Npts_fine = Ntime * Ndepth * Nlat_fine * Nlon_fine,
                    Npts_coarse = coarse_data.variables.at("coarse_field").size();
    std::vector<double> var_fine(Npts_fine);

    // Next, the coarse velocities
    int Itime, Idepth, Ilat_fine, Ilon_fine, lat_lb, lon_lb, LEFT, RIGHT, BOT, TOP;
    double target_lat, target_lon, LR_perc, TB_perc, 
           BL_val, BR_val, TL_val, TR_val, L_interp, R_interp, interp_val;
    size_t II_fine, BL_coarse, BR_coarse, TL_coarse, TR_coarse;

    for ( int Ivar = 0; Ivar < Nvars; Ivar++ ) {

        coarse_data.load_variable( "coarse_field", vars_to_refine.at(Ivar), coarse_fname, true, true );

        std::fill( var_fine.begin(), var_fine.end(), 0. );

        #pragma omp parallel \
        default(none) \
        shared( coarse_data, fine_data, var_fine, vars_to_refine, stdout ) \
        private( lat_lb, lon_lb, target_lat, target_lon, Itime, Idepth, II_fine, Ilat_fine, Ilon_fine, \
                 RIGHT, LEFT, BOT, TOP, LR_perc, TB_perc, BL_coarse, BR_coarse, TL_coarse, TR_coarse, \
                 BL_val, BR_val, TL_val, TR_val, L_interp, R_interp, interp_val )
        {
            #pragma omp for collapse(1) schedule(static)
            for (II_fine = 0; II_fine < Npts_fine; ++II_fine) {

                Index1to4( II_fine, Itime, Idepth, Ilat_fine, Ilon_fine, Ntime, Ndepth, Nlat_fine, Nlon_fine );


                // lat_lb is the smallest index such that coarse_lat(lat_lb) >= fine_lat(Ilat_fine)
                target_lat = fine_data.latitude.at(Ilat_fine);
                lat_lb =  std::lower_bound( coarse_data.latitude.begin(), coarse_data.latitude.end(), target_lat ) - coarse_data.latitude.begin();
                lat_lb = (lat_lb < 0) ? 0 : (lat_lb >= Nlat_coarse) ? Nlat_coarse - 1 : lat_lb;


                // lon_lb is the smallest index such that coarse_lon(lon_lb) >= fine_lon(Ilon_fine)
                target_lon = fine_data.longitude.at(Ilon_fine);
                lon_lb =  std::lower_bound( coarse_data.longitude.begin(), coarse_data.longitude.end(), target_lon ) - coarse_data.longitude.begin();
                lon_lb = (lon_lb < 0) ? 0 : (lon_lb >= Nlon_coarse) ? Nlon_coarse - 1 : lon_lb;

                // Get the points for the bounding box in the coarse grid
                RIGHT   = lon_lb == 0 ? 1 : lon_lb;
                LEFT    = RIGHT - 1;
                LR_perc = ( target_lon - coarse_data.longitude.at(LEFT) ) / ( coarse_data.longitude.at(RIGHT) - coarse_data.longitude.at(LEFT) );

                TOP = lat_lb == 0 ? 1 : lat_lb;
                BOT = TOP - 1;
                TB_perc = ( target_lat - coarse_data.latitude.at(BOT) ) / ( coarse_data.latitude.at(TOP) - coarse_data.latitude.at(BOT) );

                // Get the corresponding indices in the coarse grid
                BL_coarse = Index( Itime, Idepth, BOT, LEFT,   Ntime, Ndepth, Nlat_coarse, Nlon_coarse );
                BR_coarse = Index( Itime, Idepth, BOT, RIGHT,  Ntime, Ndepth, Nlat_coarse, Nlon_coarse );
                TL_coarse = Index( Itime, Idepth, TOP, LEFT,   Ntime, Ndepth, Nlat_coarse, Nlon_coarse );
                TR_coarse = Index( Itime, Idepth, TOP, RIGHT,  Ntime, Ndepth, Nlat_coarse, Nlon_coarse );

                // Pull out the values
                BL_val = coarse_data.variables.at( "coarse_field" ).at( BL_coarse );
                BR_val = coarse_data.variables.at( "coarse_field" ).at( BR_coarse );
                TL_val = coarse_data.variables.at( "coarse_field" ).at( TL_coarse );
                TR_val = coarse_data.variables.at( "coarse_field" ).at( TR_coarse );

                // Get the interpolation value
                L_interp = BL_val * (1 - TB_perc) + TL_val * TB_perc;
                R_interp = BR_val * (1 - TB_perc) + TR_val * TB_perc;

                interp_val = L_interp * (1 - LR_perc) + R_interp * LR_perc;

                // And drop into the fine grid
                var_fine.at(II_fine) = interp_val;
            }
        }
        write_field_to_output( var_fine, vars_in_output.at(Ivar), starts, counts, output_fname, NULL );
    }

    if (wRank == 0) {fprintf( stdout, "Storing seed count to file\n" ); }
    add_attr_to_file( "seed_count", full_Ntime, output_fname.c_str() );

    #if DEBUG >= 1
    fprintf(stdout, "Processor %d / %d waiting to finalize.\n", wRank + 1, wSize);
    #endif
    MPI_Finalize();
    return 0;
}
