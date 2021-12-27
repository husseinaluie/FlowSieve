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

    size_t starts[4] = { fine_data.myStarts.at(0), fine_data.myStarts.at(1), 0,                0                };
    size_t counts[4] = { fine_data.myCounts.at(0), fine_data.myCounts.at(1), coarse_data.Nlat, coarse_data.Nlon };

    // Compute the area of each 'cell' which will be necessary for creating the output file
    if (wRank == 0) { fprintf( stdout, "Computing cell areas.\n" ); }
    coarse_data.compute_cell_areas();

    // Initialize file and write out coarsened fields
    if (wRank == 0) { fprintf( stdout, "Preparing output file\n" ); }
    initialize_output_file( coarse_data, vars_in_output, output_fname.c_str() );

    #if DEBUG >= 1
    if (wRank == 0) {
        fprintf( stdout, " c(%d,%d,%d,%d) -> f(%d,%d,%d,%d)\n", full_Ntime, Ndepth, Nlat_fine,   Nlon_fine,
                                                                full_Ntime, Ndepth, Nlat_coarse, Nlon_coarse );
    }
    #endif

    // Now coarsen the velocity fields
    const size_t    Npts_coarse = Ntime * Ndepth * Nlat_coarse * Nlon_coarse,
                    Npts_fine = fine_data.variables.at("fine_field").size();
    std::vector<double> var_coarse(Npts_coarse);
    std::vector<bool> mask_coarse(Npts_coarse, true);

    // Next, the coarse velocities
    int cnt, land_cnt, Itime, Idepth, LEFT, RIGHT, BOT, TOP, Ilat_fine, Ilon_fine, Ilat_coarse, Ilon_coarse;
    double target_lat, target_lon, interp_val;
    size_t II_fine, II_coarse;

    for ( int Ivar = 0; Ivar < Nvars; Ivar++ ) {

        fine_data.load_variable( "fine_field", vars_to_refine.at(Ivar), fine_fname, true, true );

        std::fill( var_coarse.begin(), var_coarse.end(), 0. );

        #pragma omp parallel \
        default(none) \
        shared( coarse_data, fine_data, var_coarse, mask_coarse, vars_to_refine, stdout ) \
        private( target_lat, target_lon, Itime, Idepth, II_coarse, Ilat_coarse, Ilon_coarse, \
                 RIGHT, LEFT, BOT, TOP, II_fine, Ilat_fine, Ilon_fine, cnt, interp_val, land_cnt )
        {
            #pragma omp for collapse(1) schedule(static)
            for (II_coarse = 0; II_coarse < Npts_coarse; ++II_coarse) {

                Index1to4( II_coarse, Itime, Idepth, Ilat_coarse, Ilon_coarse, Ntime, Ndepth, Nlat_coarse, Nlon_coarse );


                //
                //// Get bounding box indices in the fine grid.
                //

                // bottom
                if ( Ilat_coarse == 0 ) {
                    target_lat = coarse_data.latitude.at(Ilat_coarse);
                } else {
                    target_lat = 0.5 * ( coarse_data.latitude.at(Ilat_coarse) + coarse_data.latitude.at(Ilat_coarse - 1) );
                }
                BOT =  std::lower_bound( fine_data.latitude.begin(), fine_data.latitude.end(), target_lat ) - fine_data.latitude.begin();
                BOT = (BOT < 0) ? 0 : (BOT >= Nlat_fine) ? Nlat_fine - 1 : BOT;

                // top
                if ( Ilat_coarse < Nlat_coarse - 1 ) {
                    target_lat = 0.5 * ( coarse_data.latitude.at(Ilat_coarse) + coarse_data.latitude.at(Ilat_coarse + 1) );
                } else {
                    target_lat = coarse_data.latitude.at(Ilat_coarse);
                }
                TOP =  std::lower_bound( fine_data.latitude.begin(), fine_data.latitude.end(), target_lat ) - fine_data.latitude.begin();
                TOP = (TOP < 0) ? 0 : (TOP >= Nlat_fine) ? Nlat_fine - 1 : TOP;

                // left
                if ( Ilon_coarse == 0 ) {
                    target_lon = coarse_data.longitude.at(Ilon_coarse);
                } else {
                    target_lon = 0.5 * ( coarse_data.longitude.at(Ilon_coarse) + coarse_data.longitude.at(Ilon_coarse - 1) );
                }
                LEFT =  std::lower_bound( fine_data.longitude.begin(), fine_data.longitude.end(), target_lon ) - fine_data.longitude.begin();
                LEFT = (LEFT < 0) ? 0 : (LEFT >= Nlon_fine) ? Nlon_fine - 1 : LEFT;

                // right
                if ( Ilon_coarse < Nlon_coarse - 1 ) {
                    target_lon = 0.5 * ( coarse_data.longitude.at(Ilon_coarse) + coarse_data.longitude.at(Ilon_coarse + 1) );
                } else {
                    target_lon = coarse_data.longitude.at(Ilon_coarse);
                }
                RIGHT =  std::lower_bound( fine_data.longitude.begin(), fine_data.longitude.end(), target_lon ) - fine_data.longitude.begin();
                RIGHT = (RIGHT < 0) ? 0 : (RIGHT >= Nlon_fine) ? Nlon_fine - 1 : RIGHT;

                // Loop through the fine grid and build up the coarsened value.
                interp_val = 0.;
                land_cnt = 0.;
                cnt = 0;
                for (Ilat_fine = BOT; Ilat_fine <= TOP; Ilat_fine++) {
                    for (Ilon_fine = LEFT; Ilon_fine <= RIGHT; Ilon_fine++) {
                        II_fine = Index( Itime, Idepth, Ilat_fine, Ilon_fine, Ntime, Ndepth, Nlat_fine, Nlon_fine );
                        if ( fine_data.mask.at( II_fine ) ) {
                            interp_val += fine_data.variables.at( "fine_field" ).at( II_fine );
                        } else {
                            land_cnt++;
                        }
                        cnt++;
                    }
                }

                // And drop into the coarse grid
                var_coarse.at(II_coarse) = ( cnt == 0 ) ? 0 : ( interp_val / cnt );

                // If half or more of the points were land, then the new cell is also land
                mask_coarse.at(II_coarse) = ( double(land_cnt) / double(cnt) < 0.5 );

                #if DEBUG >= 3
                fprintf( stdout, " %'zu : [ %d, %d, %d, %d ] : %g, %d : %g \n", II_coarse, LEFT, BOT, TOP, RIGHT, interp_val, cnt, var_coarse.at(II_coarse) );
                #endif
            }
        }
        write_field_to_output( var_coarse, vars_in_output.at(Ivar), starts, counts, output_fname, &mask_coarse );
    }

    if (wRank == 0) {fprintf( stdout, "Storing seed count to file\n" ); }
    add_attr_to_file( "seed_count", full_Ntime, output_fname.c_str() );

    #if DEBUG >= 1
    fprintf(stdout, "Processor %d / %d waiting to finalize.\n", wRank + 1, wSize);
    #endif
    MPI_Finalize();
    return 0;
}
