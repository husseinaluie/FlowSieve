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
    const std::string   &input_fname    = input.getCmdOption("--input_file",    "input.nc"),
                        &grid_fname     = input.getCmdOption("--grid_file",     "grid.nc"),
                        &output_fname   = input.getCmdOption("--output_file",   "output.nc");

    const std::string   &time_dim_name      = input.getCmdOption("--time",        "time"),
                        &depth_dim_name     = input.getCmdOption("--depth",       "depth"),
                        &latitude_dim_name  = input.getCmdOption("--latitude",    "latitude"),
                        &longitude_dim_name = input.getCmdOption("--longitude",   "longitude");

    const std::string &latlon_in_degrees  = input.getCmdOption("--is_degrees",   "true");

    const std::string   &Nprocs_in_time_string  = input.getCmdOption("--Nprocs_in_time",  "1"),
                        &Nprocs_in_depth_string = input.getCmdOption("--Nprocs_in_depth", "1"),
                        &Nlayers_string         = input.getCmdOption("--Nlayers", "5"),
                        &RBase_string           = input.getCmdOption("--Rlargest", "50e3");
    const int   Nprocs_in_time_input  = stoi(Nprocs_in_time_string),
                Nprocs_in_depth_input = stoi(Nprocs_in_depth_string),
                num_interp_layers     = stoi(Nlayers_string);
    const double RBase                = stod(RBase_string);

    std::vector< std::string > vars_to_interp, vars_in_output;
    input.getListofStrings( vars_to_interp, "--input_variables" );
    input.getListofStrings( vars_in_output, "--output_variables" );
    const int Nvars = vars_to_interp.size();

    // Print processor assignments
    const int max_threads = omp_get_max_threads();
    omp_set_num_threads( max_threads );

    // Print some header info, depending on debug level
    print_header_info();

    // Initialize dataset class instance
    dataset orig_data, targ_data;

    // Read in source data / get size information
    #if DEBUG >= 1
    if (wRank == 0) { fprintf(stdout, "Reading in source data.\n\n"); }
    #endif

    // Read in the grid coordinates
    orig_data.load_time(      time_dim_name,      input_fname );
    orig_data.load_depth(     depth_dim_name,     input_fname );
    // Will need to use special LLC grid loading function
    // ... will also need to write said function
    read_LLC_latlon_from_file( orig_data.latitude,  latitude_dim_name,  input_fname );
    read_LLC_latlon_from_file( orig_data.longitude, longitude_dim_name, input_fname );
    //orig_data.load_latitude(  latitude_dim_name,  input_fname );
    //orig_data.load_longitude( longitude_dim_name, input_fname );

    // Read in the grid coordinates
    targ_data.load_time(      time_dim_name,      input_fname );
    targ_data.load_depth(     depth_dim_name,     input_fname );
    targ_data.load_latitude(  latitude_dim_name,  grid_fname );
    targ_data.load_longitude( longitude_dim_name, grid_fname );

    // Apply some cleaning to the processor allotments if necessary. 
    orig_data.Nlat = 1;
    orig_data.Nlon = 1;
    orig_data.check_processor_divisions( Nprocs_in_time_input, Nprocs_in_depth_input );
    orig_data.Nprocs_in_time  = orig_data.Nprocs_in_time;
    orig_data.Nprocs_in_depth = orig_data.Nprocs_in_depth;
     
    // Convert to radians, if appropriate
    if ( latlon_in_degrees == "true" ) {
        convert_coordinates( orig_data.longitude,   orig_data.latitude );
        convert_coordinates( targ_data.longitude,   targ_data.latitude );
    }

    // Read in field to get size info
    orig_data.load_variable( "to_interp", vars_to_interp.at(0), input_fname, true, true );

    size_t starts[4] = { orig_data.myStarts.at(0), orig_data.myStarts.at(1), 0,              0              };
    size_t counts[4] = { orig_data.myCounts.at(0), orig_data.myCounts.at(1), targ_data.Nlat, targ_data.Nlon };

    // Initialize file and write out coarsened fields
    if (wRank == 0) { fprintf( stdout, "Preparing output file\n" ); }
    initialize_output_file( targ_data, vars_in_output, output_fname.c_str() );

    // Get the relevant sizes
    const size_t    Npts_in  = orig_data.latitude.size(),
                    Nlatlon_in = Npts_in / ( orig_data.myCounts.at(0) * orig_data.myCounts.at(1) ),
                    Npts_out = targ_data.Ntime * targ_data.Ndepth * targ_data.Nlat * targ_data.Nlon;

    std::vector<double> interped_field(Npts_out);

    // Next, the coarse velocities
    for ( int Ivar = 0; Ivar < Nvars; Ivar++ ) {

        // Load the data
        orig_data.load_variable( "to_interp", vars_to_interp.at(Ivar), input_fname, true, true );

        // Grid to store nearby points.
        size_t out_index, in_index, pt_index = 0;
        double local_dist = 0, proj_x, proj_y, proj_f, poleward_lat, del_lon_lim,
               src_lon, src_lat, target_lon, target_lat;
        bool near_pole;
        int Ngrid_spline;

        double lambdav = 0.000;

        size_t LL_ind, LR_ind, UL_ind, UR_ind;
        double LL_x, LL_y, LL_dist = 1e100, LL_f,
               LR_x, LR_y, LR_dist = 1e100, LR_f,
               UL_x, UL_y, UL_dist = 1e100, UL_f,
               UR_x, UR_y, UR_dist = 1e100, UR_f,
               L_f, R_f;
        double x1, x2, x3, x4, y1, y2, y3, y4,
               term1, term2, term3, term4, denom;


        // pragma omp this loop
        #pragma omp parallel \
        default(none) \
        shared( orig_data, targ_data, interped_field ) \
        private( in_index, out_index, proj_x, proj_y, proj_f, near_pole, poleward_lat, \
                 del_lon_lim, src_lon, src_lat, target_lon, target_lat, pt_index, local_dist, \
                 LL_ind, LL_x, LL_y, LL_dist, LL_f, LR_ind, LR_x, LR_y, LR_dist, LR_f, \
                 UL_ind, UL_x, UL_y, UL_dist, UL_f, UR_ind, UR_x, UR_y, UR_dist, UR_f, \
                 x1, x2, x3, x4, y1 ,y2, y3, y4, term1, term2, term3, term4, denom  ) \
        firstprivate( Npts_in, Npts_out, Nlatlon_in, RBase, stdout )
        {
            #pragma omp for collapse(1) schedule(guided)
            for ( out_index = 0; out_index < Npts_out; out_index++ ) {
                pt_index = 0;

                UR_dist = 1e100;
                UL_dist = 1e100;
                LR_dist = 1e100;
                LL_dist = 1e100;

                #if DEBUG > 0
                target_lon = targ_data.longitude.at(  out_index % targ_data.Nlon );
                target_lat = targ_data.latitude.at(  (out_index / targ_data.Nlon) % targ_data.Nlat );
                #else
                target_lon = targ_data.longitude[  out_index % targ_data.Nlon ];
                target_lat = targ_data.latitude[  (out_index / targ_data.Nlon) % targ_data.Nlat ];
                #endif

                for ( in_index = 0; in_index < Npts_in; in_index++ ) {

                    #if DEBUG > 0
                    src_lon = orig_data.longitude.at( in_index % Nlatlon_in );
                    src_lat = orig_data.latitude.at(  in_index % Nlatlon_in );
                    #else
                    src_lon = orig_data.longitude[ in_index % Nlatlon_in ];
                    src_lat = orig_data.latitude[  in_index % Nlatlon_in ];
                    #endif

                    // First, a lazy pruning by latitude to avoid excessive distance calculations
                    if ( fabs( target_lat - src_lat ) < (RBase / constants::R_earth) ) {

                        // Next, a slightly less laxy pruning by longitude
                        poleward_lat = fabs(target_lat) + (RBase / constants::R_earth);
                        near_pole = poleward_lat >= (M_PI * 89.8 / 180);
                        del_lon_lim = std::fmin( M_PI, RBase / ( cos(poleward_lat) * constants::R_earth ) );
                        if ( (near_pole) or ( cos( target_lon - src_lon ) > cos(del_lon_lim) ) ) {

                            near_sided_project(
                                    proj_x, proj_y,
                                    src_lon, src_lat,
                                    target_lon, target_lat,
                                    2e5
                                    );
                            local_dist = sqrt( proj_x*proj_x + proj_y*proj_y );

                            if ( (proj_x >= 0) and (proj_y > 0) and ( local_dist < UR_dist ) ) {
                                UR_x = proj_x;
                                UR_y = proj_y;
                                UR_ind = in_index;
                                UR_dist = local_dist;
                                #if DEBUG > 0
                                UR_f = orig_data.variables.at("to_interp").at( in_index % Nlatlon_in );
                                #else
                                UR_f = orig_data.variables["to_interp"][ in_index % Nlatlon_in ];
                                #endif
                            } else if ( (proj_x > 0) and (proj_y <= 0) and ( local_dist < UL_dist ) ) {
                                UL_x = proj_x;
                                UL_y = proj_y;
                                UL_ind = in_index;
                                UL_dist = local_dist;
                                #if DEBUG > 0
                                UL_f = orig_data.variables.at("to_interp").at( in_index % Nlatlon_in );
                                #else
                                UL_f = orig_data.variables["to_interp"][ in_index % Nlatlon_in ];
                                #endif
                            } else if ( (proj_x < 0) and (proj_y >= 0) and ( local_dist < LR_dist ) ) {
                                LR_x = proj_x;
                                LR_y = proj_y;
                                LR_ind = in_index;
                                LR_dist = local_dist;
                                #if DEBUG > 0
                                LR_f = orig_data.variables.at("to_interp").at( in_index % Nlatlon_in );
                                #else
                                LR_f = orig_data.variables["to_interp"][ in_index % Nlatlon_in ];
                                #endif
                            } else if ( (proj_x <= 0) and (proj_y < 0) and ( local_dist < LL_dist ) ) {
                                LL_x = proj_x;
                                LL_y = proj_y;
                                LL_ind = in_index;
                                LL_dist = local_dist;
                                #if DEBUG > 0
                                LL_f = orig_data.variables.at("to_interp").at( in_index % Nlatlon_in );
                                #else
                                LL_f = orig_data.variables["to_interp"][ in_index % Nlatlon_in ];
                                #endif
                            }
                        }
                    }
                }

                x1 = LL_x;
                x2 = LR_x;
                x3 = UL_x;
                x4 = UR_x;

                y1 = LL_y;
                y2 = LR_y;
                y3 = UL_y;
                y4 = UR_y;

                denom = x1 * (x2 * (y1 - y2) * (y3 - y4) - x3 * (y1 - y3) * (y2 - y4) + x4 * (y1 - y4) * (y2 - y3))
                          + x2 * x3 * (y1 - y4) * (y2 - y3) 
                          - x2 * x4 * (y1 - y3) * (y2 - y4) 
                          + x3 * x4 * (y1 - y2) * (y3 - y4);

                term1 = -x2 * x3 * y2 * y4 + x2 * x3 * y3 * y4 + x2 * x4 * y2 * y3 - x2 * x4 * y3 * y4 - x3 * x4 * y2 * y3 + x3 * x4 * y2 * y4;
                term2 =  x1 * x3 * y1 * y4 - x1 * x3 * y3 * y4 - x1 * x4 * y1 * y3 + x1 * x4 * y3 * y4 + x3 * x4 * y1 * y3 - x3 * x4 * y1 * y4;
                term3 = -x1 * x2 * y1 * y4 + x1 * x2 * y2 * y4 + x1 * x4 * y1 * y2 - x1 * x4 * y2 * y4 - x2 * x4 * y1 * y2 + x2 * x4 * y1 * y4; 
                term4 =  x1 * x2 * y1 * y3 - x1 * x2 * y2 * y3 - x1 * x3 * y1 * y2 + x1 * x3 * y2 * y3 + x2 * x3 * y1 * y2 - x2 * x3 * y1 * y3;

                /*
                fprintf( stdout, "  (%.2g, %.2g) ->  (%.2g, %.2g, %.2g, %.3g),  (%.2g, %.2g, %.2g, %.3g),  (%.2g, %.2g, %.2g, %.3g),  (%.2g, %.2g, %.2g, %.3g)\n",
                      target_lon * 180 / M_PI, target_lat * 180 / M_PI, 
                      LL_x, LL_y, LL_f, term1 / denom,
                      LR_x, LR_y, LR_f, term2 / denom,
                      UL_x, UL_y, UL_f, term3 / denom,
                      UR_x, UR_y, UR_f, term4 / denom
                      );
                */

                // is at point (0,0)
                interped_field.at(out_index) = ( term1 * LL_f + term2 * LR_f + term3 * UL_f + term4 * UR_f ) / denom;

                //fprintf( stdout, "num points in stencil: %'zu, (%.2g, %.2g) -> %.2g\n", pt_index, target_lon * 180 / M_PI, target_lat * 180 / M_PI, interped_field.at(out_index) );
            }
        }

        // Write the data
        write_field_to_output( interped_field, vars_in_output.at(Ivar), starts, counts, output_fname );
    }

    #if DEBUG >= 1
    fprintf(stdout, "Processor %d / %d waiting to finalize.\n", wRank + 1, wSize);
    #endif
    MPI_Finalize();
    return 0;
}
