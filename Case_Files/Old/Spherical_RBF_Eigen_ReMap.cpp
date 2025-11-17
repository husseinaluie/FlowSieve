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

#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/IterativeSolvers>

void refine_field(
        std::vector<double> & fine_var,
        std::vector<bool> & fine_mask,
        const std::vector<double> & coarse_var,
        const std::vector<bool> & coarse_mask,
        const std::vector<double> & coarse_lat,
        const std::vector<double> & coarse_lon,
        const std::vector<double> & fine_lat,
        const std::vector<double> & fine_lon
        ){

    int Itime, Idepth, Ilat_fine, Ilon_fine, lat_lb, lon_lb, LEFT, RIGHT, BOT, TOP;
    double target_lat, target_lon, LR_perc, TB_perc, 
           BL_val, BR_val, TL_val, TR_val, L_interp, R_interp, interp_val;
    bool BL_mask, BR_mask, TL_mask, TR_mask, L_mask, R_mask;
    size_t II_fine, BL_coarse, BR_coarse, TL_coarse, TR_coarse;

    const size_t Npts_fine = fine_lat.size() * fine_lon.size();
    const int Nlon_coarse = coarse_lon.size(),
              Nlat_coarse = coarse_lat.size(),
              Nlon_fine   = fine_lon.size(),
              Nlat_fine   = fine_lat.size(),
              Ntime       = 1,
              Ndepth      = 1;

    fine_var.resize( Npts_fine );
    fine_mask.resize( Npts_fine );
    std::fill( fine_var.begin(), fine_var.end(), 0. );
    std::fill( fine_mask.begin(), fine_mask.end(), false );

    const bool  COARSE_LAT_GRID_INCREASING = true,
                COARSE_LON_GRID_INCREASING = true;

    #pragma omp parallel \
    default(none) \
    shared( coarse_var, coarse_mask, fine_var, fine_mask, coarse_lat, coarse_lon, fine_lat, fine_lon ) \
    private( lat_lb, lon_lb, target_lat, target_lon, Itime, Idepth, II_fine, Ilat_fine, Ilon_fine, \
             RIGHT, LEFT, BOT, TOP, LR_perc, TB_perc, BL_coarse, BR_coarse, TL_coarse, TR_coarse, \
             BL_val, BR_val, TL_val, TR_val, L_interp, R_interp, interp_val, \
             BL_mask, BR_mask, TL_mask, TR_mask, L_mask, R_mask ) \
    firstprivate( Npts_fine, Nlon_fine, Nlat_fine, Nlon_coarse, Nlat_coarse, Ndepth, Ntime, \
                  COARSE_LAT_GRID_INCREASING, COARSE_LON_GRID_INCREASING )
    {
        #pragma omp for collapse(1) schedule(static)
        for (II_fine = 0; II_fine < Npts_fine; ++II_fine) {

            Index1to4( II_fine, Itime, Idepth, Ilat_fine, Ilon_fine, Ntime, Ndepth, Nlat_fine, Nlon_fine );


            target_lat = fine_lat.at(Ilat_fine);
            if ( COARSE_LAT_GRID_INCREASING ) {
                // lat_lb is the smallest index such that coarse_lat(lat_lb) >= fine_lat(Ilat_fine)
                lat_lb = std::lower_bound( coarse_lat.begin(), coarse_lat.end(), target_lat ) 
                            - coarse_lat.begin();
            } else {
                // lat_lb is the smallest index such that coarse_lat(lat_lb) < fine_lat(Ilat_fine)
                lat_lb = std::lower_bound( coarse_lat.rbegin(), coarse_lat.rend(), target_lat ) 
                            - coarse_lat.rbegin();
                lat_lb = (Nlat_coarse - 1) - lat_lb;
            }
            lat_lb = (lat_lb < 0) ? 0 : (lat_lb >= Nlat_coarse) ? Nlat_coarse - 1 : lat_lb;


            target_lon = fine_lon.at(Ilon_fine);
            if ( COARSE_LON_GRID_INCREASING ) {
                // lon_lb is the smallest index such that coarse_lon(lon_lb) >= fine_lon(Ilon_fine)
                lon_lb = std::lower_bound( coarse_lon.begin(), coarse_lon.end(), target_lon ) 
                            - coarse_lon.begin();
            } else {
                // lon_lb is the smallest index such that coarse_lon(lon_lb) < fine_lon(Ilon_fine)
                lon_lb = std::lower_bound( coarse_lon.rbegin(), coarse_lon.rend(), target_lon ) 
                            - coarse_lon.rbegin();
                lon_lb = (Nlon_coarse - 1) - lon_lb;
            }
            lon_lb = (lon_lb < 0) ? 0 : (lon_lb >= Nlon_coarse) ? Nlon_coarse - 1 : lon_lb;

            // Get the points for the bounding box in the coarse grid
            if ( COARSE_LON_GRID_INCREASING ) {
                RIGHT   = lon_lb == 0 ? 1 : lon_lb;
                LEFT    = RIGHT - 1;
            } else {
                LEFT  = lon_lb == 0 ? 1 : lon_lb;
                RIGHT = LEFT - 1;
            }
            LR_perc = ( target_lon - coarse_lon.at(LEFT) ) / ( coarse_lon.at(RIGHT) - coarse_lon.at(LEFT) );

            if ( COARSE_LAT_GRID_INCREASING ) {
                TOP = lat_lb == 0 ? 1 : lat_lb;
                BOT = TOP - 1;
            } else {
                BOT = lat_lb == 0 ? 1 : lat_lb;
                TOP = BOT - 1;
            }
            TB_perc = ( target_lat - coarse_lat.at(BOT) ) / ( coarse_lat.at(TOP) - coarse_lat.at(BOT) );

            // Get the corresponding indices in the coarse grid
            BL_coarse = Index( Itime, Idepth, BOT, LEFT,   Ntime, Ndepth, Nlat_coarse, Nlon_coarse );
            BR_coarse = Index( Itime, Idepth, BOT, RIGHT,  Ntime, Ndepth, Nlat_coarse, Nlon_coarse );
            TL_coarse = Index( Itime, Idepth, TOP, LEFT,   Ntime, Ndepth, Nlat_coarse, Nlon_coarse );
            TR_coarse = Index( Itime, Idepth, TOP, RIGHT,  Ntime, Ndepth, Nlat_coarse, Nlon_coarse );

            // Pull out the values
            BL_val = coarse_var.at( BL_coarse );
            BR_val = coarse_var.at( BR_coarse );
            TL_val = coarse_var.at( TL_coarse );
            TR_val = coarse_var.at( TR_coarse );

            // Handle different mask arrangements
            BL_mask = coarse_mask.at( BL_coarse );
            BR_mask = coarse_mask.at( BR_coarse );
            TL_mask = coarse_mask.at( TL_coarse );
            TR_mask = coarse_mask.at( TR_coarse );

            // Get the interpolation value
            L_mask = true;
            if ( BL_mask and TL_mask ) {
                L_interp = BL_val * (1 - TB_perc) + TL_val * TB_perc;
            } else if ( BL_mask ) {
                L_interp = BL_val;
            } else if ( TL_mask ) {
                L_interp = TL_val;
            } else {
                L_mask = false;
            }

            R_mask = true;
            if ( BR_mask and TR_mask ) {
                R_interp = BR_val * (1 - TB_perc) + TR_val * TB_perc;
            } else if ( BR_mask ) {
                R_interp = BR_val;
            } else if ( TR_mask ) {
                R_interp = TR_val;
            } else {
                R_mask = false;
            }

            // If there's any water points, use those to interpolate,
            // otherwise, it's land
            if ( L_mask and R_mask ) {
                interp_val = L_interp * (1 - LR_perc) + R_interp * LR_perc;
            } else if ( L_mask ) {
                interp_val = L_interp;
            } else if ( R_mask ) {
                interp_val = R_interp;
            } else {
                fine_mask.at(II_fine) = false; // otherwise, it's just land
                interp_val = constants::fill_value;
            }

            // And drop into the fine grid
            fine_var.at(II_fine) = interp_val;
        }
    }
}


void downsample_field( 
        std::vector<double> & coarse_var,
        std::vector<bool> & coarse_mask,
        const std::vector<double> & fine_var,
        const std::vector<bool> & fine_mask,
        const std::vector<double> & fine_lat,
        const std::vector<double> & fine_lon,
        const std::vector<double> & coarse_lat,
        const std::vector<double> & coarse_lon
        ){

    // Initialize fields
    const size_t Npts_coarse = coarse_lat.size() * coarse_lon.size();
    const int Nlon_coarse = coarse_lon.size(),
              Nlat_coarse = coarse_lat.size(),
              Nlon_fine   = fine_lon.size(),
              Nlat_fine   = fine_lat.size(),
              Ntime       = 1,
              Ndepth      = 1;
    coarse_var.resize(  Npts_coarse );
    coarse_mask.resize( Npts_coarse );
    std::fill( coarse_var.begin(),  coarse_var.end(), 0. );
    std::fill( coarse_mask.begin(), coarse_mask.end(), false );

    double target_lat, target_lon, interp_val;
    size_t II_fine, II_coarse;
    int cnt, land_cnt, Itime, Idepth, LEFT, RIGHT, BOT, TOP, Ilat_fine, Ilon_fine, Ilat_coarse, Ilon_coarse;
    #pragma omp parallel \
    default(none) \
    shared( coarse_var, coarse_mask, fine_var, fine_mask, \
            coarse_lat, coarse_lon, fine_lat, fine_lon ) \
    private( target_lat, target_lon, Itime, Idepth, II_coarse, Ilat_coarse, Ilon_coarse, \
             RIGHT, LEFT, BOT, TOP, II_fine, Ilat_fine, Ilon_fine, cnt, interp_val, land_cnt ) \
    firstprivate( Npts_coarse, Nlon_coarse, Nlat_coarse, Nlon_fine, Nlat_fine, Ntime, Ndepth, stderr )
    {
        #pragma omp for collapse(1) schedule(guided)
        for (II_coarse = 0; II_coarse < Npts_coarse; ++II_coarse) {

            Index1to4( II_coarse, Itime, Idepth, Ilat_coarse, Ilon_coarse, Ntime, Ndepth, Nlat_coarse, Nlon_coarse );

            // bottom
            if ( Ilat_coarse == 0 ) {
                target_lat = coarse_lat.at(Ilat_coarse);
            } else {
                target_lat = 0.5 * ( coarse_lat.at(Ilat_coarse) + coarse_lat.at(Ilat_coarse - 1) );
            }
            BOT = std::lower_bound( fine_lat.begin(), fine_lat.end(), target_lat ) 
                  - fine_lat.begin();
            BOT = (BOT < 0) ? 0 : (BOT >= Nlat_fine) ? Nlat_fine - 1 : BOT;

            // top
            if ( Ilat_coarse < Nlat_coarse - 1 ) {
                target_lat = 0.5 * ( coarse_lat.at(Ilat_coarse) + coarse_lat.at(Ilat_coarse + 1) );
            } else {
                target_lat = coarse_lat.at(Ilat_coarse);
            }
            TOP =  std::lower_bound( fine_lat.begin(), fine_lat.end(), target_lat ) 
                    - fine_lat.begin();
            TOP = (TOP < 0) ? 0 : (TOP >= Nlat_fine) ? Nlat_fine - 1 : TOP;

            // left
            if ( Ilon_coarse == 0 ) {
                target_lon = coarse_lon.at(Ilon_coarse);
            } else {
                target_lon = 0.5 * ( coarse_lon.at(Ilon_coarse) + coarse_lon.at(Ilon_coarse - 1) );
            }
            LEFT = std::lower_bound( fine_lon.begin(), fine_lon.end(), target_lon ) 
                    - fine_lon.begin();
            LEFT = (LEFT < 0) ? 0 : (LEFT >= Nlon_fine) ? Nlon_fine - 1 : LEFT;

            // right
            if ( Ilon_coarse < Nlon_coarse - 1 ) {
                target_lon = 0.5 * ( coarse_lon.at(Ilon_coarse) + coarse_lon.at(Ilon_coarse + 1) );
            } else {
                target_lon = coarse_lon.at(Ilon_coarse);
            }
            RIGHT = std::lower_bound( fine_lon.begin(), fine_lon.end(), target_lon ) 
                    - fine_lon.begin();
            RIGHT = (RIGHT < 0) ? 0 : (RIGHT >= Nlon_fine) ? Nlon_fine - 1 : RIGHT;

            // Loop through the fine grid and build up the coarsened value.
            interp_val = 0.;
            land_cnt = 0.;
            cnt = 0;
            for (Ilat_fine = std::min(BOT,TOP); Ilat_fine <= std::max(BOT,TOP); Ilat_fine++) {
                for (Ilon_fine = std::min(LEFT,RIGHT); Ilon_fine <= std::max(LEFT,RIGHT); Ilon_fine++) {
                    II_fine = Index( Itime, Idepth, Ilat_fine, Ilon_fine, Ntime, Ndepth, Nlat_fine, Nlon_fine );
                    if ( fine_mask.at( II_fine ) ) {
                        interp_val += fine_var.at( II_fine );
                        cnt++;
                    } else {
                        land_cnt++;
                    }
                }
            }

            if ( (cnt == 0) and (land_cnt == 0) ) {
                fprintf( stderr, "No points in fine-grid correspond to coarse grid (%d, %d)\n", Ilat_coarse, Ilon_coarse );
            }

            // Set to land if more than half the points in the averaging were land
            //if ( (cnt == 0) or (cnt < land_cnt) ) {
            // Set to land if no water points
            if ( cnt == 0 ) {
                coarse_mask.at(II_coarse) = false;
                coarse_var.at(II_coarse) = constants::FILTER_OVER_LAND ? 0 : constants::fill_value;
            } else {
                coarse_mask.at(II_coarse) = true;
                coarse_var.at(II_coarse) = interp_val / cnt;
            }
        }
    }
}






int main(int argc, char *argv[]) {
    
    // PERIODIC_Y implies UNIFORM_LAT_GRID
    static_assert( (constants::UNIFORM_LAT_GRID) or (not(constants::PERIODIC_Y)),
            "PERIODIC_Y requires UNIFORM_LAT_GRID.\n"
            "Please update constants.hpp accordingly.\n");

    // Cannot extend to poles AND be Cartesian
    static_assert( not( (constants::EXTEND_DOMAIN_TO_POLES) and (constants::CARTESIAN) ),
            "Cartesian implies that there are no poles, so cannot extend to poles."
            "Please update constants.hpp accordingly.");

    // We need to keep the mask, so cannot have constants::FILTER_OVER_LAND turned on
    static_assert( not( constants::FILTER_OVER_LAND ),
            "Spherical RBF requires the land mask, so FILTER_OVER_LAND must be turned OFF" );

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
                        &output_fname  = input.getCmdOption("--output_file",    "output.nc", asked_help);

    const std::string   &time_dim_name      = input.getCmdOption("--time",        "time",      asked_help),
                        &depth_dim_name     = input.getCmdOption("--depth",       "depth",     asked_help),
                        &latitude_dim_name  = input.getCmdOption("--latitude",    "latitude",  asked_help),
                        &longitude_dim_name = input.getCmdOption("--longitude",   "longitude", asked_help);

    const std::string &latlon_in_degrees  = input.getCmdOption("--is_degrees",   "true", asked_help);

    const std::string   &Nprocs_in_time_string      = input.getCmdOption("--Nprocs_in_time",  "1", asked_help),
                        &Nprocs_in_depth_string     = input.getCmdOption("--Nprocs_in_depth", "1", asked_help),
                        &num_interp_layers_string   = input.getCmdOption("--Num_interp_layers", "1", asked_help),
                        &max_interp_scale_string    = input.getCmdOption("--max_interp_scale", "2645e3", asked_help),
                        &skip_factor_string         = input.getCmdOption("--skip_factor", "10", asked_help),
                        &max_iterations_string      = input.getCmdOption("--max_iterations", "1e3", asked_help),
                        &tolerance_string           = input.getCmdOption("--tolerance", "1e-4", asked_help);
    const int   Nprocs_in_time_input  = stoi(Nprocs_in_time_string),
                Nprocs_in_depth_input = stoi(Nprocs_in_depth_string),
                num_interp_layers     = stoi(num_interp_layers_string),
                base_skip_factor      = stoi(skip_factor_string),
                max_iterations        = stod(max_iterations_string);

    const double max_interp_scale = stod(max_interp_scale_string),
                 tolerance        = stod(tolerance_string);

    // Also read in the fields to be filtered from commandline
    //   e.g. --variables "rho salinity p" (names must match with input netcdf file)
    std::vector< std::string > vars_to_interpolate;
    input.getListofStrings( vars_to_interpolate, "--variables", asked_help );
    const size_t Nvars = vars_to_interpolate.size();

    if (asked_help) { return 0; }

    // Print processor assignments
    const int max_threads = omp_get_max_threads();
    omp_set_num_threads( max_threads );

    // Print some header info, depending on debug level
    print_header_info();

    Timing_Records timing_records;

    // Initialize dataset class instance
    dataset source_data;

    int nEigen = Eigen::nbThreads( );
    fprintf( stdout, "Eigen will use %d threads.\n", nEigen );

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
    for (size_t field_ind = 0; field_ind < vars_to_interpolate.size(); field_ind++) {
        source_data.load_variable( vars_to_interpolate.at(field_ind), vars_to_interpolate.at(field_ind), input_fname, true, true );
    }

    // Get the MPI-local dimension sizes
    source_data.Ntime  = source_data.myCounts[0];
    source_data.Ndepth = source_data.myCounts[1];
    const size_t Npts = source_data.Ntime * source_data.Ndepth * source_data.Nlat * source_data.Nlon;

    const std::vector<int>  &myStarts = source_data.myStarts;

    //
    const int   Ntime   = source_data.Ntime,
                Ndepth  = source_data.Ndepth,
                Nlat    = source_data.Nlat,
                Nlon    = source_data.Nlon;

    const std::vector<double>   &latitude   = source_data.latitude,
                                &longitude  = source_data.longitude;

    const std::vector<bool> &mask = source_data.mask;

    // Initialize interpolation fields
    #if DEBUG >= 1
    if (wRank == 0) { fprintf( stdout, "Setting up %'zu interp fields.\n", Nvars ); fflush(stdout); }
    #endif
    std::vector< std::vector<double> > interpolated_fields(Nvars);
    std::vector<double> fine_interp, coarse_interp, coarse_lat, coarse_lon, coarse_var;
    std::vector<bool> fine_mask, coarse_mask;
    for (size_t field_ind = 0; field_ind < Nvars; field_ind++) {
        interpolated_fields.at(field_ind).resize( Npts, 0. );
    }


    // Coastline bool arrays [actually a bit-array under the hood]
    std::vector<bool> is_coastal_land, is_coastal_water;

    // various indices etc
    int Itime, Idepth, Ilat, Ilon, LAT, LON, curr_lat, curr_lon;
    double lat_at_curr, lat_at_ilat;
    size_t index, inner_index, i, j, i_land, i_water, Ivar;
    bool is_coast;

    double clock_on;

    for ( int interp_layer = 0; interp_layer < num_interp_layers; interp_layer++ ) {

        // Get the length-scale for this layer of interpolation
        double interp_scale = max_interp_scale / pow( 2, interp_layer );
        const int Nlon_coarse = std::min( Nlon, (int) ceil( 2 * M_PI * constants::R_earth * base_skip_factor / interp_scale ) );
        const int Nlat_coarse = std::min( Nlat, (int) ceil( M_PI * constants::R_earth * base_skip_factor / interp_scale ) - 1 );

        fprintf( stdout, "\n\nApplying RBF layer %d - %gkm, with Nlat/Nlon %d / %d\n", interp_layer+1, interp_scale/1e3, Nlat_coarse, Nlon_coarse );

        // Downsample onto the desired resolution grid
        coarse_lat.resize( Nlat_coarse );
        for ( int Itmp = 0; Itmp < Nlat_coarse; Itmp++ ) {
            coarse_lat[Itmp] = (Itmp+1) * M_PI / (Nlat_coarse+1) - (M_PI/2.);
        }
        coarse_lon.resize( Nlon_coarse );
        for ( int Itmp = 0; Itmp < Nlon_coarse; Itmp++ ) {
            coarse_lon[Itmp] = Itmp * 2 * M_PI / Nlon_coarse - M_PI;
        }
        downsample_field( coarse_var, coarse_mask, 
                          source_data.variables.at(vars_to_interpolate[0]), mask,
                          latitude, longitude,
                          coarse_lat, coarse_lon );
        size_t Npts_coarse = coarse_var.size();
        coarse_interp.resize( Npts_coarse);

        //
        //// Get the 'coastline' mask for this layer
        //
        fprintf(stdout, "Mapping the coastline.\n");
        clock_on = MPI_Wtime();
        is_coastal_land.resize( Npts_coarse );
        is_coastal_water.resize( Npts_coarse );
        std::fill( is_coastal_land.begin(), is_coastal_land.end(), 0);
        std::fill( is_coastal_water.begin(), is_coastal_water.end(), 0);

        size_t num_coastal_land = 0,
               num_coastal_water = 0;
        int LAT_lb, LAT_ub, LON_ub, LON_lb;
        #pragma omp parallel \
        default(none) \
        shared( coarse_mask, coarse_lat, coarse_lon, is_coastal_land, is_coastal_water ) \
        private( Itime, Idepth, Ilat, Ilon, is_coast, index, inner_index, \
                 LAT_lb, LAT_ub, LON_lb, LON_ub, \
                 LAT, LON, lat_at_ilat, lat_at_curr, curr_lat, curr_lon) \
        firstprivate( Ntime, Ndepth, Nlat_coarse, Nlon_coarse, interp_scale ) \
        reduction( +:num_coastal_land,num_coastal_water )
        {
            #pragma omp for collapse(2) schedule(guided)
            for ( Ilat = 0; Ilat < Nlat_coarse; Ilat++ ) {
                for ( Ilon = 0; Ilon < Nlon_coarse; Ilon++ ) {
                    Itime = 0;  // just hard-code for now
                    Idepth = 0;

                    is_coast = false;

                    index = Index(Itime, Idepth, Ilat, Ilon,
                                  Ntime, Ndepth, Nlat_coarse, Nlon_coarse);
                    get_lat_bounds(LAT_lb, LAT_ub, coarse_lat, Ilat, 
                            1.35 * interp_scale / constants::KernPad );

                    lat_at_ilat = coarse_lat[Ilat];

                    for ( LAT = LAT_lb; LAT < LAT_ub; LAT++) {

                        // Handle periodicity if necessary
                        if (constants::PERIODIC_Y) { curr_lat = ( LAT % Nlat_coarse + Nlat_coarse ) % Nlat_coarse; }
                        else                       { curr_lat = LAT; }
                        lat_at_curr = coarse_lat.at(curr_lat);

                        get_lon_bounds(LON_lb, LON_ub, coarse_lon, Ilon, lat_at_ilat, lat_at_curr, 
                                1.35 * interp_scale / constants::KernPad );
                        for ( LON = LON_lb; LON < LON_ub; LON++ ) {

                            // Handle periodicity if necessary
                            if (constants::PERIODIC_X) { curr_lon = ( LON % Nlon_coarse + Nlon_coarse ) % Nlon_coarse; }
                            else                       { curr_lon = LON; }

                            inner_index = Index(Itime, Idepth, curr_lat, curr_lon, Ntime, Ndepth, Nlat_coarse, Nlon_coarse);

                            // If mask is different between the two points, then it's coastal!
                            if ( coarse_mask.at(index) xor coarse_mask.at(inner_index) ) {
                                is_coast = true;
                                break; // break out of LON loop
                            }
                        }
                        if ( is_coast ) { break; } // break out of LAT loop
                    }
                    if (is_coast) { 
                        if ( coarse_mask.at(index) ) {
                            is_coastal_water.at(index) = true;
                            num_coastal_water++;
                        } else {
                            is_coastal_land.at(index) = true;
                            num_coastal_land++;
                        }
                    }
                }
            }
        }
        // end of computing coastal mask
        fprintf(stdout, "Mapped the coastline (%gs). %'zu land points, %'zu water points.\n",
                MPI_Wtime() - clock_on, num_coastal_land, num_coastal_water);

        fprintf(stdout, "Building utility index arrays.\n");
        clock_on = MPI_Wtime();
        // for utility, build an array of coastal water and land indices
        std::vector< size_t > coastal_water_indices( num_coastal_water, 0 );
        index = 0;
        for ( size_t II = 0; II < Npts_coarse; II++ ) {
            if ( is_coastal_water[II] ) { 
                coastal_water_indices[index] = II; 
                index++;
            }
        }
        /*
        size_t II = 0, count = 0;
        clock_on = MPI_Wtime();
        for ( index = 0; index < num_coastal_water; index++ ) {
            if ( is_coastal_water[II] ) {
                if ( count == index ) {
                    coastal_water_indices[index] = II;
                }
                count++;
            }
            II++;
        }
        fprintf( stdout, " Non-threaded (inverted): %gs\n", MPI_Wtime() - clock_on );
        */

        std::vector< size_t > coastal_land_indices( num_coastal_land, 0 );
        index = 0;
        for ( size_t II = 0; II < Npts_coarse; II++ ) {
            if ( is_coastal_land[II] ) { 
                coastal_land_indices[index] = II; 
                index++;
            }
        }
        fprintf( stdout, "  took %gs\n", MPI_Wtime() - clock_on );


        // Count number of non-zero elements going into the Aij matrix
        double lat_i, lat_j, lon_i, lon_j, kern_val;
        size_t num_nonzero = 0;
        std::vector<unsigned int> nnz( num_coastal_water, 0 );
        clock_on = MPI_Wtime();
        #pragma omp parallel default(none) \
        shared( is_coastal_water, coastal_water_indices, coarse_lat, coarse_lon, nnz ) \
        private( i, j, index, inner_index, Itime, Idepth, Ilat, Ilon, curr_lat, curr_lon, \
                 LAT, LON, LAT_lb, LAT_ub, LON_ub, LON_lb, lon_i, lat_i, lon_j, lat_j, \
                 kern_val ) \
        firstprivate( Ntime, Ndepth, Nlat_coarse, Nlon_coarse, num_coastal_water, interp_scale ) \
        reduction( +:num_nonzero )
        {
            #pragma omp for collapse(1) schedule(static)
            for ( i = 0; i < num_coastal_water; i++ ) {
                Itime = 0;
                Idepth = 0;
                index = coastal_water_indices[i];
                Index1to4( index, Itime, Idepth, Ilat,        Ilon,
                                  Ntime, Ndepth, Nlat_coarse, Nlon_coarse );
                get_lat_bounds(LAT_lb, LAT_ub, coarse_lat, Ilat, 
                        1.35 * interp_scale / constants::KernPad );

                lon_i = coarse_lon[Ilon];
                lat_i = coarse_lat[Ilat];

                for ( LAT = LAT_lb; LAT < LAT_ub; LAT++ ) {

                    // Handle periodicity if necessary
                    if (constants::PERIODIC_Y) { curr_lat = ( LAT % Nlat_coarse + Nlat_coarse ) % Nlat_coarse; }
                    else                       { curr_lat = LAT; }
                    lat_j = coarse_lat.at(curr_lat);

                    get_lon_bounds(LON_lb, LON_ub, coarse_lon, Ilon, lat_i, lat_j, 
                            1.35 * interp_scale / constants::KernPad );
                    for ( LON = LON_lb; LON < LON_ub; LON++ ) {

                        // Handle periodicity if necessary
                        if (constants::PERIODIC_X) { curr_lon = ( LON % Nlon_coarse + Nlon_coarse ) % Nlon_coarse; }
                        else                       { curr_lon = LON; }
                        lon_j = coarse_lon[curr_lon];

                        inner_index = Index(Itime, Idepth, curr_lat, curr_lon, Ntime, Ndepth, Nlat_coarse, Nlon_coarse);

                        if ( not( is_coastal_water[inner_index] ) ) { continue; }

                        //if (inner_index < index) { continue; } // only store upper triangular part

                        //kern_val = (distance( lon_i, lat_i, lon_j, lat_j ) < interp_scale / 2 ) ? 1. : 0.;
                        /*
                        kern_val = kernel( distance( lon_i, lat_i, lon_j, lat_j ), interp_scale );
                        if ( kern_val > 1e-3 ) {
                            num_nonzero++;
                            nnz[i]++;
                        }
                        */
                        num_nonzero++;
                        nnz[i]++;
                    }
                }
            }
        }
        fprintf( stdout, "Total non-zeros: %'zu (%gs)\n", num_nonzero, MPI_Wtime() - clock_on );

        clock_on = MPI_Wtime();
        std::vector<size_t> prev_nnz(num_coastal_water, 0);
        for ( j = 1; j < num_coastal_water; j++ ) {
            prev_nnz[j] = prev_nnz[j-1] + nnz[j-1];
        }
        fprintf( stdout, "Built prev_nnz (%gs)\n", MPI_Wtime() - clock_on );

        // We can't do parallel writing to Aij directly,
        //  so build list of the entries in parallel, and then
        //  just copy them over in serial
        size_t Innz = 0;
        typedef Eigen::Triplet<double> T;
        clock_on = MPI_Wtime();
        std::vector<T> Aij_triplets(num_nonzero);
        #pragma omp parallel default(none) \
        shared( is_coastal_water, coastal_water_indices, coarse_lat, coarse_lon, \
                nnz, prev_nnz, Aij_triplets, stderr ) \
        private( i, j, i_water, index, inner_index, Innz, Itime, Idepth, Ilat, Ilon, curr_lat, curr_lon, \
                 LAT, LON, LAT_lb, LAT_ub, LON_ub, LON_lb, lon_i, lat_i, lon_j, lat_j, kern_val ) \
        firstprivate( Ntime, Ndepth, Nlat_coarse, Nlon_coarse, num_coastal_water, \
                      interp_scale )
        {
            #pragma omp for collapse(1) schedule(static)
            for ( i = 0; i < num_coastal_water; i++ ) {
                Itime = 0;
                Idepth = 0;
                index = coastal_water_indices[i];
                Index1to4( index, Itime, Idepth, Ilat,        Ilon,
                                  Ntime, Ndepth, Nlat_coarse, Nlon_coarse );
                get_lat_bounds(LAT_lb, LAT_ub, coarse_lat, Ilat, 
                        1.35 * interp_scale / constants::KernPad );

                lon_i = coarse_lon[Ilon];
                lat_i = coarse_lat[Ilat];

                Innz = 0;
                for ( LAT = LAT_lb; LAT < LAT_ub; LAT++ ) {

                    // Handle periodicity if necessary
                    if (constants::PERIODIC_Y) { curr_lat = ( LAT % Nlat_coarse + Nlat_coarse ) % Nlat_coarse; }
                    else                       { curr_lat = LAT; }
                    lat_j = coarse_lat.at(curr_lat);

                    get_lon_bounds(LON_lb, LON_ub, coarse_lon, Ilon, lat_i, lat_j, 
                            1.35 * interp_scale / constants::KernPad );
                    for ( LON = LON_lb; LON < LON_ub; LON++ ) {

                        // Handle periodicity if necessary
                        if (constants::PERIODIC_X) { curr_lon = ( LON % Nlon_coarse + Nlon_coarse ) % Nlon_coarse; }
                        else                       { curr_lon = LON; }
                        lon_j = coarse_lon[curr_lon];

                        inner_index = Index(Itime, Idepth, curr_lat, curr_lon, Ntime, Ndepth, Nlat_coarse, Nlon_coarse);

                        if ( not( is_coastal_water[inner_index] ) ) { continue; }

                        //if (inner_index < index) { continue; } // only store upper triangular part

                        //kern_val = (distance( lon_i, lat_i, lon_j, lat_j ) < interp_scale / 2 ) ? 1. : 0.;
                        kern_val = kernel( distance( lon_i, lat_i, lon_j, lat_j ), interp_scale );
                        //if ( kern_val > 1e-3 ) {
                        if ( true ) {

                            i_water = std::lower_bound( coastal_water_indices.begin(), 
                                                        coastal_water_indices.end(), 
                                                        inner_index ) 
                                - coastal_water_indices.begin();

                            Aij_triplets[prev_nnz[i] + Innz] = T( i, i_water, kern_val );
                            Innz++;
                        }
                    }
                }
                if (Innz != nnz[i] ) {
                    fprintf( stderr, "(i,nnz[i],Innz) = (%'zu, %'u, %'zu)\n", i, nnz[i], Innz );
                }
                assert(Innz == nnz[i]);
            }
        }
        fprintf( stdout, "Built the Aij triplets (%gs)\n", MPI_Wtime() - clock_on );

        //
        //// Build the linear system
        //
        std::vector<double> RHS_vector( num_coastal_water, 0. );
        Eigen::SparseMatrix<double> Aij( num_coastal_water, num_coastal_water );

        //fprintf(stdout, "Building the system.\n");
        Aij.setFromTriplets( Aij_triplets.begin(), Aij_triplets.end() );

        //fprintf(stdout, "Converting the system to crs.\n");
        Aij.makeCompressed();

        fprintf(stdout, "Solving the system!\n");
        clock_on = MPI_Wtime();
        //Eigen::ConjugateGradient< Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper > solver;
        Eigen::LeastSquaresConjugateGradient< Eigen::SparseMatrix<double> > solver;
        //Eigen::BiCGSTAB< Eigen::SparseMatrix<double> > solver;
        //Eigen::GMRES< Eigen::SparseMatrix<double> > solver;
        fprintf( stdout, "  Setting solver tolerance and max iterations to %g and %d.\n", tolerance, max_iterations );
        solver.setMaxIterations(max_iterations);
        solver.setTolerance(tolerance);
        solver.compute( Aij );
        if ( solver.info() != Eigen::Success ) {
            // decomposition failed
            fprintf( stderr, "!!Eigen decomposition failed.\n" );
            return -1;
        }
        fprintf( stdout, "  Solve took %gs.\n", MPI_Wtime() - clock_on );

        for ( Ivar = 0; Ivar < Nvars; Ivar++ ) {
            fprintf(stdout, "Applying to variable %zu of %zu.\n", Ivar+1, Nvars);

            // Do a bit of extra processing if this is the first round
            if ( interp_layer == 0 ) {
                // Get the global mean value
                double val_sum = 0, area_sum = 0, dA;
                #pragma omp parallel default(none) \
                shared( mask, source_data, vars_to_interpolate, dA ) \
                private( i ) \
                firstprivate( Ivar, Npts ) \
                reduction( +:val_sum,area_sum )
                {
                    #pragma omp for collapse(1) schedule(static)
                    for ( i = 0; i < Npts; i++ ) {
                        if ( mask[i] ) {
                            dA = source_data.areas.at(i);
                            val_sum += dA * source_data.variables.at(vars_to_interpolate[Ivar])[i];
                            area_sum += dA;
                        }
                    }
                }
                const double global_mean_val = val_sum / area_sum;

                // For anything outside of the largest coastline, just copy over the original
                // Otherwise, initialize with global mean
                #pragma omp parallel default(none) \
                shared( mask, source_data, vars_to_interpolate, interpolated_fields, is_coastal_water ) \
                private( i ) \
                firstprivate( Ivar, Npts, global_mean_val )
                {
                    #pragma omp for collapse(1) schedule(static)
                    for ( i = 0; i < Npts; i++ ) {
                        if ( mask[i] and not( is_coastal_water[i] ) ) {
                            // if water, but not coastal water
                            interpolated_fields[Ivar][i] = source_data.variables.at(vars_to_interpolate[Ivar])[i];
                        } else {
                            interpolated_fields[Ivar][i] = global_mean_val;
                        }
                    }
                }

                // And finally, subtract the global mean off from the original field
                #pragma omp parallel default(none) \
                shared( mask, source_data, vars_to_interpolate ) \
                private( i ) \
                firstprivate( Ivar, Npts, global_mean_val )
                {
                    #pragma omp for collapse(1) schedule(static)
                    for ( i = 0; i < Npts; i++ ) {
                        if ( mask[i] ) {
                            source_data.variables.at(vars_to_interpolate[Ivar])[i] -= global_mean_val;
                        }
                    }
                }
            }

            // Downsample the field
            if ( (Nlat_coarse < Nlat) or (Nlon_coarse < Nlon) ) {
                downsample_field( coarse_var, coarse_mask, 
                                  source_data.variables.at(vars_to_interpolate[Ivar]), mask,
                                  latitude, longitude,
                                  coarse_lat, coarse_lon );
            } else {
                coarse_var = source_data.variables.at(vars_to_interpolate[Ivar]);
                coarse_mask.resize( mask.size() );
            }
            std::fill( coarse_mask.begin(), coarse_mask.end(), true );


            // Set up RHS vector
            #pragma omp parallel default(none) \
            shared(RHS_vector, coarse_var, coastal_water_indices) \
            private( i ) \
            firstprivate( Ivar, num_coastal_water )
            {
                #pragma omp for collapse(1) schedule(static)
                for ( i = 0; i < num_coastal_water; i++ ) {
                    RHS_vector[i] = coarse_var[coastal_water_indices[i]];
                }
            }

            Eigen::Map<Eigen::VectorXd> RHS( &RHS_vector[0], num_coastal_water );

            // Apply the solved linear system
            clock_on = MPI_Wtime();
            fprintf(stdout, "  Applying the solved system.\n");
            Eigen::VectorXd F_Eigen = solver.solve( RHS );
            /*
            if ( solver.info() != Eigen::Success ) {
                // solving failed
                fprintf( stderr, "Eigen solving failed.\n" );
                return -1;
            } else {
                fprintf( stdout, "Solver converged after %ld iterations to error %g.\n", 
                         solver.iterations(), solver.error() );
            }
            */
            // We want the solution, even if it didn't full converge.
            fprintf( stdout, "    Solver converged after %ld iterations to error %g (%gs).\n", 
                    solver.iterations(), solver.error(), MPI_Wtime() - clock_on );
            std::vector<double> wj( F_Eigen.data(), F_Eigen.data() + num_coastal_water);
    
            // And now apply the interpolation weights to fill in the water parts
            fprintf(stdout, "  Applying interpolation weights.\n");
            clock_on = MPI_Wtime();
            double lat_land, lon_land, lat_water, lon_water;
            #pragma omp parallel \
            default(none) \
            shared( coarse_lat, coarse_lon, wj, coarse_interp, \
                    coastal_water_indices, coastal_land_indices, is_coastal_water ) \
            private( i_land, i_water, lat_land, lon_land, lat_water, lon_water, kern_val, \
                     Itime, Idepth, Ilat, Ilon, LAT, LON, curr_lat, curr_lon, \
                     LAT_ub, LAT_lb, LON_ub, LON_lb, index, inner_index ) \
            firstprivate( num_coastal_water, num_coastal_land, interp_scale, Ivar,\
                          Ntime, Ndepth, Nlat_coarse, Nlon_coarse )
            {
                #pragma omp for collapse(1) schedule(static)
                for ( i_land = 0; i_land < num_coastal_land; i_land++ ) {
                    index = coastal_land_indices[i_land];
                    Index1to4( index, Itime, Idepth, Ilat,        Ilon,
                                      Ntime, Ndepth, Nlat_coarse, Nlon_coarse );
                    get_lat_bounds(LAT_lb, LAT_ub, coarse_lat, Ilat, interp_scale );

                    lon_land = coarse_lon[Ilon];
                    lat_land = coarse_lat[Ilat];

                    for ( LAT = LAT_lb; LAT < LAT_ub; LAT++) {

                        // Handle periodicity if necessary
                        if (constants::PERIODIC_Y) { curr_lat = ( LAT % Nlat_coarse + Nlat_coarse ) % Nlat_coarse; }
                        else                       { curr_lat = LAT; }
                        lat_water = coarse_lat.at(curr_lat);

                        get_lon_bounds(LON_lb, LON_ub, coarse_lon, Ilon, lat_land, lat_water, interp_scale );
                        for ( LON = LON_lb; LON < LON_ub; LON++ ) {

                            // Handle periodicity if necessary
                            if (constants::PERIODIC_X) { curr_lon = ( LON % Nlon_coarse + Nlon_coarse ) % Nlon_coarse; }
                            else                       { curr_lon = LON; }
                            lon_water = coarse_lon[curr_lon];

                            inner_index = Index(Itime, Idepth, curr_lat, curr_lon, Ntime, Ndepth, Nlat_coarse, Nlon_coarse);

                            // If inner_index isn't on the coast, skip it
                            if (not(is_coastal_water[inner_index])) { continue; }

                            i_water = std::lower_bound( coastal_water_indices.begin(), 
                                                        coastal_water_indices.end(), 
                                                        inner_index ) 
                                - coastal_water_indices.begin();
    
                            kern_val = kernel( distance( lon_land,  lat_land, 
                                                         lon_water, lat_water ),
                                               interp_scale );
                            /*
                            kern_val = ( distance( lon_land,  lat_land, lon_water, lat_water ) < interp_scale / 2 )
                                        ? 1.
                                        : 0.;
                             */
    
                            coarse_interp[coastal_land_indices[i_land]] += wj[i_water] * kern_val;

                        }
                    }
                }
            }
            fprintf( stdout, "     took %gs.\n", MPI_Wtime() - clock_on );

            fprintf(stdout, "  Computing residual.\n");
            clock_on = MPI_Wtime();
            size_t j_water;
            // Get the residual of the fields for the next interpolation layer
            //  and apply the interpolator to coastal water
            #pragma omp parallel \
            default(none) \
            shared( coarse_lat, coarse_lon, wj, coarse_interp, coarse_var, \
                    Aij, coastal_water_indices, is_coastal_water ) \
            private( i_water, j_water, lat_i, lat_j, kern_val, \
                     Itime, Idepth, Ilat, Ilon, LAT, LON, curr_lat, curr_lon, \
                     LAT_ub, LAT_lb, LON_ub, LON_lb, index, inner_index ) \
            firstprivate( num_coastal_water, interp_scale, Ivar, Ntime, Ndepth, Nlat_coarse, Nlon_coarse )
            {
                #pragma omp for collapse(1) schedule(static)
                for ( i_water = 0; i_water < num_coastal_water; i_water++ ) {
                    index = coastal_water_indices[i_water];
                    Index1to4( index, Itime, Idepth, Ilat, Ilon,
                                      Ntime, Ndepth, Nlat_coarse, Nlon_coarse );
                    get_lat_bounds(LAT_lb, LAT_ub, coarse_lat, Ilat, interp_scale );
                    lat_i = coarse_lat[Ilat];

                    for ( LAT = LAT_lb; LAT < LAT_ub; LAT++) {

                        // Handle periodicity if necessary
                        if (constants::PERIODIC_Y) { curr_lat = ( LAT % Nlat_coarse + Nlat_coarse ) % Nlat_coarse; }
                        else                       { curr_lat = LAT; }
                        lat_j = coarse_lat[curr_lat];

                        get_lon_bounds(LON_lb, LON_ub, coarse_lon, Ilon, lat_i, lat_j, interp_scale );
                        for ( LON = LON_lb; LON < LON_ub; LON++ ) {

                            // Handle periodicity if necessary
                            if (constants::PERIODIC_X) { curr_lon = ( LON % Nlon_coarse + Nlon_coarse ) % Nlon_coarse; }
                            else                       { curr_lon = LON; }

                            inner_index = Index(Itime, Idepth, curr_lat, curr_lon, Ntime, Ndepth, Nlat_coarse, Nlon_coarse);

                            // If inner_index isn't on the coast, skip it
                            if (not(is_coastal_water[inner_index])) { continue; }

                            j_water = std::lower_bound( coastal_water_indices.begin(), 
                                                        coastal_water_indices.end(), 
                                                        inner_index ) 
                                - coastal_water_indices.begin();
    
                            kern_val = Aij.coeff( i_water, j_water );
    
                            coarse_interp[coastal_water_indices[i_water]] += wj[j_water] * kern_val;

                        }
                    }
                }
            }
            fprintf( stdout, "    took %gs\n", MPI_Wtime() - clock_on );

            // Refine the latest interp onto the full grid, and apply it
            clock_on = MPI_Wtime();
            if ( (Nlat_coarse < Nlat) or (Nlon_coarse < Nlon) ) {
                refine_field( fine_interp, fine_mask, coarse_interp, coarse_mask,
                              coarse_lat, coarse_lon, latitude, longitude );
            } else {
                fine_interp = coarse_interp;
            }
            #pragma omp parallel default(none) \
            shared( interpolated_fields, fine_interp, source_data, vars_to_interpolate, \
                    fine_mask, is_coastal_water ) \
            private( i ) \
            firstprivate( Ivar, Npts )
            {
                #pragma omp for collapse(1) schedule(static)
                for ( i = 0; i < Npts; i++ ) {
                    interpolated_fields[Ivar][i] += fine_interp[i];
                    if ( is_coastal_water[i] ) { // If it's coastal water, get the residual
                        source_data.variables.at(vars_to_interpolate[Ivar])[i] -= fine_interp[i];
                    }
                }
            }
            fprintf( stdout, "Subtracting for residual took %gs\n", MPI_Wtime() - clock_on );



        } // end of Ivar loop

    } // End of interp layers loop


    //
    //// Create output file
    //
    if (not(constants::NO_FULL_OUTPUTS)) {
        initialize_output_file( source_data, vars_to_interpolate, output_fname, -1 );

        add_attr_to_file( "max_interp_scale",   max_interp_scale,   output_fname );
        add_attr_to_file( "num_interp_layers",  num_interp_layers,  output_fname );
        add_attr_to_file( "max_iterations",     max_iterations,     output_fname );
        add_attr_to_file( "tolerance",          tolerance,          output_fname );
        add_attr_to_file( "base_skip_factor",   base_skip_factor,   output_fname );

        const int ndims = 4;
        size_t starts[ndims];
        starts[0] = size_t(myStarts.at(0));
        starts[1] = size_t(myStarts.at(1));
        starts[2] = size_t(myStarts.at(2));
        starts[3] = size_t(myStarts.at(3));
        size_t counts[ndims] = { size_t(Ntime), size_t(Ndepth), size_t(Nlat), size_t(Nlon) };

        for ( Ivar = 0; Ivar < Nvars; Ivar++ ) {
            write_field_to_output( interpolated_fields.at(Ivar), 
                                   vars_to_interpolate.at(Ivar), 
                                   starts, counts, output_fname, NULL );
        }
    }


    // DONE!
    #if DEBUG >= 1
    fprintf(stdout, "Processor %d / %d waiting to finalize.\n", wRank + 1, wSize);
    #endif
    MPI_Finalize();

    fprintf( stdout, "\n\nProcess Complete\n\n" );
    return 0;

}
