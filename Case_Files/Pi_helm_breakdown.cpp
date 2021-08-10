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
#include "../differentiation_tools.hpp"

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

    static_assert( constants::FILTER_OVER_LAND, "Must filter over land" );

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
    const std::string   &tor_input_fname   = input.getCmdOption("--toroidal_input_file",        "toroidal_projection.nc"),
                        &pot_input_fname   = input.getCmdOption("--potential_input_file",       "potential_projection.nc"),
                        &vel_input_fname   = input.getCmdOption("--velocity_input_file",        "toroidal_projection.nc");

    const std::string   &time_dim_name      = input.getCmdOption("--time",        "time"),
                        &depth_dim_name     = input.getCmdOption("--depth",       "depth"),
                        &latitude_dim_name  = input.getCmdOption("--latitude",    "latitude"),
                        &longitude_dim_name = input.getCmdOption("--longitude",   "longitude");

    const std::string &latlon_in_degrees  = input.getCmdOption("--is_degrees",   "true");

    const std::string   &Nprocs_in_time_string  = input.getCmdOption("--Nprocs_in_time",  "1"),
                        &Nprocs_in_depth_string = input.getCmdOption("--Nprocs_in_depth", "1");
    const int   Nprocs_in_time_input  = stoi(Nprocs_in_time_string),
                Nprocs_in_depth_input = stoi(Nprocs_in_depth_string);

    const std::string   &tor_field_var_name     = input.getCmdOption("--tor_field",     "F"),
                        &pot_field_var_name     = input.getCmdOption("--pot_field",     "F"),
                        &vel_field_var_name     = input.getCmdOption("--vel_field",     "u_lat");

    // Also read in the filter scales from the commandline
    //   e.g. --filter_scales "10.e3 150.76e3 1000e3" (units are in metres)
    std::vector<double> filter_scales;
    input.getFilterScales( filter_scales, "--filter_scales" );
    assert( filter_scales.size() == 1 );
    const double filter_scale = filter_scales[0];

    // Print processor assignments
    const int max_threads = omp_get_max_threads();
    omp_set_num_threads( max_threads );

    // Print some header info, depending on debug level
    print_header_info();

    // Initialize dataset class instance
    dataset source_data;

    // Read in source data / get size information
    #if DEBUG >= 1
    if (wRank == 0) { fprintf(stdout, "Reading in source data.\n\n"); }
    #endif

    // Read in the grid coordinates
    //   implicitely assume coordinates are the same between input files
    source_data.load_time(      time_dim_name,      tor_input_fname );
    source_data.load_depth(     depth_dim_name,     tor_input_fname );
    source_data.load_latitude(  latitude_dim_name,  tor_input_fname );
    source_data.load_longitude( longitude_dim_name, tor_input_fname );

    // Apply some cleaning to the processor allotments if necessary. 
    source_data.check_processor_divisions( Nprocs_in_time_input, Nprocs_in_depth_input );
     
    // Convert to radians, if appropriate
    if ( latlon_in_degrees == "true" ) {
        convert_coordinates( source_data.longitude, source_data.latitude );
    }

    // Compute the area of each 'cell' which will be necessary for integration
    source_data.compute_cell_areas();

    // Read in the toroidal and potential fields
    source_data.load_variable( "F_potential", pot_field_var_name, pot_input_fname, false, true );
    source_data.load_variable( "F_toroidal",  tor_field_var_name, tor_input_fname, false, true );

    // Get the MPI-local dimension sizes
    source_data.Ntime  = source_data.myCounts[0];
    source_data.Ndepth = source_data.myCounts[1];

    const int   Ntime   = source_data.Ntime,
                Ndepth  = source_data.Ndepth,
                Nlat    = source_data.Nlat,
                Nlon    = source_data.Nlon;

    const int ndims = 4;
    size_t starts[ndims] = { 0, 0, 0,            0            };
    size_t counts[ndims] = { 1, 1, size_t(Nlat), size_t(Nlon) };

    const size_t num_pts = source_data.Nlat * source_data.Nlon;

    // read in velocity to get the mask
    source_data.load_variable( "sample_velocity", vel_field_var_name, vel_input_fname, true, true);

    // Mask out the pole, if necessary (i.e. set lat = 90 to land)
    mask_out_pole( source_data.latitude, source_data.mask, source_data.Ntime, source_data.Ndepth, source_data.Nlat, source_data.Nlon );

    const std::vector<double>   &time       = source_data.time,
                                &depth      = source_data.depth,
                                &latitude   = source_data.latitude,
                                &longitude  = source_data.longitude,
                                &dAreas     = source_data.areas;

    const std::vector<double>   &F_potential    = source_data.variables.at("F_potential"),
                                &F_toroidal     = source_data.variables.at("F_toroidal");

    const std::vector<bool> &mask = source_data.mask;

    // Get cartesian versions
    std::vector<double> u_r_zero( num_pts, 0. ),
        u_x_tor( num_pts, 0. ), u_y_tor( num_pts, 0. ), u_z_tor( num_pts, 0. ),
        u_x_pot( num_pts, 0. ), u_y_pot( num_pts, 0. ), u_z_pot( num_pts, 0. ),
        u_x_tot( num_pts, 0. ), u_y_tot( num_pts, 0. ), u_z_tot( num_pts, 0. ),
        u_x_tor_bar( num_pts, 0. ), u_y_tor_bar( num_pts, 0. ), u_z_tor_bar( num_pts, 0. ),
        u_x_pot_bar( num_pts, 0. ), u_y_pot_bar( num_pts, 0. ), u_z_pot_bar( num_pts, 0. ),
        u_lon_tor( num_pts, 0. ), u_lat_tor( num_pts, 0. ), u_lon_pot( num_pts, 0. ), u_lat_pot( num_pts, 0. ),
        u_lon_tor_bar( num_pts, 0. ), u_lat_tor_bar( num_pts, 0. ), u_lon_pot_bar( num_pts, 0. ), u_lat_pot_bar( num_pts, 0. );

    // Get pot and tor velocities
    toroidal_vel_from_F(  u_lon_tor, u_lat_tor, F_toroidal,  longitude, latitude, Ntime, Ndepth, Nlat, Nlon, mask );
    potential_vel_from_F( u_lon_pot, u_lat_pot, F_potential, longitude, latitude, Ntime, Ndepth, Nlat, Nlon, mask );

    vel_Spher_to_Cart( u_x_tor, u_y_tor, u_z_tor, u_r_zero, u_lon_tor, u_lat_tor, source_data );
    vel_Spher_to_Cart( u_x_pot, u_y_pot, u_z_pot, u_r_zero, u_lon_pot, u_lat_pot, source_data );
    std::vector<double*> filtered_vals(6), x_deriv_vals, y_deriv_vals, z_deriv_vals;

    // Filtering variables
    std::vector<bool> filt_use_mask;
    std::vector<const std::vector<double>*> filter_fields, deriv_fields(4);

    size_t index;
    int ii, jj, Ilat, Ilon, LAT_lb, LAT_ub;
    const int Itime = 0, Idepth = 0;

    //
    //// First, filtering Phi and Psi
    //
    filter_fields.resize(2);
    filt_use_mask.resize(2, false);

    filter_fields[0] = &F_potential;
    filter_fields[1] = &F_toroidal;

    // Now filter Phi, Psi (to get ui, uj), and uiuj
    std::vector<double> Psi_coarse(num_pts), Phi_coarse(num_pts),
        uxux_tor(num_pts), uxuy_tor(num_pts), uxuz_tor(num_pts), uyuy_tor(num_pts), uyuz_tor(num_pts), uzuz_tor(num_pts),
        uxux_pot(num_pts), uxuy_pot(num_pts), uxuz_pot(num_pts), uyuy_pot(num_pts), uyuz_pot(num_pts), uzuz_pot(num_pts);
    std::vector<double> local_kernel(num_pts);
    double F_tor_tmp, F_pot_tmp, uxux_tmp, uxuy_tmp, uxuz_tmp, uyuy_tmp, uyuz_tmp, uzuz_tmp,
        vort_ux_tmp, vort_uy_tmp, vort_uz_tmp;
    for (Ilat = 0; Ilat < Nlat; Ilat++) {

        get_lat_bounds(LAT_lb, LAT_ub, source_data.latitude,  Ilat, filter_scale); 

        std::fill(local_kernel.begin(), local_kernel.end(), 0);
        compute_local_kernel( local_kernel, filter_scale, source_data, Ilat, 0, LAT_lb, LAT_ub );

        #pragma omp parallel default(none) \
        shared( filter_fields, source_data, Ilat, LAT_lb, LAT_ub, filt_use_mask, local_kernel, \
                Psi_coarse, Phi_coarse, mask, ii, jj, latitude, longitude, \
                uxux_tor, uxuy_tor, uxuz_tor, uyuy_tor, uyuz_tor, uzuz_tor, \
                uxux_pot, uxuy_pot, uxuz_pot, uyuy_pot, uyuz_pot, uzuz_pot, \
                u_x_tor, u_y_tor, u_z_tor, u_x_pot, u_y_pot, u_z_pot, u_r_zero ) \
        private( filtered_vals, Ilon, F_tor_tmp, F_pot_tmp, vort_ux_tmp, vort_uy_tmp, vort_uz_tmp, \
                uxux_tmp, uxuy_tmp, uxuz_tmp, uyuy_tmp, uyuz_tmp, uzuz_tmp, index )
        {

            filtered_vals.resize(2);

            filtered_vals[0] = &F_pot_tmp;
            filtered_vals[1] = &F_tor_tmp;

            #pragma omp for collapse(1) schedule(static)
            for (Ilon = 0; Ilon < Nlon; Ilon++) {

                index = Index( Itime, Idepth, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon );

                // Psi, Phi
                apply_filter_at_point( filtered_vals, filter_fields, source_data, Itime, Idepth, Ilat, Ilon, 
                        LAT_lb, LAT_ub, filter_scale, filt_use_mask, local_kernel );

                Psi_coarse.at( index ) = F_tor_tmp;
                Phi_coarse.at( index ) = F_pot_tmp;

                // uiuj (tor)
                apply_filter_at_point_for_quadratics(
                        uxux_tmp, uxuy_tmp, uxuz_tmp, uyuy_tmp, uyuz_tmp, uzuz_tmp, vort_ux_tmp, vort_uy_tmp, vort_uz_tmp,
                        u_x_tor,  u_y_tor,  u_z_tor, u_r_zero, source_data, Itime, Idepth, Ilat, Ilon,
                        LAT_lb, LAT_ub, filter_scale, local_kernel);

                uxux_tor.at(index) = uxux_tmp;
                uxuy_tor.at(index) = uxuy_tmp;
                uxuz_tor.at(index) = uxuz_tmp;
                uyuy_tor.at(index) = uyuy_tmp;
                uyuz_tor.at(index) = uyuz_tmp;
                uzuz_tor.at(index) = uzuz_tmp;

                // uiuj (pot)
                apply_filter_at_point_for_quadratics(
                        uxux_tmp, uxuy_tmp, uxuz_tmp, uyuy_tmp, uyuz_tmp, uzuz_tmp, vort_ux_tmp, vort_uy_tmp, vort_uz_tmp,
                        u_x_pot,  u_y_pot,  u_z_pot, u_r_zero, source_data, Itime, Idepth, Ilat, Ilon,
                        LAT_lb, LAT_ub, filter_scale, local_kernel);

                uxux_pot.at(index) = uxux_tmp;
                uxuy_pot.at(index) = uxuy_tmp;
                uxuz_pot.at(index) = uxuz_tmp;
                uyuy_pot.at(index) = uyuy_tmp;
                uyuz_pot.at(index) = uyuz_tmp;
                uzuz_pot.at(index) = uzuz_tmp;

            }
        }
    }

    // Get coarse pot and tor velocities, and Cartesian equivalents
    toroidal_vel_from_F(  u_lon_tor_bar, u_lat_tor_bar, Psi_coarse, longitude, latitude, Ntime, Ndepth, Nlat, Nlon, mask );
    potential_vel_from_F( u_lon_pot_bar, u_lat_pot_bar, Phi_coarse, longitude, latitude, Ntime, Ndepth, Nlat, Nlon, mask );

    vel_Spher_to_Cart( u_x_tor_bar, u_y_tor_bar, u_z_tor_bar, u_r_zero, u_lon_tor_bar, u_lat_tor_bar, source_data );
    vel_Spher_to_Cart( u_x_pot_bar, u_y_pot_bar, u_z_pot_bar, u_r_zero, u_lon_pot_bar, u_lat_pot_bar, source_data );

    std::vector< std::string > vars_to_write;
    vars_to_write.push_back( "Psi" );
    vars_to_write.push_back( "Phi" );
    vars_to_write.push_back( "ulon_tor" );
    vars_to_write.push_back( "ulat_tor" );
    vars_to_write.push_back( "ulon_pot" );
    vars_to_write.push_back( "ulat_pot" );
    vars_to_write.push_back( "uxux_tor" );
    vars_to_write.push_back( "uxuy_tor" );
    vars_to_write.push_back( "uxuz_tor" );
    vars_to_write.push_back( "uyuy_tor" );
    vars_to_write.push_back( "uyuz_tor" );
    vars_to_write.push_back( "uzuz_tor" );
    vars_to_write.push_back( "uxux_pot" );
    vars_to_write.push_back( "uxuy_pot" );
    vars_to_write.push_back( "uxuz_pot" );
    vars_to_write.push_back( "uyuy_pot" );
    vars_to_write.push_back( "uyuz_pot" );
    vars_to_write.push_back( "uzuz_pot" );


    std::string output_filename = "Pi_breakdown.nc";
    initialize_output_file( source_data, vars_to_write, output_filename.c_str(), filter_scale );

    write_field_to_output( Psi_coarse, "Psi", starts, counts, output_filename, &mask);
    write_field_to_output( Phi_coarse, "Phi", starts, counts, output_filename, &mask);

    write_field_to_output( u_lon_tor_bar, "ulon_tor", starts, counts, output_filename, &mask);
    write_field_to_output( u_lat_tor_bar, "ulat_tor", starts, counts, output_filename, &mask);
    write_field_to_output( u_lon_pot_bar, "ulon_pot", starts, counts, output_filename, &mask);
    write_field_to_output( u_lat_pot_bar, "ulat_pot", starts, counts, output_filename, &mask);

    write_field_to_output( uxux_tor, "uxux_tor", starts, counts, output_filename, &mask);
    write_field_to_output( uxuy_tor, "uxuy_tor", starts, counts, output_filename, &mask);
    write_field_to_output( uxuz_tor, "uxuz_tor", starts, counts, output_filename, &mask);
    write_field_to_output( uyuy_tor, "uyuy_tor", starts, counts, output_filename, &mask);
    write_field_to_output( uyuz_tor, "uyuz_tor", starts, counts, output_filename, &mask);
    write_field_to_output( uzuz_tor, "uzuz_tor", starts, counts, output_filename, &mask);

    write_field_to_output( uxux_pot, "uxux_pot", starts, counts, output_filename, &mask);
    write_field_to_output( uxuy_pot, "uxuy_pot", starts, counts, output_filename, &mask);
    write_field_to_output( uxuz_pot, "uxuz_pot", starts, counts, output_filename, &mask);
    write_field_to_output( uyuy_pot, "uyuy_pot", starts, counts, output_filename, &mask);
    write_field_to_output( uyuz_pot, "uyuz_pot", starts, counts, output_filename, &mask);
    write_field_to_output( uzuz_pot, "uzuz_pot", starts, counts, output_filename, &mask);

    //
    //// Now prepare to do that actual ui uj loops
    //

    double uu_tor_tor_tmp, uu_tor_pot_tmp, uu_pot_tor_tmp, uu_pot_pot_tmp;

    //
    vars_to_write.clear();
    vars_to_write.push_back( "ui_j_tor" );
    vars_to_write.push_back( "ui_j_pot" );
    vars_to_write.push_back( "uj_i_tor" );
    vars_to_write.push_back( "uj_i_pot" );
    vars_to_write.push_back( "ui_tor_coarse" );
    vars_to_write.push_back( "uj_tor_coarse" );
    vars_to_write.push_back( "ui_pot_coarse" );
    vars_to_write.push_back( "uj_pot_coarse" );
    vars_to_write.push_back( "uu_tor_tor_coarse" );
    vars_to_write.push_back( "uu_tor_pot_coarse" );
    vars_to_write.push_back( "uu_pot_tor_coarse" );
    vars_to_write.push_back( "uu_pot_pot_coarse" );

    const std::vector<double> *ui_tor, *ui_pot, *uj_tor, *uj_pot;
    std::vector<double> uu_tor_tor( num_pts, 0. ), uu_tor_pot( num_pts, 0. ), uu_pot_tor( num_pts, 0. ), uu_pot_pot( num_pts, 0. ), 
        ui_tor_coarse(num_pts), ui_pot_coarse(num_pts), uj_tor_coarse(num_pts), uj_pot_coarse(num_pts),
        ui_j_tor(num_pts), ui_j_pot(num_pts), uj_i_tor(num_pts), uj_i_pot(num_pts),
        uu_tor_tor_coarse(num_pts), uu_tor_pot_coarse(num_pts), uu_pot_tor_coarse(num_pts), uu_pot_pot_coarse(num_pts);
    double ui_j_tor_tmp, ui_j_pot_tmp, uj_i_tor_tmp, uj_i_pot_tmp;

    filter_fields.resize(4);
    filt_use_mask.resize(4, false);

    for (ii = 0; ii < 3; ii++) {

        fprintf( stdout, "ii = %d\n", ii );

        for (jj = 0; jj < 3; jj++) {

            fprintf( stdout, "  jj = %d\n", jj );

            // uiuj
            for ( index = 0; index < num_pts; index++ ) {
                double ui_T = (ii == 0) ? u_x_tor[index] : (ii == 1) ? u_y_tor[index] : u_z_tor[index];
                double ui_P = (ii == 0) ? u_x_pot[index] : (ii == 1) ? u_y_pot[index] : u_z_pot[index];
                double uj_T = (jj == 0) ? u_x_tor[index] : (jj == 1) ? u_y_tor[index] : u_z_tor[index];
                double uj_P = (jj == 0) ? u_x_pot[index] : (jj == 1) ? u_y_pot[index] : u_z_pot[index];
                uu_tor_tor.at(index) = ui_T * uj_T;
                uu_tor_pot.at(index) = ui_T * uj_P;
                uu_pot_tor.at(index) = ui_P * uj_T;
                uu_pot_pot.at(index) = ui_P * uj_P;
            }

            switch (ii) {
                case 0 :
                    switch (jj) {
                        case 0 : output_filename = "ux_ux.nc"; break;
                        case 1 : output_filename = "ux_uy.nc"; break;
                        case 2 : output_filename = "ux_uz.nc"; break;
                    };
                    break;
                case 1 :
                    switch (jj) {
                        case 0 : output_filename = "uy_ux.nc"; break;
                        case 1 : output_filename = "uy_uy.nc"; break;
                        case 2 : output_filename = "uy_uz.nc"; break;
                    };
                    break;
                case 2 :
                    switch (jj) {
                        case 0 : output_filename = "uz_ux.nc"; break;
                        case 1 : output_filename = "uz_uy.nc"; break;
                        case 2 : output_filename = "uz_uz.nc"; break;
                    };
                    break;
            }

            filter_fields[0] = &uu_tor_tor;
            filter_fields[1] = &uu_tor_pot;
            filter_fields[2] = &uu_pot_tor;
            filter_fields[3] = &uu_pot_pot;

            // Now filter Phi, Psi (to get ui, uj), and uiuj
            for (Ilat = 0; Ilat < Nlat; Ilat++) {

                get_lat_bounds(LAT_lb, LAT_ub, source_data.latitude,  Ilat, filter_scale); 

                std::fill(local_kernel.begin(), local_kernel.end(), 0);
                compute_local_kernel( local_kernel, filter_scale, source_data, Ilat, 0, LAT_lb, LAT_ub );

                #pragma omp parallel \
                default(none) \
                shared( filter_fields, deriv_fields, source_data, Ilat, LAT_lb, LAT_ub, filt_use_mask, local_kernel, \
                        uu_tor_tor_coarse, uu_tor_pot_coarse, uu_pot_tor_coarse, uu_pot_pot_coarse, \
                        mask, ii, jj, latitude, longitude ) \
                private( filtered_vals, Ilon, uu_tor_tor_tmp, uu_tor_pot_tmp, uu_pot_tor_tmp, uu_pot_pot_tmp )
                {

                    filtered_vals.resize(4);

                    filtered_vals[0] = &uu_tor_tor_tmp;
                    filtered_vals[1] = &uu_tor_pot_tmp;
                    filtered_vals[2] = &uu_pot_tor_tmp;
                    filtered_vals[3] = &uu_pot_pot_tmp;

                    #pragma omp for collapse(1) schedule(static)
                    for (Ilon = 0; Ilon < Nlon; Ilon++) {

                        apply_filter_at_point( filtered_vals, filter_fields, source_data, Itime, Idepth, Ilat, Ilon, 
                                LAT_lb, LAT_ub, filter_scale, filt_use_mask, local_kernel );

                        uu_tor_tor_coarse.at( Ilat*Nlon + Ilon ) = uu_tor_tor_tmp;
                        uu_tor_pot_coarse.at( Ilat*Nlon + Ilon ) = uu_tor_pot_tmp;
                        uu_pot_tor_coarse.at( Ilat*Nlon + Ilon ) = uu_pot_tor_tmp;
                        uu_pot_pot_coarse.at( Ilat*Nlon + Ilon ) = uu_pot_pot_tmp;

                    }
                }
            }

            // Now get derivatives of the velocities 
            deriv_fields[0] = (ii == 0) ? &u_x_tor_bar : (ii == 1) ? &u_y_tor_bar : &u_z_tor_bar;
            deriv_fields[1] = (ii == 0) ? &u_x_pot_bar : (ii == 1) ? &u_y_pot_bar : &u_z_pot_bar;
            deriv_fields[2] = (jj == 0) ? &u_x_tor_bar : (jj == 1) ? &u_y_tor_bar : &u_z_tor_bar;
            deriv_fields[3] = (jj == 0) ? &u_x_pot_bar : (jj == 1) ? &u_y_pot_bar : &u_z_pot_bar;

            for (Ilat = 0; Ilat < Nlat; Ilat++) {

                #pragma omp parallel \
                default(none) \
                shared( deriv_fields, Ilat, ui_j_tor, ui_j_pot, uj_i_tor, uj_i_pot, mask, ii, jj, latitude, longitude ) \
                private( Ilon, x_deriv_vals, y_deriv_vals, z_deriv_vals, ui_j_tor_tmp, ui_j_pot_tmp, uj_i_tor_tmp, uj_i_pot_tmp )
                {

                    x_deriv_vals.resize(4);
                    x_deriv_vals.at(0) = (jj == 0) ? &ui_j_tor_tmp : NULL;
                    x_deriv_vals.at(1) = (jj == 0) ? &ui_j_pot_tmp : NULL;
                    x_deriv_vals.at(2) = (ii == 0) ? &uj_i_tor_tmp : NULL;
                    x_deriv_vals.at(3) = (ii == 0) ? &uj_i_pot_tmp : NULL;

                    y_deriv_vals.resize(4);
                    y_deriv_vals.at(0) = (jj == 1) ? &ui_j_tor_tmp : NULL;
                    y_deriv_vals.at(1) = (jj == 1) ? &ui_j_pot_tmp : NULL;
                    y_deriv_vals.at(2) = (ii == 1) ? &uj_i_tor_tmp : NULL;
                    y_deriv_vals.at(3) = (ii == 1) ? &uj_i_pot_tmp : NULL;

                    z_deriv_vals.resize(4);
                    z_deriv_vals.at(0) = (jj == 2) ? &ui_j_tor_tmp : NULL;
                    z_deriv_vals.at(1) = (jj == 2) ? &ui_j_pot_tmp : NULL;
                    z_deriv_vals.at(2) = (ii == 2) ? &uj_i_tor_tmp : NULL;
                    z_deriv_vals.at(3) = (ii == 2) ? &uj_i_pot_tmp : NULL;


                    #pragma omp for collapse(1) schedule(static)
                    for (Ilon = 0; Ilon < Nlon; Ilon++) {

                        // ui_j
                        Cart_derivatives_at_point(
                                x_deriv_vals, y_deriv_vals, z_deriv_vals, deriv_fields,
                                latitude, longitude, Itime, Idepth, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon,
                                mask);

                        ui_j_tor.at( Ilat*Nlon + Ilon ) = ui_j_tor_tmp;
                        ui_j_pot.at( Ilat*Nlon + Ilon ) = ui_j_pot_tmp;
                        uj_i_tor.at( Ilat*Nlon + Ilon ) = uj_i_tor_tmp;
                        uj_i_pot.at( Ilat*Nlon + Ilon ) = uj_i_pot_tmp;
                    }
                }
            }

            // Outputs
            initialize_output_file( source_data, vars_to_write, output_filename.c_str(), filter_scale );

            if ( ii == 0 ) {
                write_field_to_output( u_x_tor_bar, "ui_tor_coarse", starts, counts, output_filename, &mask);
                write_field_to_output( u_x_pot_bar, "ui_pot_coarse", starts, counts, output_filename, &mask);
            } else if ( ii == 1 ) {
                write_field_to_output( u_y_tor_bar, "ui_tor_coarse", starts, counts, output_filename, &mask);
                write_field_to_output( u_y_pot_bar, "ui_pot_coarse", starts, counts, output_filename, &mask);
            } else if ( ii == 2 ) {
                write_field_to_output( u_z_tor_bar, "ui_tor_coarse", starts, counts, output_filename, &mask);
                write_field_to_output( u_z_pot_bar, "ui_pot_coarse", starts, counts, output_filename, &mask);
            }

            if ( jj == 0 ) {
                write_field_to_output( u_x_tor_bar, "uj_tor_coarse", starts, counts, output_filename, &mask);
                write_field_to_output( u_x_pot_bar, "uj_pot_coarse", starts, counts, output_filename, &mask);
            } else if ( jj == 1 ) {
                write_field_to_output( u_y_tor_bar, "uj_tor_coarse", starts, counts, output_filename, &mask);
                write_field_to_output( u_y_pot_bar, "uj_pot_coarse", starts, counts, output_filename, &mask);
            } else if ( jj == 2 ) {
                write_field_to_output( u_z_tor_bar, "uj_tor_coarse", starts, counts, output_filename, &mask);
                write_field_to_output( u_z_pot_bar, "uj_pot_coarse", starts, counts, output_filename, &mask);
            }

            write_field_to_output( ui_j_tor, "ui_j_tor", starts, counts, output_filename, &mask);
            write_field_to_output( ui_j_pot, "ui_j_pot", starts, counts, output_filename, &mask);

            write_field_to_output( uj_i_tor, "uj_i_tor", starts, counts, output_filename, &mask);
            write_field_to_output( uj_i_pot, "uj_i_pot", starts, counts, output_filename, &mask);

            write_field_to_output( uu_tor_tor_coarse, "uu_tor_tor_coarse", starts, counts, output_filename, &mask);
            write_field_to_output( uu_tor_pot_coarse, "uu_tor_pot_coarse", starts, counts, output_filename, &mask);
            write_field_to_output( uu_pot_tor_coarse, "uu_pot_tor_coarse", starts, counts, output_filename, &mask);
            write_field_to_output( uu_pot_pot_coarse, "uu_pot_pot_coarse", starts, counts, output_filename, &mask);

        }
    }

    MPI_Finalize();
    return 0;
}
