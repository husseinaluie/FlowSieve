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
#include "../differentiation_tools.hpp"
#include "../constants.hpp"

int main(int argc, char *argv[]) {
    
    static_assert ( not(constants::FILTER_OVER_LAND), "Cannot have FILTER_OVER_LAND on when computing vonStorch" );

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
    assert(wSize == 1);

    //
    //// Parse command-line arguments
    //
    InputParser input(argc, argv);
    if(input.cmdOptionExists("--version")){
        if (wRank == 0) { print_compile_info(NULL); } 
        return 0;
    }

    // first argument is the flag, second argument is default value (for when flag is not present)
    const std::string   &input_fname        = input.getCmdOption("--input_file",    "input.nc"),
                        &output_fname       = input.getCmdOption("--output_file",   "vonStorch.nc");

    const std::string   &time_dim_name      = input.getCmdOption("--time",        "time"),
                        &depth_dim_name     = input.getCmdOption("--depth",       "depth"),
                        &latitude_dim_name  = input.getCmdOption("--latitude",    "latitude"),
                        &longitude_dim_name = input.getCmdOption("--longitude",   "longitude");

    const std::string &latlon_in_degrees  = input.getCmdOption("--is_degrees",   "true");

    const std::string   &Nprocs_in_time_string  = input.getCmdOption("--Nprocs_in_time",  "1"),
                        &Nprocs_in_depth_string = input.getCmdOption("--Nprocs_in_depth", "1");
    const int   Nprocs_in_time_input  = stoi(Nprocs_in_time_string),
                Nprocs_in_depth_input = stoi(Nprocs_in_depth_string);

    const std::string   &zonal_vel_name = input.getCmdOption("--zonal_vel", "uo"),
                        &merid_vel_name = input.getCmdOption("--merid_vel", "vo"),
                        &tau_uu_name    = input.getCmdOption("--tau_uu",    "tau_uu"),
                        &tau_uv_name    = input.getCmdOption("--tau_uv",    "tau_uv"),
                        &tau_vv_name    = input.getCmdOption("--tau_vv",    "tau_vv");

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
    source_data.load_time(      time_dim_name,      input_fname );
    source_data.load_depth(     depth_dim_name,     input_fname );
    source_data.load_latitude(  latitude_dim_name,  input_fname );
    source_data.load_longitude( longitude_dim_name, input_fname );

    // Apply some cleaning to the processor allotments if necessary. 
    source_data.check_processor_divisions( Nprocs_in_time_input, Nprocs_in_depth_input );
     
    // Convert to radians, if appropriate
    if ( latlon_in_degrees == "true" ) { convert_coordinates( source_data.longitude, source_data.latitude ); }

    // Compute the area of each 'cell' which will be necessary for integration
    source_data.compute_cell_areas();

    // Read in the toroidal and potential fields
    source_data.load_variable( "zonal_vel", zonal_vel_name, input_fname, true, true );
    source_data.load_variable( "merid_vel", merid_vel_name, input_fname, true, true );

    source_data.load_variable( "tau_uu", tau_uu_name, input_fname, false, true );
    source_data.load_variable( "tau_uv", tau_uv_name, input_fname, false, true );
    source_data.load_variable( "tau_vv", tau_vv_name, input_fname, false, true );

    const std::vector<double>   &latitude  = source_data.latitude,
                                &longitude = source_data.longitude;

    const std::vector<double>   &u_lon  = source_data.variables.at("zonal_vel"),
                                &u_lat  = source_data.variables.at("merid_vel"),
                                &tau_uu = source_data.variables.at("tau_uu"),
                                &tau_uv = source_data.variables.at("tau_uv"),
                                &tau_vu = source_data.variables.at("tau_uv"),
                                &tau_vv = source_data.variables.at("tau_vv");

    // Get the MPI-local dimension sizes
    source_data.Ntime  = source_data.myCounts[0];
    source_data.Ndepth = source_data.myCounts[1];

    //
    const int   Ntime   = source_data.Ntime,
                Ndepth  = source_data.Ndepth,
                Nlat    = source_data.Nlat,
                Nlon    = source_data.Nlon;

    // Mask out the pole, if necessary (i.e. set lat = 90 to land)
    mask_out_pole( source_data.latitude, source_data.mask, Ntime, Ndepth, Nlat, Nlon );

    // Compute vonStorch point-wise
    std::vector<double> vonStorch( source_data.variables.at("zonal_vel").size() );
    std::vector<const std::vector<double>*> deriv_fields { &u_lon, &u_lat };

    double ulon_lon, ulon_lat, ulat_lon, ulat_lat;
    std::vector<double*>    lon_deriv_vals { &ulon_lon, &ulat_lon },
                            lat_deriv_vals { &ulon_lat, &ulat_lat };

    for (int Idepth = 0; Idepth < source_data.Ndepth; Idepth++) {
        for (int Ilat = 0; Ilat < source_data.Nlat; Ilat++) {
            for (int Ilon = 0; Ilon < source_data.Nlon; Ilon++) {

                size_t index = Index( 0, Idepth, Ilat, Ilon, 1, Ndepth, Nlat, Nlon);

                spher_derivative_at_point( lat_deriv_vals, deriv_fields, latitude,  "lat", 0, Idepth, Ilat, Ilon, 1, Ndepth, Nlat, Nlon, source_data.mask );
                spher_derivative_at_point( lon_deriv_vals, deriv_fields, longitude, "lon", 0, Idepth, Ilat, Ilon, 1, Ndepth, Nlat, Nlon, source_data.mask );

                // rho0 * [ tau( u, vec(u) ) dot grad( u )  +  tau( v, vec(u) ) dot grad ( v ) ]
                //      where tau( a, b ) = bar(ab) - bar(a) bar(b)
                //      and bar(.) is a time mean
                vonStorch.at(index) = ( constants::rho0 / constants::R_earth ) * (
                              tau_uu.at(index) * ulon_lon / cos(latitude.at(Ilat))
                            + tau_uv.at(index) * ulon_lat 
                            + tau_vu.at(index) * ulat_lon / cos(latitude.at(Ilat))
                            + tau_vv.at(index) * ulat_lat
                        );
            }
        }
    }

    //
    //// Write the output
    //
    const int ndims = 4;
    size_t starts[ndims] = { 0, 0,              0,            0            };
    size_t counts[ndims] = { 1, size_t(Ndepth), size_t(Nlat), size_t(Nlon) };

    std::vector<std::string> vars_to_write;
    vars_to_write.push_back("C_Km_Ke");

    initialize_output_file( source_data, vars_to_write, output_fname.c_str(), -1);
    write_field_to_output( vonStorch, "C_Km_Ke", starts, counts, output_fname.c_str(), &(source_data.mask) );

    // Done!
    MPI_Finalize();
    return 0;
}
