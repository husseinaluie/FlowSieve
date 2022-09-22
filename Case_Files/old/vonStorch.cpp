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
#include "../postprocess.hpp"

int main(int argc, char *argv[]) {
    
    static_assert ( not(constants::FILTER_OVER_LAND), "Cannot have FILTER_OVER_LAND on when computing vonStorch" );

    static_assert( not( constants::DO_OKUBOWEISS_ANALYSIS ), "No OkuboWeiss available to pass to post-processing." );

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
                        &output_fname       = input.getCmdOption("--output_file",   "vonStorch.nc"),
                        &postproc_fname     = input.getCmdOption("--postproc_file", "postprocess");

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
                        &uu_name    = input.getCmdOption("--uu",    "uu"),
                        &uv_name    = input.getCmdOption("--uv",    "uv"),
                        &vv_name    = input.getCmdOption("--vv",    "vv");

    const std::string   &region_defs_fname    = input.getCmdOption("--region_definitions_file",    "region_definitions.nc"),
                        &region_defs_dim_name = input.getCmdOption("--region_definitions_dim",     "region"),
                        &region_defs_var_name = input.getCmdOption("--region_definitions_var",     "region_definition");

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

    source_data.load_variable( "uu", uu_name, input_fname, true, true );
    source_data.load_variable( "uv", uv_name, input_fname, true, true );
    source_data.load_variable( "vv", vv_name, input_fname, true, true );

    const std::vector<double>   &latitude  = source_data.latitude,
                                &longitude = source_data.longitude;

    const std::vector<double>   &u_lon  = source_data.variables.at("zonal_vel"),
                                &u_lat  = source_data.variables.at("merid_vel"),
                                &uu = source_data.variables.at("uu"),
                                &uv = source_data.variables.at("uv"),
                                &vu = source_data.variables.at("uv"),
                                &vv = source_data.variables.at("vv");

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
    std::vector<double> vonStorch( source_data.variables.at("zonal_vel").size(), 0. );
    std::vector<const std::vector<double>*> deriv_fields { &u_lon, &u_lat };

    double ulon_lon, ulon_lat, ulat_lon, ulat_lat;
    std::vector<double*> lon_deriv_vals, lat_deriv_vals;

    const size_t Npts = Ntime * Ndepth * Nlat * Nlon;
    size_t index;
    int Itime, Idepth, Ilat, Ilon;

    const int OMP_chunksize = get_omp_chunksize(Nlat,Nlon);

    double u_loc, v_loc, tau_uu, tau_uv, tau_vu, tau_vv;
    #pragma omp parallel \
    default(none) \
    shared( source_data, latitude, longitude, deriv_fields, u_lon, u_lat, uu, uv, vu, vv, vonStorch )\
    private( Itime, Idepth, Ilat, Ilon, index, ulon_lon, ulon_lat, ulat_lon, ulat_lat, lon_deriv_vals, lat_deriv_vals,\
             u_loc, v_loc, tau_uu, tau_uv, tau_vu, tau_vv )
    {

        lon_deriv_vals.push_back( &ulon_lon );
        lon_deriv_vals.push_back( &ulat_lon );

        lat_deriv_vals.push_back( &ulon_lat );
        lat_deriv_vals.push_back( &ulat_lat );

        #pragma omp for collapse(1) schedule(dynamic, OMP_chunksize)
        for ( index = 0; index < Npts; index++ ) {

            if (source_data.mask.at(index)) {
                Index1to4( index, Itime, Idepth, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon );

                spher_derivative_at_point( lat_deriv_vals, deriv_fields, latitude,  "lat", Itime, Idepth, Ilat, Ilon, 
                                                                                           Ntime, Ndepth, Nlat, Nlon, source_data.mask );
                spher_derivative_at_point( lon_deriv_vals, deriv_fields, longitude, "lon", Itime, Idepth, Ilat, Ilon, 
                                                                                           Ntime, Ndepth, Nlat, Nlon, source_data.mask );

                // rho0 * [ tau( u, vec(u) ) dot grad( u )  +  tau( v, vec(u) ) dot grad ( v ) ]
                //      where tau( a, b ) = bar(ab) - bar(a) bar(b)
                //      and bar(.) is a time mean
                
                u_loc = u_lon.at(index);
                v_loc = u_lat.at(index);
                tau_uu = uu.at(index) - u_loc * u_loc;
                tau_uv = uv.at(index) - u_loc * v_loc;
                tau_vu = vu.at(index) - v_loc * u_loc;
                tau_vv = vv.at(index) - v_loc * v_loc;
                      
                vonStorch.at(index) = ( constants::rho0 / constants::R_earth ) * (
                              tau_uu * ulon_lon / cos(latitude.at(Ilat))
                            + tau_uv * ulon_lat 
                            + tau_vu * ulat_lon / cos(latitude.at(Ilat))
                            + tau_vv * ulat_lat
                        );
            }
        }
    }

    //
    //// Write the output
    //
    const int ndims = 4;
    size_t starts[ndims] = { 0,     0,      0,    0    };
    size_t counts[ndims] = { Ntime, Ndepth, Nlat, Nlon };

    std::vector<std::string> vars_to_write;
    vars_to_write.push_back("C_Km_Ke");

    if ( not( constants::NO_FULL_OUTPUTS ) ) {
        initialize_output_file( source_data, vars_to_write, output_fname.c_str(), -1);
        write_field_to_output( vonStorch, "C_Km_Ke", starts, counts, output_fname.c_str(), &(source_data.mask) );
    }

    //
    //// Also do postprocess for region summing
    //

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

    std::vector<const std::vector<double>*> postprocess_fields;
    std::vector<std::string> postprocess_names;

    postprocess_names.push_back( "C_Km_Ke" );
    postprocess_fields.push_back( &vonStorch );

    std::vector<double> OkuboWeiss;
    Apply_Postprocess_Routines( source_data, postprocess_fields, postprocess_names, OkuboWeiss, -1, postproc_fname );

    // Done!
    MPI_Finalize();
    return 0;
}
