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
    const std::string   &input_fname  = input.getCmdOption("--input_file",   "input.nc"),
                        &output_fname = input.getCmdOption("--output_file",  "coarse_vel.nc");

    const std::string   &time_dim_name      = input.getCmdOption("--time",        "time"),
                        &depth_dim_name     = input.getCmdOption("--depth",       "depth"),
                        &latitude_dim_name  = input.getCmdOption("--latitude",    "latitude"),
                        &longitude_dim_name = input.getCmdOption("--longitude",   "longitude");

    const std::string &latlon_in_degrees  = input.getCmdOption("--is_degrees",   "true");

    const std::string   &Nprocs_in_time_string  = input.getCmdOption("--Nprocs_in_time",  "1"),
                        &Nprocs_in_depth_string = input.getCmdOption("--Nprocs_in_depth", "1");
    const int   Nprocs_in_time_input  = stoi(Nprocs_in_time_string),
                Nprocs_in_depth_input = stoi(Nprocs_in_depth_string);

    const std::string   &zonal_vel_name    = input.getCmdOption("--zonal_vel",   "uo"),
                        &merid_vel_name    = input.getCmdOption("--merid_vel",   "vo");

    const std::string   &Nlat_reduce_factor_string = input.getCmdOption("--Nlat_reduce_factor", "1"),
                        &Nlon_reduce_factor_string = input.getCmdOption("--Nlon_reduce_factor", "1");
    const int   Nlat_reduce_factor = stoi(Nlat_reduce_factor_string),
                Nlon_reduce_factor = stoi(Nlon_reduce_factor_string);  

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
    source_data.load_time(      time_dim_name,      input_fname );
    source_data.load_depth(     depth_dim_name,     input_fname );
    source_data.load_latitude(  latitude_dim_name,  input_fname );
    source_data.load_longitude( longitude_dim_name, input_fname );

    // Apply some cleaning to the processor allotments if necessary. 
    source_data.check_processor_divisions( Nprocs_in_time_input, Nprocs_in_depth_input );
     
    // Convert to radians, if appropriate
    if ( latlon_in_degrees == "true" ) {
        convert_coordinates( source_data.longitude, source_data.latitude );
    }

    // Read in velocity fields
    source_data.load_variable( "u_lon", zonal_vel_name, input_fname, true, true );
    source_data.load_variable( "u_lat", merid_vel_name, input_fname, true, true );

    const int   full_Ntime  = source_data.full_Ntime,
                Ntime       = source_data.myCounts[0],
                Ndepth      = source_data.myCounts[1],
                Nlat        = source_data.Nlat,
                Nlon        = source_data.Nlon,
                Nlat_coarse = Nlat / Nlat_reduce_factor,
                Nlon_coarse = Nlon / Nlon_reduce_factor;

    // Now coarsen the velocity fields
    const unsigned int Npts_coarse = Ntime * Ndepth * Nlat_coarse * Nlon_coarse;
    std::vector<double> u_lon_coarse(Npts_coarse), u_lat_coarse(Npts_coarse),
                        longitude_coarse(Nlon_coarse), latitude_coarse(Nlat_coarse);
    std::vector<bool> mask_coarse(Npts_coarse);

    int Itime, Idepth, Ilat, Ilon;

    // First, the coarse grid
    for (Ilat = 0; Ilat < Nlat_coarse; ++Ilat) { latitude_coarse.at( Ilat) = source_data.latitude.at( Ilat * Nlat_reduce_factor); }
    for (Ilon = 0; Ilon < Nlon_coarse; ++Ilon) { longitude_coarse.at(Ilon) = source_data.longitude.at(Ilon * Nlon_reduce_factor); }

    // Next, the coarse velocities
    size_t II_fine, II;
    #pragma omp parallel \
    default(none) \
    shared( source_data, mask_coarse, u_lon_coarse, u_lat_coarse ) \
    private( II, II_fine, Itime, Idepth, Ilat, Ilon )
    {
        #pragma omp for collapse(1) schedule(static)
        for (II = 0; II < Npts_coarse; ++II) {

            Index1to4( II, Itime, Idepth, Ilat,        Ilon,
                           Ntime, Ndepth, Nlat_coarse, Nlon_coarse );

            II_fine = Index( Itime, Idepth, Ilat * Nlat_reduce_factor, Ilon * Nlon_reduce_factor,
                             Ntime, Ndepth, Nlat,                      Nlon );

            if ( source_data.mask.at(II_fine) ) {
                mask_coarse.at(II) = true;
                u_lon_coarse.at(II) = source_data.variables.at("u_lon").at(II_fine);
                u_lat_coarse.at(II) = source_data.variables.at("u_lat").at(II_fine);
            } else {
                mask_coarse.at(II) = false;
            }
        }
    }

    // Initialize file and write out coarsened fields
    size_t starts[4] = { source_data.myStarts.at(0), source_data.myStarts.at(1), source_data.myStarts.at(2), source_data.myStarts.at(3)};
    size_t counts[4] = { Ntime,                      Ndepth,                     Nlat_coarse,                Nlon_coarse };

    // Set up grid vectors for coarse data
    dataset coarse_data;
    coarse_data.time      = source_data.time;
    coarse_data.depth     = source_data.depth;
    coarse_data.latitude  = longitude_coarse;
    coarse_data.longitude = latitude_coarse;
    coarse_data.compute_cell_areas();

    std::vector<std::string> vars_to_write = { zonal_vel_name, merid_vel_name };
    initialize_output_file( coarse_data, vars_to_write, output_fname.c_str() );

    write_field_to_output( u_lon_coarse, zonal_vel_name, starts, counts, output_fname, &mask_coarse );
    write_field_to_output( u_lat_coarse, merid_vel_name, starts, counts, output_fname, &mask_coarse );

    add_attr_to_file( "seed_count", full_Ntime, output_fname.c_str() );

    #if DEBUG >= 1
    fprintf(stdout, "Processor %d / %d waiting to finalize.\n", wRank + 1, wSize);
    #endif
    MPI_Finalize();
    return 0;
}
