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

    // Currently cannot be Cartesian
    static_assert( not(constants::CARTESIAN),
            "Toroidal projection not set to handle Cartesian coordinates.\n"
            );

    // Enable all floating point exceptions but FE_INEXACT
    //feenableexcept(FE_ALL_EXCEPT & ~FE_INEXACT);

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
    const std::string   &input_fname        = input.getCmdOption("--input_file",    "input.nc"),
                        &output_fname       = input.getCmdOption("--output_file",   "projection_Helmholtz_uiuj.nc"),
                        &seed_fname         = input.getCmdOption("--seed_file",     "seed.nc");

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
                        &merid_vel_name    = input.getCmdOption("--merid_vel",   "vo"),
                        &v_r_seed_name     = input.getCmdOption("--v_r_seed",    "v_r_seed"),
                        &Phi_v_seed_name   = input.getCmdOption("--Phi_v_seed",  "Phi_v_seed"),
                        &Psi_v_seed_name   = input.getCmdOption("--Psi_v_seed",  "Psi_v_seed");

    const std::string &Tikhov_Lap_string = input.getCmdOption("--Tikhov_Laplace", "1.");
    const double Tikhov_Laplace = stod(Tikhov_Lap_string);  

    const std::string &Tikhov_string = input.getCmdOption("--Tikhov_Lambda", "0");
    const double Tikhov_Lambda = stod(Tikhov_string);  

    const std::string &tolerance_string = input.getCmdOption("--tolerance", "5e-3");
    const double tolerance = stod(tolerance_string);  

    const std::string &iteration_string = input.getCmdOption("--max_iterations", "100000");
    const int max_iterations = stod(iteration_string);  

    const std::string &use_mask_string = input.getCmdOption("--use_mask", "false");
    const bool use_mask = string_to_bool(use_mask_string);

    const std::string &use_area_weight_string = input.getCmdOption("--use_area_weight", "true");
    const bool use_area_weight = string_to_bool(use_area_weight_string);

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
    if ( (latlon_in_degrees == "true") and (not(constants::CARTESIAN)) ) {
        convert_coordinates( source_data.longitude, source_data.latitude );
    }

    // Compute the area of each 'cell' which will be necessary for integration
    source_data.compute_cell_areas();

    // Read in the velocity fields
    source_data.load_variable( "u_lon", zonal_vel_name, input_fname, true, true );
    source_data.load_variable( "u_lat", merid_vel_name, input_fname, true, true );

    // Mask out the pole, if necessary (i.e. set lat = 90 to land)
    mask_out_pole( source_data.latitude, source_data.mask, source_data.Ntime, source_data.Ndepth, source_data.Nlat, source_data.Nlon );

    // Read in the seed
    double seed_count;
    bool single_seed;
    std::vector<double> seed_v_r, seed_Phi_v, seed_Psi_v;
    if (seed_fname == "zero") {
        seed_count = 1.;
        seed_v_r.resize(   source_data.Nlat * source_data.Nlon, 0. );
        seed_Phi_v.resize( source_data.Nlat * source_data.Nlon, 0. );
        seed_Psi_v.resize( source_data.Nlat * source_data.Nlon, 0. );
    } else {
        read_attr_from_file(seed_count, "seed_count", seed_fname);
        const int Nprocs_in_time  = source_data.Nprocs_in_time,
                  Nprocs_in_depth = source_data.Nprocs_in_depth;
        read_var_from_file( seed_v_r,   v_r_seed_name,   seed_fname, NULL, NULL, NULL, Nprocs_in_time, Nprocs_in_depth, not(single_seed) );
        read_var_from_file( seed_Phi_v, Phi_v_seed_name, seed_fname, NULL, NULL, NULL, Nprocs_in_time, Nprocs_in_depth, not(single_seed) );
        read_var_from_file( seed_Psi_v, Psi_v_seed_name, seed_fname, NULL, NULL, NULL, Nprocs_in_time, Nprocs_in_depth, not(single_seed) );
    }
    single_seed = (seed_count == 1);

    // Apply to projection routine
    Apply_Helmholtz_Projection_uiuj( output_fname, source_data, seed_v_r, seed_Phi_v, seed_Psi_v, 
            single_seed, 3 * source_data.Nlat * source_data.Nlon * Tikhov_Lambda, Tikhov_Laplace, 
            tolerance, max_iterations, use_area_weight, use_mask );

    // Done!
    #if DEBUG >= 0
    if (wRank == 0) {
        fprintf(stdout, "\n\n");
        fprintf(stdout, "Process completed.\n");
        fprintf(stdout, "\n");
    }
    #endif

    #if DEBUG >= 1
    fprintf(stdout, "Processor %d / %d waiting to finalize.\n", wRank + 1, wSize);
    #endif 
    MPI_Finalize();
    return 0;
}
