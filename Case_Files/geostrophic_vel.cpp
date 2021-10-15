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
    const std::string   &input_fname       = input.getCmdOption("--input_file",  "input.nc"),
                        &output_fname      = input.getCmdOption("--output_file", "output.nc");

    const std::string   &time_dim_name      = input.getCmdOption("--time",        "time"),
                        &depth_dim_name     = input.getCmdOption("--depth",       "depth"),
                        &latitude_dim_name  = input.getCmdOption("--latitude",    "latitude"),
                        &longitude_dim_name = input.getCmdOption("--longitude",   "longitude");

    const std::string &latlon_in_degrees  = input.getCmdOption("--is_degrees",   "true");

    const std::string &ssh_var_name = input.getCmdOption("--ssh_var",   "zos");

    const std::string   &Nprocs_in_time_string  = input.getCmdOption("--Nprocs_in_time",  "1"),
                        &Nprocs_in_depth_string = input.getCmdOption("--Nprocs_in_depth", "1");
    const int   Nprocs_in_time_input  = stoi(Nprocs_in_time_string),
                Nprocs_in_depth_input = stoi(Nprocs_in_depth_string);

    // Set OpenMP thread number
    const int max_threads = omp_get_max_threads();
    omp_set_num_threads( max_threads );

    // Print some header info, depending on debug level
    print_header_info();

    std::vector<double> longitude, latitude, time, depth;
    std::vector<double> ssh, u_lon, u_lat;
    std::vector<bool> mask;
    std::vector<int> myCounts, myStarts;

    // Read in source data / get size information
    #if DEBUG >= 1
    if (wRank == 0) { fprintf(stdout, "Reading in source data.\n\n"); }
    #endif

    dataset source_data;

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
    source_data.compute_cell_areas();

    // Read in the velocity fields
    read_var_from_file(ssh, ssh_var_name, input_fname, &mask, &myCounts, &myStarts, source_data.Nprocs_in_time, source_data.Nprocs_in_depth);
    u_lon.resize( ssh.size() );
    u_lat.resize( ssh.size() );

    std::vector<size_t> starts(myStarts.size()), counts(myCounts.size());
    for (size_t index = 0; index < starts.size(); ++index) {
        starts.at(index) = myStarts.at(index);
        counts.at(index) = myCounts.at(index);
    }

    const int   Ntime  = myCounts[0],
                Ndepth = myCounts[1],
                Nlat   = source_data.Nlat,
                Nlon   = source_data.Nlon;

    // Initialize output file
    std::vector<std::string> vars_to_write;
    vars_to_write.push_back("u_geos");
    vars_to_write.push_back("v_geos");

    initialize_output_file( source_data, vars_to_write, output_fname.c_str() );

    // Compute geostrophic velocity
    toroidal_vel_from_F( u_lon, u_lat, ssh, source_data.longitude, source_data.latitude, Ntime, Ndepth, Nlat, Nlon, mask );

    // Now mask out the equator before writing the geostrophic velocities
    int Itime, Idepth, Ilat, Ilon;
    const double Omega = 2. * M_PI / (24. * 60. * 60.);
    const double eqr_cut = 3;
    const double arg_eps = 1e-5;
    double signum, curr_lat, f_Coriolis;
    for (size_t index = 0; index < ssh.size(); ++index) {

        Index1to4(index, Itime, Idepth, Ilat, Ilon,
                         Ntime, Ndepth, Nlat, Nlon);

        curr_lat = source_data.latitude.at(Ilat);

        if ( fabs( curr_lat * 180. / M_PI ) < eqr_cut ) {
            signum = (curr_lat >= 0) - (curr_lat < 0);
            mask.at(index) = false;
        } else {
            signum = 0.;
        }

        f_Coriolis = 2 * Omega * sin( curr_lat + signum * arg_eps );
        u_lon.at(index) = u_lon.at(index) * constants::g / f_Coriolis; 
        u_lat.at(index) = u_lat.at(index) * constants::g / f_Coriolis; 

        // If bad things happen outside of the equator, throw a warning
        if ( ( ( fabs( u_lon.at(index) ) >= 1000. ) or ( fabs( u_lat.at(index) ) >= 1000. ) ) and ( signum == 0 ) ) {
            fprintf(stdout, "Bad point found!\n");
        }
    }


    // Write geostrophic velocity
    write_field_to_output( u_lon, "u_geos", &starts[0], &counts[0], output_fname.c_str(), &mask );
    write_field_to_output( u_lat, "v_geos", &starts[0], &counts[0], output_fname.c_str(), &mask );


    fprintf(stdout, "Processor %d / %d waiting to finalize.\n", wRank + 1, wSize);
    MPI_Finalize();
    return 0;
}
