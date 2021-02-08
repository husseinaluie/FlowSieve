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
    const std::string &input_fname       = input.getCmdOption("--input_file",  "input.nc");
    const std::string &output_fname      = input.getCmdOption("--output_file", "output.nc");

    const std::string &time_dim_name      = input.getCmdOption("--time",        "time");
    const std::string &depth_dim_name     = input.getCmdOption("--depth",       "depth");
    const std::string &latitude_dim_name  = input.getCmdOption("--latitude",    "latitude");
    const std::string &longitude_dim_name = input.getCmdOption("--longitude",   "longitude");

    const std::string &ssh_vel_name = input.getCmdOption("--ssh_vel",   "zos");

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

    // Read in the grid coordinates
    read_var_from_file(time,      time_dim_name,      input_fname);
    read_var_from_file(depth,     depth_dim_name,     input_fname);
    read_var_from_file(latitude,  latitude_dim_name,  input_fname);
    read_var_from_file(longitude, longitude_dim_name, input_fname);
     
    convert_coordinates(longitude, latitude);

    const int Ntime  = time.size();
    const int Ndepth = depth.size();
    const int Nlon   = longitude.size();
    const int Nlat   = latitude.size();

    std::vector<double> areas(Nlon * Nlat);
    compute_areas(areas, longitude, latitude);

    // Read in the velocity fields
    read_var_from_file(ssh, ssh_vel_name, input_fname, &mask, &myCounts, &myStarts);
    u_lon.resize( ssh.size() );
    u_lat.resize( ssh.size() );

    std::vector<size_t> starts(myStarts.size()), counts(myCounts.size());
    for (size_t index = 0; index < starts.size(); ++index) {
        starts.at(index) = myStarts.at(index);
        counts.at(index) = myCounts.at(index);
    }

    // Initialize output file
    std::vector<std::string> vars_to_write;
    vars_to_write.push_back("u_geos");
    vars_to_write.push_back("v_geos");

    initialize_output_file(
        time, depth, longitude, latitude, areas,
        vars_to_write, output_fname.c_str(), 0.
        );

    // Compute geostrophic velocity
    toroidal_vel_from_F( u_lon, u_lat, ssh, 
            longitude, latitude,
            Ntime, Ndepth, Nlat, Nlon, mask);

    // Now mask out the equator before writing the geostrophic velocities
    int Itime, Idepth, Ilat, Ilon;
    const double Omega = 2. * M_PI / (24. * 60. * 60.);
    const double eqr_cut = 3;
    const double arg_eps = 1e-5;
    double signum, curr_lat, f_Coriolis;
    for (size_t index = 0; index < ssh.size(); ++index) {

        Index1to4(index, Itime, Idepth, Ilat, Ilon,
                         Ntime, Ndepth, Nlat, Nlon);

        curr_lat = latitude.at(Ilat);

        if ( fabs( curr_lat * 180. / M_PI ) < eqr_cut ) {
            signum = (curr_lat >= 0) - (curr_lat < 0);
            mask.at(index) = false;
        } else {
            signum = 0.;
        }

        f_Coriolis = 2 * Omega * sin( curr_lat + signum * arg_eps );
        u_lon.at(index) = u_lon.at(index) * constants::g / f_Coriolis; 
        u_lat.at(index) = u_lat.at(index) * constants::g / f_Coriolis; 

        // Also mask out garbage
        if ( ( fabs( u_lon.at(index) ) >= 1000. ) or ( fabs( u_lat.at(index) ) >= 1000. ) ) {
            //mask.at(index) = false;
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
