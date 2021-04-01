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
#include "../particles.hpp"

int main(int argc, char *argv[]) {
    
    // PERIODIC_Y implies UNIFORM_LAT_GRID
    static_assert( (constants::UNIFORM_LAT_GRID) or (not(constants::PERIODIC_Y)),
            "PERIODIC_Y requires UNIFORM_LAT_GRID.\n"
            "Please update constants.hpp accordingly.\n");

    static_assert( ( (constants::PERIODIC_X) and (not(constants::PERIODIC_Y)) ),
            "The particles routine currently requires globe-like periodicity.\n"
            "Please update constants.hpp accordingly.\n");


    // Specify the number of OpenMP threads
    //   and initialize the MPI world
    int thread_safety_provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &thread_safety_provided);
    //MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI::ERRORS_THROW_EXCEPTIONS);

    int wRank=-1, wSize=-1;
    MPI_Comm_rank( MPI_COMM_WORLD, &wRank );
    MPI_Comm_size( MPI_COMM_WORLD, &wSize );

    assert(wSize == 1);

    // Default input file name
    const int str_len = 50;
    char    input_fname_1[ str_len], flag_input_1[str_len],
            input_fname_2[ str_len], flag_input_2[str_len],
            output_fname[  str_len], flag_output[ str_len];

    snprintf(input_fname_1, str_len, "input_1.nc");
    snprintf(input_fname_2, str_len, "input_2.nc");
    snprintf(output_fname,  str_len, "particle_comparison.nc");

    snprintf(flag_input_1,  str_len, "--input_1");
    snprintf(flag_input_2,  str_len, "--input_2");
    snprintf(flag_output,   str_len, "--output");


    // Parse command-line flags
    int ii = 1;
    while(ii < argc) {  
        if (wRank == 0) {
            fprintf(stdout, "Argument %d : %s\n", ii, argv[ii]);
        }

        if (strcmp(argv[ii], flag_input_1) == 0) {
            snprintf(input_fname_1, 50, argv[ii+1]);
            ++ii;
        } else if (strcmp(argv[ii], flag_input_2) == 0) {
            snprintf(input_fname_2, 50, argv[ii+1]);
            ++ii;
        } else if (strcmp(argv[ii], flag_output) == 0) {
            snprintf(output_fname, 50, argv[ii+1]);
            ++ii;
        } else {
            // Otherwise, the flag is unrecognized
            if (wRank == 0) {
                fprintf(stderr, "Flag %s not recognized.\n", argv[ii]);
            }
            return -1;
        }
        ++ii;
    }

    #if DEBUG >= 0
    if (wRank == 0) {
        fprintf(stdout, "\n\n");
        fprintf(stdout, "Compiled at %s on %s.\n", __TIME__, __DATE__);
        fprintf(stdout, "  Version %d.%d.%d \n", 
                MAJOR_VERSION, MINOR_VERSION, PATCH_VERSION);
        fprintf(stdout, "\n");
        fflush(stdout);
    }
    #endif

    #if DEBUG >= 0
    if (wRank == 0) {
        if (constants::CARTESIAN) { 
            fprintf(stdout, "Using Cartesian coordinates.\n");
        } else {
            fprintf(stdout, "Using spherical coordinates.\n");
        }
        fflush(stdout);
    }
    #endif
    MPI_Barrier(MPI_COMM_WORLD);

    // Print processor assignments
    const int max_threads = omp_get_max_threads();
    omp_set_num_threads( max_threads );
    #if DEBUG >= 2
    int tid, nthreads;
    #pragma omp parallel default(none) private(tid, nthreads) \
        shared(stdout) firstprivate(wRank, wSize)
    {
        tid = omp_get_thread_num();
        nthreads = omp_get_num_threads();
        fprintf(stdout, "Hello from thread %d of %d on processor %d of %d.\n", 
                tid+1, nthreads, wRank+1, wSize);
    }
    fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD);
    #endif

    std::vector<double> time, trajectory, lon_1, lat_1, lon_2, lat_2;

    // Read in source data / get size information
    #if DEBUG >= 1
    if (wRank == 0) { fprintf(stdout, "Reading in source data.\n\n"); }
    #endif

    // Read in the grid coordinates
    read_var_from_file(time,       "time",       input_fname_1);
    read_var_from_file(trajectory, "trajectory", input_fname_1);

    const int Ntime  = time.size();
    const int Nparts = trajectory.size();

    read_var_from_file(lon_1, "longitude", input_fname_1);
    read_var_from_file(lat_1,  "latitude", input_fname_1);

    read_var_from_file(lon_2, "longitude", input_fname_2);
    read_var_from_file(lat_2,  "latitude", input_fname_2);

    assert( lon_1.size() == lat_1.size() );
    assert( lon_2.size() == lat_2.size() );

    const size_t Npts = std::min(lon_1.size(), lon_2.size());
    size_t index;

    std::vector<double> traj_dists(Npts, constants::fill_value);
    std::vector<double> out_mask(traj_dists.size());

    double lon1, lat1, lon2, lat2;

    #pragma omp parallel \
    default(none) \
    shared( lon_1, lat_1, lon_2, lat_2, traj_dists, out_mask ) \
    private( index, lon1, lat1, lon2, lat2 )
    {
        #pragma omp for collapse(1) schedule(dynamic)
        for (index = 0; index < Npts; ++index) {

            lon1 = lon_1.at(index);
            lat1 = lat_1.at(index);
            lon2 = lon_2.at(index);
            lat2 = lat_2.at(index);

            if (     ( lon1 != constants::fill_value )
                 and ( lat1 != constants::fill_value )  
                 and ( lon2 != constants::fill_value )  
                 and ( lat2 != constants::fill_value )  
               ) {
                traj_dists.at(index) = distance( lon1, lat1, lon2, lat2 );
            }

            out_mask.at(index) = traj_dists.at(index) == constants::fill_value ? 0 : 1;
        } 
    }

    //
    //// Create output file
    //


    // Open the NETCDF file
    int FLAG = NC_NETCDF4 | NC_CLOBBER | NC_MPIIO;
    int ncid=0, retval;
    char buffer [50];
    snprintf(buffer, 50, output_fname);
    retval = nc_create_par(buffer, FLAG, MPI_COMM_WORLD, MPI_INFO_NULL, &ncid);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

    // Define the dimensions
    int time_dimid, traj_dimid;
    retval = nc_def_dim(ncid, "time",       Ntime,          &time_dimid);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }
    retval = nc_def_dim(ncid, "trajectory", Nparts * wSize, &traj_dimid);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

    // Define coordinate variables
    int time_varid, traj_varid;
    retval = nc_def_var(ncid, "time",       NC_DOUBLE, 1, &time_dimid, &time_varid);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }
    retval = nc_def_var(ncid, "trajectory", NC_DOUBLE, 1, &traj_dimid, &traj_varid);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

    // Write the coordinate variables
    size_t start[1], count[1];
    start[0] = 0;
    count[0] = Ntime;
    retval = nc_put_vara_double(ncid, time_varid, start, count, &time[0]);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

    start[0] = wRank * Nparts;
    count[0] = Nparts;
    retval = nc_put_vara_double(ncid, traj_varid, start, count, &trajectory[0]);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

    // Close the file
    retval = nc_close(ncid);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

    // Add the trajectory
    const char* dim_names[] = {"time", "trajectory"};
    const int ndims = 2;
    add_var_to_file("trajectory_dists", dim_names, ndims, buffer);

    // Write
    size_t starts[2], counts[2];

    starts[0] = 0;
    counts[0] = Ntime;

    starts[1] = 0;
    counts[1] = Nparts;

    write_field_to_output(traj_dists, "trajectory_dists", starts, counts, output_fname, &out_mask);

    MPI_Finalize();
    return 0;
}
