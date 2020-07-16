#include <math.h>
#include <vector>
#include <mpi.h>
#include "../netcdf_io.hpp"
#include "../constants.hpp"

void initialize_particle_file(
        const std::vector<double> & time,
        const std::vector<double> & trajectory,
        std::vector<std::string> & vars,
        const std::string & filename,
        const MPI_Comm comm
        ) {

    int wRank=-1, wSize=-1;
    MPI_Comm_rank( MPI_COMM_WORLD, &wRank );
    MPI_Comm_size( MPI_COMM_WORLD, &wSize );

    // Open the NETCDF file
    int FLAG = NC_NETCDF4 | NC_CLOBBER | NC_MPIIO;
    int ncid=0, retval;
    char buffer [50];
    snprintf(buffer, 50, filename.c_str());
    retval = nc_create_par(buffer, FLAG, comm, MPI_INFO_NULL, &ncid);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

    // Record coordinate type
    if (constants::CARTESIAN) {
        retval = nc_put_att_text(ncid, NC_GLOBAL, "coord-type", 10, "cartesian");
    } else {
        retval = nc_put_att_text(ncid, NC_GLOBAL, "coord-type", 10, "spherical");
    }
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

    // Extract dimension sizes
    const int Ntime  = time.size();
    const int Nparts = trajectory.size();

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

    #if DEBUG >= 2
    if (wRank == 0) { fprintf(stdout, "\nOutput file (%s) initialized.\n", buffer); }
    #endif

    vars.push_back("longitude");
    vars.push_back("latitude");

    #if DEBUG >= 1
    vars.push_back("rev_longitude");
    vars.push_back("rev_latitude");
    vars.push_back("fore_back_dists");
    #endif

    if (wRank == 0) {
        // Loop through and add the desired variables
        // Dimension names (in order!)
        const char* dim_names[] = {"time", "trajectory"};
        const int ndims = 2;
        for (size_t varInd = 0; varInd < vars.size(); ++varInd) {
            add_var_to_file(vars.at(varInd), dim_names, ndims, buffer);
        }
    }

    #if DEBUG >= 2
    if (wRank == 0) { fprintf(stdout, "\n"); }
    #endif
}
