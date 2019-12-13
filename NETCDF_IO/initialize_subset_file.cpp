#include <math.h>
#include <vector>
#include <mpi.h>
#include "../netcdf_io.hpp"
#include "../constants.hpp"

void initialize_subset_file(
        const std::vector<double> & time,       /**< [in] time vector (1D) */
        const std::vector<double> & depth,      /**< [in] longitude vector (1D) */
        const std::vector<double> & windows,    /**< [in] windowing sizes */
        const int & Nsamples,                   /**< [in] number of samples */
        const std::vector<std::string> & vars,  /**< [in] name of variables to write */
        const char * filename,                  /**< [in] name for the output file */
        const double filter_scale,              /**< [in] lengthscale used in the filter */
        const MPI_Comm comm                     /**< [in] MPI Communicator */
        ) {

    int wRank=-1, wSize=-1;
    MPI_Comm_rank( MPI_COMM_WORLD, &wRank );
    MPI_Comm_size( MPI_COMM_WORLD, &wSize );

    // Open the NETCDF file
    int FLAG = NC_NETCDF4 | NC_CLOBBER | NC_MPIIO;
    int ncid=0, retval;
    char buffer [50];
    snprintf(buffer, 50, filename);
    retval = nc_create_par(buffer, FLAG, comm, MPI_INFO_NULL, &ncid);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

    retval = nc_put_att_double(ncid, NC_GLOBAL, "filter_scale", NC_FLOAT, 1, &filter_scale);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

    retval = nc_put_att_double(ncid, NC_GLOBAL, "rho0", NC_FLOAT, 1, &constants::rho0);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

    // Extract dimension sizes
    const int Ntime    = time.size();
    const int Ndepth   = depth.size();
    const int Nwindows = windows.size();

    // Define the dimensions
    int time_dimid, depth_dimid, wind_dimid, samp_dimid;
    retval = nc_def_dim(ncid, "time",   Ntime,    &time_dimid);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }
    retval = nc_def_dim(ncid, "depth",  Ndepth,   &depth_dimid);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }
    retval = nc_def_dim(ncid, "window", Nwindows, &wind_dimid);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }
    retval = nc_def_dim(ncid, "sample", Nsamples, &samp_dimid);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

    // Define coordinate variables
    int time_varid, depth_varid, wind_varid, samp_varid;
    retval = nc_def_var(ncid, "time",   NC_FLOAT, 1, &time_dimid,  &time_varid);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }
    retval = nc_def_var(ncid, "depth",  NC_FLOAT, 1, &depth_dimid, &depth_varid);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }
    retval = nc_def_var(ncid, "window", NC_FLOAT, 1, &wind_dimid,  &wind_varid);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }
    retval = nc_def_var(ncid, "sample", NC_FLOAT, 1, &samp_dimid,  &samp_varid);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

    // Write the coordinate variables
    size_t start[1], count[1];
    start[0] = 0;
    count[0] = Ntime;
    retval = nc_put_vara_double(ncid, time_varid,  start, count, &time[0]);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

    count[0] = Ndepth;
    retval = nc_put_vara_double(ncid, depth_varid, start, count, &depth[0]);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

    count[0] = Nwindows;
    retval = nc_put_vara_double(ncid, wind_varid, start, count, &windows[0]);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

    // Close the file
    retval = nc_close(ncid);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

    #if DEBUG >= 2
    fprintf(stdout, "Output file (%s) initialized.\n\n", buffer);
    #endif

    if (wRank == 0) {
        // Loop through and add the desired variables
        // Dimension names (in order!)
        const char* dim_names[] = {"time", "depth", "window", "sample"};
        const int ndims = 4;
        for (size_t varInd = 0; varInd < vars.size(); ++varInd) {
            add_var_to_file(vars.at(varInd), dim_names, ndims, buffer);
        }
    }
}
