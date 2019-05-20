#include <math.h>
#include <vector>
#include <mpi.h>
#include "../netcdf_io.hpp"
#include "../constants.hpp"

void initialize_output_file(
        const std::vector<double> & time,       /**< [in] time vector (1D) */
        const std::vector<double> & depth,      /**< [in] depth vector (1D) */
        const std::vector<double> & longitude,  /**< [in] longitude vector (1D) */
        const std::vector<double> & latitude,   /**< [in] longitude vector (1D) */
        const std::vector<double> & mask,       /**< [in] masking (land vs water, 2D) */
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
    if (( retval = nc_create_par(buffer, FLAG, comm, MPI_INFO_NULL, &ncid) ))
        NC_ERR(retval, __LINE__, __FILE__);

    retval = nc_put_att_double(ncid, NC_GLOBAL, "filter_scale", NC_FLOAT, 1, &filter_scale);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

    retval = nc_put_att_double(ncid, NC_GLOBAL, "rho0", NC_FLOAT, 1, &constants::rho0);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

    // Record coordinate type
    if (constants::CARTESIAN) {
        retval = nc_put_att_text(ncid, NC_GLOBAL, "coord-type", 10, "cartesian");
    } else {
        retval = nc_put_att_text(ncid, NC_GLOBAL, "coord-type", 10, "spherical");
    }
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

    // Extract dimension sizes
    const int Ntime   = time.size();
    const int Ndepth  = depth.size();
    const int Nlat    = latitude.size();
    const int Nlon    = longitude.size();

    // Define the dimensions
    int time_dimid, depth_dimid, lat_dimid, lon_dimid;
    if ((retval = nc_def_dim(ncid, "time",      Ntime,     &time_dimid)))
        NC_ERR(retval, __LINE__, __FILE__);
    if ((retval = nc_def_dim(ncid, "depth",     Ndepth,    &depth_dimid)))
        NC_ERR(retval, __LINE__, __FILE__);
    if ((retval = nc_def_dim(ncid, "latitude",  Nlat,      &lat_dimid)))
        NC_ERR(retval, __LINE__, __FILE__);
    if ((retval = nc_def_dim(ncid, "longitude", Nlon,      &lon_dimid)))
        NC_ERR(retval, __LINE__, __FILE__);

    // Define coordinate variables
    int time_varid, depth_varid, lat_varid, lon_varid;
    if ((retval = nc_def_var(ncid, "time",      NC_FLOAT, 1, &time_dimid,  &time_varid)))
        NC_ERR(retval, __LINE__, __FILE__);
    if ((retval = nc_def_var(ncid, "depth",     NC_FLOAT, 1, &depth_dimid, &depth_varid)))
        NC_ERR(retval, __LINE__, __FILE__);
    if ((retval = nc_def_var(ncid, "latitude",  NC_FLOAT, 1, &lat_dimid,   &lat_varid)))
        NC_ERR(retval, __LINE__, __FILE__);
    if ((retval = nc_def_var(ncid, "longitude", NC_FLOAT, 1, &lon_dimid,   &lon_varid)))
        NC_ERR(retval, __LINE__, __FILE__);

    if (not(constants::CARTESIAN)) {
        const double rad_to_degree = 180. / M_PI;
        retval = nc_put_att_double(ncid, lon_varid, "scale_factor", NC_FLOAT, 1, &rad_to_degree);
        if (retval) { NC_ERR(retval, __LINE__, __FILE__); }
        retval = nc_put_att_double(ncid, lat_varid, "scale_factor", NC_FLOAT, 1, &rad_to_degree);
        if (retval) { NC_ERR(retval, __LINE__, __FILE__); }
    }

    int mask_dimids[2];
    mask_dimids[0] = lat_dimid;
    mask_dimids[1] = lon_dimid;
    int mask_varid;
    if ((retval = nc_def_var(ncid, "mask", NC_FLOAT, 2, mask_dimids, &mask_varid)))
        NC_ERR(retval, __LINE__, __FILE__);

    // Write the coordinate variables
    size_t start[1], count[1];
    start[0] = 0;
    count[0] = Ntime;
    if ((retval = nc_put_vara_double(ncid, time_varid,  start, count, &time[0])))
        NC_ERR(retval, __LINE__, __FILE__);

    count[0] = Ndepth;
    if ((retval = nc_put_vara_double(ncid, depth_varid, start, count, &depth[0])))
        NC_ERR(retval, __LINE__, __FILE__);

    count[0] = Nlat;
    if ((retval = nc_put_vara_double(ncid, lat_varid,   start, count, &latitude[0])))
        NC_ERR(retval, __LINE__, __FILE__);

    count[0] = Nlon;
    if ((retval = nc_put_vara_double(ncid, lon_varid,   start, count, &longitude[0])))
        NC_ERR(retval, __LINE__, __FILE__);

    size_t mask_start[2], mask_count[2];
    mask_start[0] = 0;
    mask_start[1] = 0;
    mask_count[0] = Nlat;
    mask_count[1] = Nlon;
    if ((retval = nc_put_vara_double(ncid, mask_varid, mask_start, mask_count, &mask[0])))
        NC_ERR(retval, __LINE__, __FILE__);

    // Close the file
    if ((retval = nc_close(ncid))) { NC_ERR(retval, __LINE__, __FILE__); }

    #if DEBUG >= 2
    fprintf(stdout, "Output file (%s) initialized.\n\n", buffer);
    #endif

    if (wRank == 0) {
        // Loop through and add the desired variables
        // Dimension names (in order!)
        const char* dim_names[] = {"time", "depth", "latitude", "longitude"};
        const int ndims = 4;
        for (size_t varInd = 0; varInd < vars.size(); ++varInd) {
            add_var_to_file(vars.at(varInd), dim_names, ndims, buffer);
        }
    }
}
