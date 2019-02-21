#include "../netcdf_io.hpp"

void initialize_output_file(
        const int Ntime, const int Ndepth, const int Nlon, 
                         const int Nlat, const int Nscales,
        const double * time,      const double * depth, 
        const double * longitude, const double * latitude, 
        const double * scales) {

    // Open the NETCDF file
    int FLAG = NC_NETCDF4 | NC_CLOBBER;
    int ncid=0, retval;
    char buffer [50];
    snprintf(buffer, 50, "filter_output.nc");
    if (( retval = nc_create(buffer, FLAG, &ncid) ))
        NC_ERR(retval, __LINE__, __FILE__);

    // Define the dimensions
    int scale_dimid, time_dimid, depth_dimid, lat_dimid, lon_dimid;
    if ((retval = nc_def_dim(ncid, "scale",     Nscales+1, &scale_dimid)))
        NC_ERR(retval, __LINE__, __FILE__);
    if ((retval = nc_def_dim(ncid, "time",      Ntime,     &time_dimid)))
        NC_ERR(retval, __LINE__, __FILE__);
    if ((retval = nc_def_dim(ncid, "depth",     Ndepth,    &depth_dimid)))
        NC_ERR(retval, __LINE__, __FILE__);
    if ((retval = nc_def_dim(ncid, "latitude",  Nlat,      &lat_dimid)))
        NC_ERR(retval, __LINE__, __FILE__);
    if ((retval = nc_def_dim(ncid, "longitude", Nlon,      &lon_dimid)))
        NC_ERR(retval, __LINE__, __FILE__);

    // Define coordinate variables
    int scale_varid, time_varid, depth_varid, lat_varid, lon_varid;
    if ((retval = nc_def_var(ncid, "scale",     NC_DOUBLE, 1, &scale_dimid, &scale_varid)))
        NC_ERR(retval, __LINE__, __FILE__);
    if ((retval = nc_def_var(ncid, "time",      NC_DOUBLE, 1, &time_dimid,  &time_varid)))
        NC_ERR(retval, __LINE__, __FILE__);
    if ((retval = nc_def_var(ncid, "depth",     NC_DOUBLE, 1, &depth_dimid, &depth_varid)))
        NC_ERR(retval, __LINE__, __FILE__);
    if ((retval = nc_def_var(ncid, "latitude",  NC_DOUBLE, 1, &lat_dimid,   &lat_varid)))
        NC_ERR(retval, __LINE__, __FILE__);
    if ((retval = nc_def_var(ncid, "longitude", NC_DOUBLE, 1, &lon_dimid,   &lon_varid)))
        NC_ERR(retval, __LINE__, __FILE__);

    // CF ordering (with scale added at the beginning)
    const int ndims = 5;
    int dimids[ndims];
    dimids[0] = scale_dimid;
    dimids[1] = time_dimid;
    dimids[2] = depth_dimid;
    dimids[3] = lat_dimid;
    dimids[4] = lon_dimid;

    // Declare variables
    int u_r_varid, u_lon_varid, u_lat_varid;
    if ((retval = nc_def_var(ncid, "u_r",   NC_DOUBLE, ndims, dimids, &u_r_varid)))
        NC_ERR(retval, __LINE__, __FILE__);
    if ((retval = nc_def_var(ncid, "u_lon", NC_DOUBLE, ndims, dimids, &u_lon_varid)))
        NC_ERR(retval, __LINE__, __FILE__);
    if ((retval = nc_def_var(ncid, "u_lat", NC_DOUBLE, ndims, dimids, &u_lat_varid)))
        NC_ERR(retval, __LINE__, __FILE__);

    // Write the coordinate variables
    size_t start[1], count[1];
    start[0] = 0;
    count[0] = Ntime;
    if ((retval = nc_put_vara_double(ncid, time_varid,  start, count, time)))
        NC_ERR(retval, __LINE__, __FILE__);

    count[0] = Ndepth;
    if ((retval = nc_put_vara_double(ncid, depth_varid, start, count, depth)))
        NC_ERR(retval, __LINE__, __FILE__);

    count[0] = Nlat;
    if ((retval = nc_put_vara_double(ncid, lat_varid,   start, count, latitude)))
        NC_ERR(retval, __LINE__, __FILE__);

    count[0] = Nlon;
    if ((retval = nc_put_vara_double(ncid, lon_varid,   start, count, longitude)))
        NC_ERR(retval, __LINE__, __FILE__);

    count[0] = Nscales;
    if ((retval = nc_put_vara_double(ncid, scale_varid, start, count, scales)))
        NC_ERR(retval, __LINE__, __FILE__);

    // Close the file
    if ((retval = nc_close(ncid))) { NC_ERR(retval, __LINE__, __FILE__); }

}
