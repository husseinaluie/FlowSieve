#include "../netcdf_io.hpp"

#ifndef DEBUG
    #define DEBUG false
#endif

/*
 * Variables to read: uo, vo
 *
 * Dimensions: time, depth, longitude, latitude
 *
 */


// Write to netcdf file
void read_source(
        double ** longitude, double ** latitude,
        double ** time,      double ** depth,
        double ** u_lon,     double ** u_lat) {

    const bool debug = DEBUG;

    // Open the NETCDF file
    int FLAG = NC_NETCDF4;
    int ncid=0, retval;
    char buffer [50];
    snprintf(buffer, 50, "input.nc");
    if (( retval = nc_open(buffer, FLAG, &ncid) )) { NC_ERR(retval, __LINE__, __FILE__); }
    
    // Define the dimensions
    int lon_dimid, lat_dimid, time_dimid, depth_dimid;
    if ((retval = nc_inq_dimid(ncid, "time",      &time_dimid  ))) { NC_ERR(retval, __LINE__, __FILE__); }
    if ((retval = nc_inq_dimid(ncid, "depth",     &depth_dimid ))) { NC_ERR(retval, __LINE__, __FILE__); }
    if ((retval = nc_inq_dimid(ncid, "latitude",  &lat_dimid   ))) { NC_ERR(retval, __LINE__, __FILE__); }
    if ((retval = nc_inq_dimid(ncid, "longitude", &lon_dimid   ))) { NC_ERR(retval, __LINE__, __FILE__); }

    size_t Nlon, Nlat, Ntime, Ndepth;
    if ((retval = nc_inq_dim(ncid, time_dimid , NULL, &Ntime  ))) { NC_ERR(retval, __LINE__, __FILE__); }
    if ((retval = nc_inq_dim(ncid, depth_dimid, NULL, &Ndepth ))) { NC_ERR(retval, __LINE__, __FILE__); }
    if ((retval = nc_inq_dim(ncid, lat_dimid  , NULL, &Nlon   ))) { NC_ERR(retval, __LINE__, __FILE__); }
    if ((retval = nc_inq_dim(ncid, lon_dimid  , NULL, &Nlat   ))) { NC_ERR(retval, __LINE__, __FILE__); }

    if (debug) {
        fprintf(stdout, "\n");
        fprintf(stdout, "Nlon   = %zu\n", Nlon);
        fprintf(stdout, "Nlat   = %zu\n", Nlat);
        fprintf(stdout, "Ntime  = %zu\n", Ntime);
        fprintf(stdout, "Ndepth = %zu\n", Ndepth);
        fprintf(stdout, "\n");
    }

    if ( (Nlon > 1e4) or (Nlat > 1e4) or (Ntime > 1e2) or (Ndepth > 1e2) ) {
        if ((retval = nc_close(ncid))) { NC_ERR(retval, __LINE__, __FILE__); }
        fprintf(stdout, "Data dimensions too large to continue.\n");
        return;
    }

    //
    //// Allocate memory for the fields
    //

    longitude[0] = new double[Nlon];
    latitude[0]  = new double[Nlat];
    time[0]      = new double[Ntime];
    depth[0]     = new double[Ndepth];

    u_lon[0] = new double[Ntime * Ndepth * Nlat * Nlon];
    u_lat[0] = new double[Ntime * Ndepth * Nlat * Nlon];

    //
    //// Get fields from IC file
    //

    // Define coordinate variables
    int lon_varid, lat_varid, time_varid, depth_varid;
    if ((retval = nc_inq_varid(ncid, "longitude", &lon_varid   ))) { NC_ERR(retval, __LINE__, __FILE__); }
    if ((retval = nc_inq_varid(ncid, "latitude",  &lat_varid   ))) { NC_ERR(retval, __LINE__, __FILE__); }
    if ((retval = nc_inq_varid(ncid, "time",      &time_varid  ))) { NC_ERR(retval, __LINE__, __FILE__); }
    if ((retval = nc_inq_varid(ncid, "depth",     &depth_varid ))) { NC_ERR(retval, __LINE__, __FILE__); }

    // Declare variables
    int ulon_varid, ulat_varid;
    if ((retval = nc_inq_varid(ncid, "uo", &ulon_varid))) { NC_ERR(retval, __LINE__, __FILE__); }
    if ((retval = nc_inq_varid(ncid, "vo", &ulat_varid))) { NC_ERR(retval, __LINE__, __FILE__); }

    // Get the coordinate variables
    size_t start_lon[1], count_lon[1],
           start_lat[1], count_lat[1];

    start_lon[0] = 0;
    count_lon[0] = Nlon;

    start_lat[0] = 0;
    count_lat[0] = Nlat;

    if ((retval = nc_get_vara_double(ncid, lon_varid, start_lon, count_lon, longitude[0]))) { NC_ERR(retval, __LINE__, __FILE__); }
    if ((retval = nc_get_vara_double(ncid, lat_varid, start_lat, count_lat, latitude[0] ))) { NC_ERR(retval, __LINE__, __FILE__); }

    // Get u, v, and h
    size_t start[4], count[4];
    start[0] = 0;
    start[1] = 0;
    start[2] = 0;
    start[3] = 0;
    count[0] = Ntime;
    count[1] = Ndepth;
    count[2] = Nlat;
    count[3] = Nlon;

    if ((retval = nc_get_vara_double(ncid, ulon_varid, start, count, u_lon[0]))) { NC_ERR(retval, __LINE__, __FILE__); }
    if ((retval = nc_get_vara_double(ncid, ulat_varid, start, count, u_lat[0]))) { NC_ERR(retval, __LINE__, __FILE__); }

    // Close the file
    if ((retval = nc_close(ncid))) { NC_ERR(retval, __LINE__, __FILE__); }

}
