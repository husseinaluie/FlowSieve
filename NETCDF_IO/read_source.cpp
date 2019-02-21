/*
 *
 * Read source data (from input.nc) for coarse graining.
 *
 * Current assumptions:
 *    Variables to read: uo (as u_lon), vo (as u_lat)
 *    Dimensions: time, depth, longitude, latitude (in that order)
 *
 */

#include "../netcdf_io.hpp"

#ifndef DEBUG
    #define DEBUG false
#endif


// Write to netcdf file
void read_source(
        int & Nlon,          int & Nlat,
        int & Ntime,         int & Ndepth,
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

    size_t Ntime_st, Ndepth_st, Nlon_st, Nlat_st;
    if ((retval = nc_inq_dim(ncid, time_dimid , NULL, &Ntime_st  ))) { NC_ERR(retval, __LINE__, __FILE__); }
    if ((retval = nc_inq_dim(ncid, depth_dimid, NULL, &Ndepth_st ))) { NC_ERR(retval, __LINE__, __FILE__); }
    if ((retval = nc_inq_dim(ncid, lat_dimid  , NULL, &Nlon_st   ))) { NC_ERR(retval, __LINE__, __FILE__); }
    if ((retval = nc_inq_dim(ncid, lon_dimid  , NULL, &Nlat_st   ))) { NC_ERR(retval, __LINE__, __FILE__); }


    // Cast the sizes to integers (to resolve some compile errors)
    //   Unless we're dealing with truly massive grids, this
    //   shouldn't be an issue.
    Ntime  = static_cast<int>(Ntime_st);
    Ndepth = static_cast<int>(Ndepth_st);
    Nlon   = static_cast<int>(Nlon_st);
    Nlat   = static_cast<int>(Nlat_st);
    if (debug) {
        fprintf(stdout, "\n");
        fprintf(stdout, "Nlon   = %d\n", Nlon);
        fprintf(stdout, "Nlat   = %d\n", Nlat);
        fprintf(stdout, "Ntime  = %d\n", Ntime);
        fprintf(stdout, "Ndepth = %d\n", Ndepth);
        fprintf(stdout, "\n");
    }

    // For the moment, as a precaution stop if we hit something too large.
    if ( (Nlon > 1e4) or (Nlat > 1e4) or (Ntime > 1e2) or (Ndepth > 1e2) ) {
        if ((retval = nc_close(ncid))) { NC_ERR(retval, __LINE__, __FILE__); }
        fprintf(stdout, "Data dimensions too large to continue. (Line %d of %s)\n", __LINE__, __FILE__);
        return;
    }

    //
    //// Allocate memory for the fields
    //

    time[0]      = new double[Ntime];
    depth[0]     = new double[Ndepth];
    longitude[0] = new double[Nlon];
    latitude[0]  = new double[Nlat];

    u_lon[0] = new double[Ntime * Ndepth * Nlat * Nlon];
    u_lat[0] = new double[Ntime * Ndepth * Nlat * Nlon];

    //
    //// Get fields from IC file
    //

    // Define coordinate variables
    int lon_varid, lat_varid, time_varid, depth_varid;
    if ((retval = nc_inq_varid(ncid, "time",      &time_varid  ))) { NC_ERR(retval, __LINE__, __FILE__); }
    if ((retval = nc_inq_varid(ncid, "depth",     &depth_varid ))) { NC_ERR(retval, __LINE__, __FILE__); }
    if ((retval = nc_inq_varid(ncid, "longitude", &lon_varid   ))) { NC_ERR(retval, __LINE__, __FILE__); }
    if ((retval = nc_inq_varid(ncid, "latitude",  &lat_varid   ))) { NC_ERR(retval, __LINE__, __FILE__); }

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

    // Get u_lon (uo) and u_lat (vo)
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
