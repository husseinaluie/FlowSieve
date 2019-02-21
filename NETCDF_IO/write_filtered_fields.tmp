#include "../netcdf_io.hpp"

// Write to netcdf file
void write_filtered_fields(int iter, double curr_time, 
        double * my_x, double * my_y,
        double * my_u, double * my_v, double * my_h, double * my_vort,
        int bRank, int my_rank_x, int my_rank_y, 
        int my_Nx, int my_Ny, int full_Nx, int full_Ny) {

    // Open the NETCDF file
    int FLAG = NC_NETCDF4 | NC_CLOBBER;
    int ncid=0, retval;
    char buffer [50];
    snprintf(buffer, 50, "output_%04d.nc", iter);
    if (( retval = nc_create(buffer, FLAG, &ncid) ))
        NC_ERR(retval, __LINE__, __FILE__);

    // Define the dimensions
    int xdimid, ydimid, tdimid;
    if ((retval = nc_def_dim(ncid, "x",    full_Nx, &xdimid)))
        NC_ERR(retval, __LINE__, __FILE__);
    if ((retval = nc_def_dim(ncid, "y",    full_Ny, &ydimid)))
        NC_ERR(retval, __LINE__, __FILE__);
    if ((retval = nc_def_dim(ncid, "time", 1,       &tdimid)))
        NC_ERR(retval, __LINE__, __FILE__);

    // Define coordinate variables
    int xvarid, yvarid, tvarid;
    if ((retval = nc_def_var(ncid, "x",    NC_DOUBLE, 1, &xdimid, &xvarid)))
        NC_ERR(retval, __LINE__, __FILE__);
    if ((retval = nc_def_var(ncid, "y",    NC_DOUBLE, 1, &ydimid, &yvarid)))
        NC_ERR(retval, __LINE__, __FILE__);
    if ((retval = nc_def_var(ncid, "time", NC_DOUBLE, 1, &tdimid, &tvarid)))
        NC_ERR(retval, __LINE__, __FILE__);

    // Transpose
    const int ndims = 2;
    int dimids[ndims];
    dimids[0] = ydimid;
    dimids[1] = xdimid;

    // Declare variables
    int uvarid, vvarid, hvarid, vortvarid;
    if ((retval = nc_def_var(ncid, "u", NC_DOUBLE, ndims, dimids, &uvarid)))
        NC_ERR(retval, __LINE__, __FILE__);
    if ((retval = nc_def_var(ncid, "v", NC_DOUBLE, ndims, dimids, &vvarid)))
        NC_ERR(retval, __LINE__, __FILE__);
    if ((retval = nc_def_var(ncid, "h", NC_DOUBLE, ndims, dimids, &hvarid)))
        NC_ERR(retval, __LINE__, __FILE__);
    if ((retval = nc_def_var(ncid, "vort", NC_DOUBLE, ndims, dimids, &vortvarid)))
        NC_ERR(retval, __LINE__, __FILE__);

    // Put the coordinate variables
    size_t start_x[1], count_x[1];
    start_x[0] = my_Nx * my_rank_x;
    count_x[0] = my_Nx;
    if (my_rank_y == 0) { 
        if ((retval = nc_put_vara_double(ncid, xvarid, start_x, count_x, my_x)))
            NC_ERR(retval, __LINE__, __FILE__);
    }
    size_t start_y[1], count_y[1];
    start_y[0] = my_Ny * my_rank_y;
    count_y[0] = my_Ny;
    if (my_rank_x == 0) { 
        if ((retval = nc_put_vara_double(ncid, yvarid, start_y, count_y, my_y)))
            NC_ERR(retval, __LINE__, __FILE__);
    }
    size_t start_t[1], count_t[1];
    start_t[0] = 0;
    count_t[0] = 1;
    if (bRank == 0) {
        if ((retval = nc_put_vara_double(ncid, tvarid, start_t, count_t, &curr_time)))
            NC_ERR(retval, __LINE__, __FILE__);
    }

    size_t start[2], count[2];
    start[0] = my_Ny * my_rank_y;
    start[1] = my_Nx * my_rank_x;
    count[0] = my_Ny;
    count[1] = my_Nx;

    if ((retval = nc_put_vara_double(ncid, uvarid, start, count, my_u)))
        NC_ERR(retval, __LINE__, __FILE__);
    if ((retval = nc_put_vara_double(ncid, vvarid, start, count, my_v)))
        NC_ERR(retval, __LINE__, __FILE__);
    if ((retval = nc_put_vara_double(ncid, hvarid, start, count, my_h)))
        NC_ERR(retval, __LINE__, __FILE__);
    if ((retval = nc_put_vara_double(ncid, vortvarid, start, count, my_vort)))
        NC_ERR(retval, __LINE__, __FILE__);

    // Close the file
    if ((retval = nc_close(ncid))) { NC_ERR(retval, __LINE__, __FILE__); }

}
