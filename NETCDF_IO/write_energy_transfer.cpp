#include "../netcdf_io.hpp"

#ifndef DEBUG
    #define DEBUG 0
#endif

void write_energy_transfer(
        const double * energy_transfer, /**< [in] vort_lat to be written to the file*/
        const int Iscale,               /**< [in] Index positioning this output in the filter dimension */
        const int Ntime,                /**< [in] Length of the time dimension */
        const int Ndepth,               /**< [in] Length of the depth dimension */
        const int Nlat,                 /**< [in] Length of the latitude dimension */  
        const int Nlon                  /**< [in] Length of the longitude dimension */
        ) {

    // Open the NETCDF file
    int FLAG = NC_WRITE;
    int ncid=0, retval;
    char buffer [50];
    snprintf(buffer, 50, "filter_output.nc");
    if (( retval = nc_open(buffer, FLAG, &ncid) ))
        NC_ERR(retval, __LINE__, __FILE__);

    // Get the variable IDs for the outputs
    int ener_transfer_varid;
    if ((retval = nc_inq_varid(ncid, "energy_transfer",   &ener_transfer_varid ))) { NC_ERR(retval, __LINE__, __FILE__); }

    // Write the current scale to the output
    size_t start[5], count[5];

    start[0] = Iscale;
    start[1] = 0;
    start[2] = 0;
    start[3] = 0;
    start[4] = 0;

    count[0] = 1;
    count[1] = Ntime;
    count[2] = Ndepth;
    count[3] = Nlat;
    count[4] = Nlon;

    if ((retval = nc_put_vara_double(ncid, ener_transfer_varid, start, count, energy_transfer)))
        NC_ERR(retval, __LINE__, __FILE__);

    // Close the file
    if ((retval = nc_close(ncid))) { NC_ERR(retval, __LINE__, __FILE__); }

    #if DEBUG >= 2
    fprintf(stdout, "-Energy transfer written-\n");
    #endif
}
