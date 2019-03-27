#include <vector>
#include "../netcdf_io.hpp"
#include "../constants.hpp"

void write_field_to_output(
        const std::vector<double> & field,  /**< [in] transfer to be written to the file*/
        const char * field_name,            /**< [in] name of the variable in the netcdf file */
        const size_t * start,               /**< [in] starting indices for the write */
        const size_t * count,               /**< [in] size of the write in each dimension */
        const char * filename               /**< [in] filename */
        ) {

    // Open the NETCDF file
    int FLAG = NC_WRITE;
    int ncid=0, retval;
    char buffer [50];
    snprintf(buffer, 50, filename);
    if (( retval = nc_open(buffer, FLAG, &ncid) ))
        NC_ERR(retval, __LINE__, __FILE__);

    // Get the variable ID for the field
    int field_varid;
    if ((retval = nc_inq_varid(ncid, field_name, &field_varid ))) { NC_ERR(retval, __LINE__, __FILE__); }

    // Write the current scale to the output
    if ((retval = nc_put_vara_double(ncid, field_varid, start, count, &field[0])))
        NC_ERR(retval, __LINE__, __FILE__);

    // Close the file
    if ((retval = nc_close(ncid))) { NC_ERR(retval, __LINE__, __FILE__); }

    #if DEBUG >= 2
    fprintf(stdout, "  - wrote %s to file -\n", field_name);
    #endif
}
