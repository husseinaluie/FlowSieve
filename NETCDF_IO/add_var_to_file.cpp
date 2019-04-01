
#include <vector>
#include "../netcdf_io.hpp"
#include "../constants.hpp"

void add_var_to_file(
        const char * var_name,   /**< [in] variable name */
        const char ** dim_list,  /**< [in] list of dimensions (in order!) */
        const int num_dims,      /**< [in] number of dimensions */
        const char * filename    /**< [in] file name */
        ) {

    // Open the NETCDF file
    int FLAG = NC_WRITE;
    int ncid=0, retval;
    char buffer [50];
    snprintf(buffer, 50, filename);
    retval = nc_open(buffer, FLAG, &ncid);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

    // Extract dimension ids sizes
    int dim_ids[num_dims];
    for (int dim_ind = 0; dim_ind < num_dims; dim_ind++) {
        retval = nc_inq_dimid(ncid, dim_list[dim_ind], &dim_ids[dim_ind] );
        if (retval) { NC_ERR(retval, __LINE__, __FILE__); }
    }

    // Declare the variable
    int var_id;
    char varname [50];
    double fill_value = -32767;
    snprintf(varname, 50, var_name);
    retval = nc_def_var(ncid, varname, NC_FLOAT, num_dims, dim_ids, &var_id);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }
    retval = nc_put_att_double(ncid, var_id, "_FillValue", NC_FLOAT, 1, &fill_value);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

    // Close the file
    retval = nc_close(ncid);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

    #if DEBUG >= 2
    fprintf(stdout, "  - added %s to %s -\n", varname, buffer);
    #endif

}
