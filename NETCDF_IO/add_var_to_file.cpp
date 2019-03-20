
#include <vector>
#include "../netcdf_io.hpp"

#ifndef DEBUG
    #define DEBUG 0
#endif

void add_var_to_file(
        const char * var_name,   /**< [in] depth vector (1D) */
        const char * filename,   /**< [in] longitude vector (1D) */
        const char ** dim_list,                 /**< [in] list of dimensions (in order!) */
        const int num_dims                      /**< [in] number of dimensions */
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
    snprintf(varname, 50, var_name);
    retval = nc_def_var(ncid, varname, NC_DOUBLE, num_dims, dim_ids, &var_id);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

    // Close the file
    retval = nc_close(ncid);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

    #if DEBUG >= 2
    fprintf(stdout, "   added %s to %s.\n", varname, buffer);
    #endif

}