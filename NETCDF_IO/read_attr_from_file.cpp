
#include "../netcdf_io.hpp"
#include "../constants.hpp"
#include <math.h>

// Write to netcdf file
void read_attr_from_file(
        double &attr,           /**< [in] var into which to store the variable */
        const char * attr_name, /**< [in] Name of the attribute */
        const char * filename,  /**< [in] Name of the file */
        const char * var_name,  /**< [in] Name of associated variable */
        const MPI_Comm comm     /**< [in] MPI Communicator */
        ) {

    int wRank, wSize;
    MPI_Comm_rank( comm, &wRank );
    MPI_Comm_size( comm, &wSize );

    #if DEBUG >= 1
    if (wRank == 0) {
        fprintf(stdout, "Attempting to read %s from %s\n", attr_name, filename);
    }
    #endif

    // Open the NETCDF file
    //int FLAG = NC_NETCDF4 | NC_NOWRITE | NC_MPIIO;
    int FLAG = NC_NETCDF4 | NC_MPIIO;
    int ncid=0, retval;
    char buffer [50];
    snprintf(buffer, 50, filename);
    retval = nc_open_par(buffer, FLAG, comm, MPI_INFO_NULL, &ncid);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

    // Get information about the variable
    int var_id = NC_GLOBAL, num_dims;
    int dim_ids[NC_MAX_VAR_DIMS];
    if (var_name == NULL) {
        // If no var_name given, then assume a global attribute
        var_id = NC_GLOBAL;
    } else {
        // Otherwise, get the appropriate variable id
        retval = nc_inq_var(ncid, var_id, NULL, NULL, &num_dims, dim_ids, NULL );
        if (retval) { NC_ERR(retval, __LINE__, __FILE__); }
    }

    char attrname [50];
    snprintf(attrname, 50, attr_name);

    nc_get_att(ncid, var_id, attrname, &attr);

    MPI_Barrier(comm);
    retval = nc_close(ncid);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }
}
