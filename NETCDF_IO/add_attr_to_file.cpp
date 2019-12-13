#include <vector>
#include <string>
#include <mpi.h>
#include "../netcdf_io.hpp"
#include "../constants.hpp"

void add_attr_to_file(
        const char * varname,
        const double value,
        const char * filename,
        const MPI_Comm comm                     /**< [in] MPI Communicator */
        ) {

    int wRank=-1, wSize=-1;
    MPI_Comm_rank( MPI_COMM_WORLD, &wRank );
    MPI_Comm_size( MPI_COMM_WORLD, &wSize );

    if (wRank == 0) {
        // Open the NETCDF file
        int FLAG = NC_WRITE;
        int ncid=0, retval;
        char buffer [50];
        snprintf(buffer, 50, filename);
        retval = nc_open(buffer, FLAG, &ncid);
        if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

        retval = nc_put_att_double(ncid, NC_GLOBAL, varname, NC_DOUBLE, 1, &value);
        if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

        // Close the file
        retval = nc_close(ncid);
        if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

        fprintf(stdout, "  - added %s to %s -\n", varname, buffer);
    }

}
