#include <vector>
#include <string>
#include <mpi.h>
#include <math.h>
#include "../netcdf_io.hpp"
#include "../constants.hpp"

void write_integral_to_post(
        const std::vector<
            std::vector<double> > & field,  /**< [in] data to be written to the file*/
        std::string field_name,             /**< [in] name of the variable in the netcdf file */
        std::string field_suffix,           /**< [in] name of the variable in the netcdf file */
        size_t * start,                     /**< [in] starting indices for the write */
        size_t * count,                     /**< [in] size of the write in each dimension */
        const char * filename,              /**< [in] name of the netcdf file */
        const MPI_Comm comm                 /**< [in] MPI Communicator */
        ) {

    int wRank, wSize;
    MPI_Comm_rank( comm, &wRank );
    MPI_Comm_size( comm, &wSize );

    // Open the NETCDF file
    int FLAG = NC_NETCDF4 | NC_WRITE | NC_MPIIO;
    int ncid=0, retval;
    char buffer [50];
    snprintf(buffer, 50, filename);

    MPI_Barrier(comm);
    retval = nc_open_par(buffer, FLAG, comm, MPI_INFO_NULL, &ncid);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

    // Get the variable ID for the field
    int field_varid;
    retval = nc_inq_varid(ncid, (field_name + field_suffix).c_str(), &field_varid );
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

    // Write the integral for each region
    for (size_t Iregion = 0; Iregion < field.size(); ++Iregion) {
        start[2] = Iregion;

        retval = nc_put_vara_double(ncid, field_varid, start, count, 
                &(field.at(Iregion)[0]));
        if (retval) { NC_ERR(retval, __LINE__, __FILE__); }
    }

    // Close the file
    MPI_Barrier(comm);
    retval = nc_close(ncid);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

    #if DEBUG >= 1
    if (wRank == 0) { fprintf(stdout, "  - wrote %s to %s -\n", 
            (field_name + field_suffix).c_str(), filename); }
    #endif
}
