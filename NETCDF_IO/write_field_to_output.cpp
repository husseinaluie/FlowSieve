#include <vector>
#include <string>
#include <mpi.h>
#include "../netcdf_io.hpp"
#include "../constants.hpp"

void write_field_to_output(
        const std::vector<double> & field,  /**< [in] data to be written to the file*/
        const char * field_name,            /**< [in] name of the variable in the netcdf file */
        const size_t * start,               /**< [in] starting indices for the write */
        const size_t * count,               /**< [in] size of the write in each dimension */
        const char * filename,              /**< [in] name of the netcdf file */
        MPI_Comm comm                       /**< [in] MPI Communicator */
        ) {

    int wRank, wSize;
    MPI_Comm_rank( comm, &wRank );
    MPI_Comm_size( comm, &wSize );

    // Open the NETCDF file
    int FLAG = NC_NETCDF4 | NC_WRITE | NC_MPIIO;
    int ncid=0, retval;
    char buffer [50];
    snprintf(buffer, 50, filename);
    //if (( retval = nc_open(buffer, FLAG, &ncid) ))
    if (( retval = nc_open_par(buffer, FLAG, comm, MPI_INFO_NULL, &ncid) ))
        NC_ERR(retval, __LINE__, __FILE__);

    // Get the variable ID for the field
    int field_varid;
    if (( retval = nc_inq_varid(ncid, field_name, &field_varid ) )) 
        NC_ERR(retval, __LINE__, __FILE__);

    // Write the current scale to the output
    if (( retval = nc_put_vara_double(ncid, field_varid, start, count, &field[0]) ))
        NC_ERR(retval, __LINE__, __FILE__);

    // Close the file
    MPI_Barrier(comm);
    if (( retval = nc_close(ncid) )) 
        NC_ERR(retval, __LINE__, __FILE__);

    #if DEBUG >= 2
    if (wRank == 0) { fprintf(stdout, "  - wrote %s to %s -\n", field_name, filename); }
    #endif
}
