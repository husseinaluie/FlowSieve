#include <vector>
#include <string>
#include <mpi.h>
#include <math.h>
#include "../netcdf_io.hpp"
#include "../constants.hpp"

void write_time_average_to_post(
        const std::vector< double > & field,
        std::string field_name,
        std::string field_suffix,
        size_t * start,
        size_t * count,
        const char * filename,
        const MPI_Comm comm
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

    // Write the time average
    retval = nc_put_vara_double(ncid, field_varid, start, count, &(field[0]));
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

    // Close the file
    MPI_Barrier(comm);
    retval = nc_close(ncid);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

    #if DEBUG >= 1
    if (wRank == 0) { fprintf(stdout, "  - wrote %s to %s -\n", 
            (field_name + field_suffix).c_str(), filename); }
    #endif
}

