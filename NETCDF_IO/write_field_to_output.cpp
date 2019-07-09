#include <vector>
#include <string>
#include <mpi.h>
#include <math.h>
#include "../netcdf_io.hpp"
#include "../constants.hpp"

void write_field_to_output(
        const std::vector<double> & field,  /**< [in] data to be written to the file*/
        const char * field_name,            /**< [in] name of the variable in the netcdf file */
        const size_t * start,               /**< [in] starting indices for the write */
        const size_t * count,               /**< [in] size of the write in each dimension */
        const char * filename,              /**< [in] name of the netcdf file */
        const std::vector<double> * mask,   /**< [in] (pointer to) mask */
        MPI_Comm comm                       /**< [in] MPI Communicator */
        ) {

    int wRank, wSize;
    MPI_Comm_rank( comm, &wRank );
    MPI_Comm_size( comm, &wSize );

    #if DEBUG >= 2
    fprintf(stdout, "  Rank %d: starts = %zu %zu %zu %zu\n", wRank, 
            start[0], start[1], start[2], start[3]);
    fprintf(stdout, "  Rank %d: counts = %zu %zu %zu %zu\n", wRank,
            count[0], count[1], count[2], count[3]);
    #endif

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
    retval = nc_inq_varid(ncid, field_name, &field_varid );
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

    std::vector<signed short> reduced_field; 
    double add_offset, scale_factor;
    if (constants::CAST_TO_INT) {
        // If we want to reduce output size, pack into short ints
        //   floats are 32bit, short ints are 16bit, so we can cut
        //   file size in half. Of course, this is at the cost
        //   of precision.

        #if DEBUG >= 2
        if (wRank == 0) { fprintf(stdout, "  Initializing reduced field storage.\n"); }
        fflush(stdout);
        #endif
        reduced_field.resize(field.size());

        #if DEBUG >= 2
        if (wRank == 0) { fprintf(stdout, "  Preparing to package the field.\n"); }
        fflush(stdout);
        #endif
        package_field(reduced_field, scale_factor, add_offset, field, mask);

        // We need to record the scale and translation used to encode in signed shorts
        retval = nc_put_att_double(
                ncid, field_varid, "scale_factor", NC_FLOAT, 1, &scale_factor);
        if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

        retval = nc_put_att_double(
                ncid, field_varid, "add_offset", NC_FLOAT, 1, &add_offset);
        if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

        retval = nc_put_vara_short(ncid, field_varid, start, count, &reduced_field[0]);
        if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

    } else {

        // Otherwise, just write the 32-bit float field to the file
        retval = nc_put_vara_double(ncid, field_varid, start, count, &field[0]);
        if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

    }

    // Close the file
    MPI_Barrier(comm);
    retval = nc_close(ncid);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

    #if DEBUG >= 1
    if (wRank == 0) { fprintf(stdout, "  - wrote %s to %s -\n", field_name, filename); }
    #endif
}
