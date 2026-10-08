#include <fenv.h>
#include <vector>
#include <string>
#include <mpi.h>
#include <math.h>
#include <cassert>
#include <fenv.h>
#include "../netcdf_io.hpp"
#include "../constants.hpp"

void write_mask_to_output(
        const std::vector<short int> & mask,
        const std::string & field_name,
        const size_t * start,
        const size_t * count,
        const std::string & filename,
        MPI_Comm comm
        ) {

    // During writing, ignore floating point exceptions
    std::fenv_t fe_env;
    feholdexcept( &fe_env );

    int wRank, wSize;
    MPI_Comm_rank( comm, &wRank );
    MPI_Comm_size( comm, &wSize );

    #if DEBUG >= 2
    if (wRank == 0) { fprintf(stdout, "  Preparing to write %s to %s (%'zu points).\n", field_name.c_str(), filename.c_str(), mask.size()); }
    fflush(stdout);
    #endif

    // Open the NETCDF file
    int FLAG = NC_NETCDF4 | NC_WRITE | NC_MPIIO;
    int ncid=0, retval;

    // Some netcdf functions [in some netcdf versions] cause floating-point errors
    //  so, we need to disable floating point exceptions when we try to open files.
    fedisableexcept( FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW );
    retval = nc_open_par( filename.c_str(), FLAG, comm, MPI_INFO_NULL, &ncid);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

    // Now we can restore fp-exception handling, after clearing out any that were raised
    feclearexcept( FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW ); // erase whatever exceptions were raised
    feenableexcept( FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW ); // re-enable exceptions

    // Get the variable ID for the field
    int field_varid;
    retval = nc_inq_varid(ncid, field_name.c_str(), &field_varid );
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

    // Get information about the variable
    int num_dims;
    int dim_ids[NC_MAX_VAR_DIMS];
    retval = nc_inq_var(ncid, field_varid, NULL, NULL, &num_dims, dim_ids, NULL );
    if (retval != NC_NOERR ) { NC_ERR(retval, __LINE__, __FILE__); }
    assert( num_dims > 0 );
    
    /*
    size_t num_pts = 1;
    for (int Idim = 0; Idim < num_dims; Idim++) {
        num_pts *= count[Idim];
    }
    */

    retval = nc_put_vara_short( ncid, field_varid, start, count, &(mask[0]) );
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

    // Close the file
    MPI_Barrier(comm);
    retval = nc_close(ncid);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

    #if DEBUG >= 1
    if (wRank == 0) { 
        fprintf(stdout, "    wrote %s to %s \n", field_name.c_str(), filename.c_str());
        fflush(stdout);
    }
    #endif

    fesetenv( &fe_env );
}
