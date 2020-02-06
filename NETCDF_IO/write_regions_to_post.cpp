#include <vector>
#include <string>
#include <mpi.h>
#include <math.h>
#include "../netcdf_io.hpp"
#include "../constants.hpp"
#include "../postprocess.hpp"

void write_regions_to_post(
        const char * filename,
        const MPI_Comm comm
        ) {

    int wRank, wSize;
    MPI_Comm_rank( comm, &wRank );
    MPI_Comm_size( comm, &wSize );

    // Unfortunately, writing string variables doesn't seem to work
    //    under parallel IO. It's a small dimension, so we'll just
    //    do it serially.
    if (wRank == 0) {

        // Open the NETCDF file
        int FLAG = NC_NETCDF4 | NC_WRITE;
        int ncid=0, retval;
        char buffer [50];
        snprintf(buffer, 50, filename);

        retval = nc_open(buffer, FLAG, &ncid);
        if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

        // Get the variable ID for the field
        int region_varid;
        retval = nc_inq_varid(ncid, "region", &region_varid );
        if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

        const size_t num_regions = RegionTest::region_names.size();

        // This is hackey and awful, but it seems to be necessary
        //    to get from the "vector of strings" to "const char **"
        //    that netcdf seems to want.
        // We're also going to have to go through and write them
        //    one at a time. Woo.
        char curr_name[50];
        const char *ptr_to_name = &curr_name[0];

        size_t start[1], count[1];
        count[0] = 1;
        for (size_t Iregion = 0; Iregion < num_regions; ++Iregion) {
            start[0] = Iregion;

            sprintf(curr_name, "%20s", RegionTest::region_names.at(Iregion).c_str());

            retval = nc_put_vara_string(ncid, region_varid, start, count, &ptr_to_name);
            if (retval) { NC_ERR(retval, __LINE__, __FILE__); }
        }

        // Close the file
        retval = nc_close(ncid);
        if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

        #if DEBUG >= 1
        if (wRank == 0) { fprintf(stdout, "  - wrote region names to %s -\n", filename); }
        #endif
    }
}
