
#include "../netcdf_io.hpp"
#include "../constants.hpp"
#include <math.h>

// Write to netcdf file
void read_var_from_file(
        std::vector<double> &var,   /**< [in] Vector into which to store the variable */
        const char * var_name,      /**< [in] Name of the variable */
        const char * filename,      /**< [in] Name of the file */
        std::vector<double> *mask,  /**< [in] Pointer to mask array to be created */
        std::vector<int> *myCounts, /**< [in] Vector into which to store sizes (if not NULL) */
        std::vector<int> *myStarts, /**< [in] Vector into which to store sizes (if not NULL) */
        const MPI_Comm comm         /**< [in] MPI Communicator */
        ) {

    int wRank, wSize;
    MPI_Comm_rank( comm, &wRank );
    MPI_Comm_size( comm, &wSize );

    #if DEBUG >= 1
    if (wRank == 0) {
        fprintf(stdout, "Attempting to read %s from %s\n", var_name, filename);
    }
    #endif

    // Open the NETCDF file
    //int FLAG = NC_NETCDF4 | NC_NOWRITE | NC_MPIIO;
    int FLAG = NC_NETCDF4 | NC_MPIIO;
    int ncid=0, retval;
    char buffer [50];
    snprintf(buffer, 50, filename);
    //if (( retval = nc_open(buffer, FLAG, &ncid) )) 
    if (wRank == 0) { fprintf(stdout, "About to open file.\n");  }
    if (( retval = nc_open_par(buffer, FLAG, comm, MPI_INFO_NULL, &ncid) ))
        NC_ERR(retval, __LINE__, __FILE__); 
    if (wRank == 0) { fprintf(stdout, "File opened.\n");  }

    char varname [50];
    snprintf(varname, 50, var_name);

    int var_id, num_dims;
    int dim_ids[NC_MAX_VAR_DIMS];

    // Get the ID for the variable
    if ((retval = nc_inq_varid(ncid, varname, &var_id ))) 
        NC_ERR(retval, __LINE__, __FILE__);

    // Get information about the variable
    if ((retval = nc_inq_var(ncid, var_id, NULL, NULL, &num_dims, dim_ids, NULL ))) 
        NC_ERR(retval, __LINE__, __FILE__);
    #if DEBUG >= 2
    if (wRank == 0) {
        if (num_dims == 1) {
            fprintf(stdout, "  has %d dimension of size ", num_dims);
        } else {
            fprintf(stdout, "  has %d dimensions of size ", num_dims);
        }
    }
    #endif

    // Get the size of each dimension
    size_t start[num_dims], count[num_dims];
    size_t num_pts = 1;
    int my_count;
    if (myCounts != NULL) {
        myCounts->resize(num_dims);
        myStarts->resize(num_dims);
    }
    for (int II = 0; II < num_dims; II++) {
        start[II] = 0;
        if ((retval = nc_inq_dim(ncid, dim_ids[II] , NULL, &count[II]  ))) 
            NC_ERR(retval, __LINE__, __FILE__);
        #if DEBUG >= 2
        fprintf(stdout, "%zu ", count[II]);
        #endif

        // If we're split on multiple MPI procs and have > 2 dimensions, 
        //   then divide all but the last two we don't split the last 
        //   two because those are assumed to be lat/lon
        if ( (num_dims > 2) and (wSize > 1) and (II == 0) ) {
            #if DEBUG >= 1
            if (wRank == 0) {
                fprintf(stdout, "More than two dimensions and more than 1 MPI proc, "
                        "so will split the first %d dimensions.\n", num_dims - 2);
            }
            #endif
            // For now, just split in time (assumed to be the first dimension)
            my_count = ((int)count[0]) / wSize;
            start[II] = (size_t) (wRank * my_count);
            if (wRank == wSize - 1) { my_count = (int)(count[II] - start[II]); }
            if (myCounts != NULL) { myCounts->at(II) = my_count;  }
            if (myStarts != NULL) { myStarts->at(II) = (int) start[II];  }
            #if DEBUG >= 2
            fprintf(stdout, "   proc %d of %d taking %d of %zu points in dimension %d.\n",
                    wRank+1, wSize, my_count, count[0], II);
            #endif
            count[II] = (size_t) my_count;
        }
        num_pts *= count[II];
    }
    #if DEBUG >= 2
    fprintf(stdout, "\n");
    #endif

    // Now resize the vector to the appropriate size
    var.resize(num_pts);

    // Now read in the data
    if ((retval = nc_get_vara_double(ncid, var_id, start, count, &var[0]))) 
        NC_ERR(retval, __LINE__, __FILE__);

    // Apply scale factor if appropriate
    double scale = 1.;
    if ((retval = nc_get_att_double(ncid, var_id, "scale_factor", &scale))) 
        NC_ERR(retval, __LINE__, __FILE__);
    if (scale != 1.) {
        #if DEBUG >= 2
        if (wRank == 0) { fprintf(stdout, "  scale factor = %g\n", scale); }
        #endif
        for (size_t II = 0; II < num_pts; II++) {
            var.at(II) = var.at(II) * scale;
        }
    }

    // Apply offset if appropriate
    double offset = 0.;
    if ((retval = nc_get_att_double(ncid, var_id, "add_offset", &offset))) 
        NC_ERR(retval, __LINE__, __FILE__);
    if (offset != 0.) {
        #if DEBUG >= 2
        if (wRank == 0) { fprintf(stdout, "  additive offset = %g\n", offset); }
        #endif
        for (size_t II = 0; II < num_pts; II++) {
            var.at(II) = var.at(II) + offset;
        }
    }

    // Determine masking, if desired
    double fill_val = 1e100;
    int num_land = 0;
    int num_water = 0;
    size_t Nlat, Nlon;
    if ( (mask != NULL) and (num_dims == 4) ) {
        // Assuming CF dimension order: time - depth - lat - lon
        Nlat = count[2];
        Nlon = count[3];
        mask->resize(Nlat*Nlon);
        if ((retval = nc_get_att_double(ncid, var_id, "_FillValue", &fill_val))) 
            NC_ERR(retval, __LINE__, __FILE__);
        #if DEBUG >= 2
        if (wRank == 0) { fprintf(stdout, "  fill value = %g\n", fill_val); }
        #endif
        for (size_t II = 0; II < Nlat * Nlon; II++) {
            if (fabs(var.at(II)) > 0.95 * fabs(fill_val*scale)) {
                mask->at(II) = 0;
                num_land++;
            } else {
                mask->at(II) = 1;
                num_water++;
            }
        }
        #if DEBUG >= 1
        if (wRank == 0) {
            fprintf(stdout, "  Land cover = %.4g%%\n", 
                    100 * ((double)num_land) / (num_land + num_water));
        }
        #endif
    }
}
