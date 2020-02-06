
#include "../netcdf_io.hpp"
#include "../constants.hpp"
#include <math.h>

// Write to netcdf file
void read_var_from_file(
        std::vector<double> &var,
        const char * var_name,
        const char * filename,
        std::vector<double> *mask,
        std::vector<int> *myCounts,
        std::vector<int> *myStarts,
        const bool do_splits,
        const MPI_Comm comm
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
    retval = nc_open_par(buffer, FLAG, comm, MPI_INFO_NULL, &ncid);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

    char varname [50];
    snprintf(varname, 50, var_name);

    int var_id, num_dims;
    int dim_ids[NC_MAX_VAR_DIMS];

    // Get the ID for the variable
    retval = nc_inq_varid(ncid, varname, &var_id );
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

    // Get information about the variable
    retval = nc_inq_var(ncid, var_id, NULL, NULL, &num_dims, dim_ids, NULL );
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }
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
    int my_count, overflow;
    if (myCounts != NULL) {
        myCounts->resize(num_dims);
        myStarts->resize(num_dims);
    }
    for (int II = 0; II < num_dims; II++) {
        start[II] = 0;
        retval = nc_inq_dim(ncid, dim_ids[II] , NULL, &count[II]);
        if (retval) { NC_ERR(retval, __LINE__, __FILE__); }
        #if DEBUG >= 2
        if (wRank == 0) { fprintf(stdout, "%zu ", count[II]); }
        #endif

        if (do_splits) {
            // If we're split on multiple MPI procs and have > 2 dimensions, 
            //   then divide all but the last two 
            //
            //   we don't split the last two because those 
            //   are assumed to be lat/lon
            if ( (num_dims > 2) and (wSize > 1) and (II == 0) ) {
                // For now, just split in time (assumed to be the first dimension)
                my_count = ((int)count[II]) / wSize;
                overflow = (int)( count[II] - my_count * wSize );

                start[II] = (size_t) (   
                          std::min(wRank,            overflow) * (my_count + 1)
                        + std::max(wRank - overflow, 0       ) *  my_count
                        );

                // Distribute the remainder over the first chunk of processors
                if (wRank < overflow) { my_count++; }
                count[II] = (size_t) my_count;
            }
        }
        num_pts *= count[II];

        if (myCounts != NULL) { myCounts->at(II) = (int) count[II]; }
        if (myStarts != NULL) { myStarts->at(II) = (int) start[II]; }
    }
    #if DEBUG >= 2
    if (wRank == 0) { fprintf(stdout, "\n"); }
    #endif

    #if DEBUG >= 2
    if (num_dims > 1) {
        if (myCounts != NULL) {
            fprintf(stdout, "  Rank %d: starts = %d %d %d %d\n", wRank, 
                    myStarts->at(0), myStarts->at(1), myStarts->at(2), myStarts->at(3));
            fprintf(stdout, "  Rank %d: counts = %d %d %d %d\n", wRank,
                    myCounts->at(0), myCounts->at(1), myCounts->at(2), myCounts->at(3));
        } else {
            fprintf(stdout, "  Rank %d: starts = %zu %zu %zu %zu\n", wRank, 
                    start[0], start[1], start[2], start[3]);
            fprintf(stdout, "  Rank %d: counts = %zu %zu %zu %zu\n", wRank,
                    count[0], count[1], count[2], count[3]);
        }
    }
    #endif

    // Now resize the vector to the appropriate size
    var.resize(num_pts);

    // Now read in the data
    retval = nc_get_vara_double(ncid, var_id, start, count, &var[0]);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

    // Apply scale factor if appropriate
    double scale = 1.;
    retval = nc_get_att_double(ncid, var_id, "scale_factor", &scale);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }
    if (scale != 1.) {
        #if DEBUG >= 2
        if (wRank == 0) { fprintf(stdout, "  scale factor = %g\n", scale); }
        #endif
        for (size_t II = 0; II < num_pts; II++) { var.at(II) = var.at(II) * scale; }
    }

    // Apply offset if appropriate
    double offset = 0.;
    retval = nc_get_att_double(ncid, var_id, "add_offset", &offset);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }
    if (offset != 0.) {
        #if DEBUG >= 2
        if (wRank == 0) { fprintf(stdout, "  additive offset = %g\n", offset); }
        #endif
        for (size_t II = 0; II < num_pts; II++) { var.at(II) = var.at(II) + offset; }
    }

    // Determine masking, if desired
    double fill_val = 1e100;  // backup value
    int num_land = 0;
    int num_water = 0;
    if ( (mask != NULL) and (num_dims == 4) ) {

        mask->resize(var.size());

        // Get the relevant fill value
        nc_get_att_double(ncid, var_id, "_FillValue", &fill_val);
        if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

        #if DEBUG >= 2
        if (wRank == 0) { fprintf(stdout, "  fill value = %g\n", fill_val); }
        #endif

        // If we're at 99% of the fill_val, call it land
        for (size_t II = 0; II < mask->size(); II++) {
            if (fabs(var.at(II)) > 0.99 * fabs(fill_val*scale)) {
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

    #if DEBUG >= 2
    if (wRank == 0) { fprintf(stdout, "\n\n"); }
    #endif

    MPI_Barrier(comm);
    retval = nc_close(ncid);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }
}
