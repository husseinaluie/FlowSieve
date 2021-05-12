
#include "../netcdf_io.hpp"
#include "../constants.hpp"
#include "../functions.hpp"
#include <cassert>
#include <math.h>

/*!
 *  \brief Read a specific variable from a specific file.
 *
 *  Accounts for variable attributes 'scale_factor' and 'add_offset'.
 *
 *  If mask != NULL, then determine the mask based on variable attribute '_FillValue'
 *
 *  @param[in,out]  var                 vector into which to store the loaded variable
 *  @param[in]      var_name            name of the variable to be read
 *  @param[in]      filename            name of the file from which to load the variable
 *  @param[in,out]  mask                point to where a mask array should be stored (if not NULL) (true = water, false = land)
 *  @param[in,out]  myCounts            the sizes of each dimension (on this MPI process) if not NULL
 *  @param[in,out]  myStarts            the starting index for each dimension, if not NULL
 *  @param[in]      Nprocs_in_time      Number of MPI processors in time (for dividing data)
 *  @param[in]      Nprocs_in_depth     Number of MPI processors in depth (for dividing data)
 *  @param[in]      do_splits           boolean indicating if the arrays should be split over MPI procs.
 *  @param[in]      force_split_dim     Dimension along which data splitting should be force
 *  @param[in]      land_fill_value     Value to place at 'land' areas, if needed
 *  @param[in]      comm                the MPI communicator world
 *
 */

void read_var_from_file(
        std::vector<double> &var,
        const std::string & var_name,
        const std::string & filename,
        std::vector<bool> *mask,
        std::vector<int> *myCounts,
        std::vector<int> *myStarts,
        const int Nprocs_in_time,
        const int Nprocs_in_depth,
        const bool do_splits,
        const int force_split_dim,
        const double land_fill_value,
        const MPI_Comm comm
        ) {

    assert( check_file_existence( filename.c_str() ) );

    int wRank, wSize, Nprocs_in_dim, Iproc_in_dim;
    MPI_Comm_rank( comm, &wRank );
    MPI_Comm_size( comm, &wSize );

    // Open the NETCDF file
    const int str_len = 100;
    int FLAG = NC_NETCDF4 | NC_MPIIO;
    int ncid=0, retval;
    char buffer [str_len];
    snprintf(buffer, str_len, filename.c_str());

    #if DEBUG >= 1
    if (wRank == 0) {
        fprintf(stdout, "Attempting to read %s from %s\n", var_name.c_str(), buffer);
        fflush(stdout);
    }
    #endif

    retval = nc_open_par(buffer, FLAG, comm, MPI_INFO_NULL, &ncid);
    if (retval != NC_NOERR ) { NC_ERR(retval, __LINE__, __FILE__); }

    // Check if netcdf-4 format
    int input_nc_format;
    retval = nc_inq_format( ncid, &input_nc_format );
    if (retval != NC_NOERR ) { NC_ERR(retval, __LINE__, __FILE__); }
    //assert( input_nc_format == NC_FORMAT_NETCDF4 ); // input file must be netCDF-4 format

    char varname [str_len];
    snprintf(varname, str_len, var_name.c_str());

    int var_id, num_dims;
    int dim_ids[NC_MAX_VAR_DIMS];

    // Get the ID for the variable
    retval = nc_inq_varid(ncid, varname, &var_id );
    if (retval != NC_NOERR ) { NC_ERR(retval, __LINE__, __FILE__); }

    // This should return an error if the variable doesn't exist
    retval = nc_inq_var(ncid, var_id, NULL, NULL, NULL, NULL, NULL);
    if (retval != NC_NOERR ) { NC_ERR(retval, __LINE__, __FILE__); }
    if (retval == NC_ENOTVAR ) { NC_ERR(NC_ENOTVAR, __LINE__, __FILE__); }

    // Get information about the variable
    retval = nc_inq_var(ncid, var_id, NULL, NULL, &num_dims, dim_ids, NULL );
    if (retval != NC_NOERR ) { NC_ERR(retval, __LINE__, __FILE__); }
    assert( num_dims > 0 );
    #if DEBUG >= 2
    if (wRank == 0) {
        if (num_dims == 1) {
            fprintf(stdout, "  has %'d dimension of size ", num_dims);
        } else {
            fprintf(stdout, "  has %'d dimensions of size ", num_dims);
        }
    }
    #endif

    // Get the size of each dimension
    size_t start[num_dims], count[num_dims];
    size_t num_pts = 1;
    int my_count, overflow,
        Itime_proc, Idepth_proc, Ilat_proc, Ilon_proc;
    if (myCounts != NULL) {
        myCounts->resize(num_dims);
        myStarts->resize(num_dims);
    }
    for (int II = 0; II < num_dims; II++) {
        start[II] = 0;
        retval = nc_inq_dim(ncid, dim_ids[II] , NULL, &count[II]);
        if (retval != NC_NOERR ) { NC_ERR(retval, __LINE__, __FILE__); }
        #if DEBUG >= 2
        if (wRank == 0) { fprintf(stdout, "%'zu ", count[II]); }
        #endif

        if (do_splits) {
            // If we're split on multiple MPI procs and have > 2 dimensions, 
            //   then divide all but the last two 
            //
            //   we don't split the last two because those 
            //   are assumed to be lat/lon

            if ( ( (num_dims > 2) and (wSize > 1) and (II <= 1) )
                 or
                 ( II == force_split_dim )
               ) {

                assert( Nprocs_in_time > 0 ); // Must specify the number of processors used in time
                assert( Nprocs_in_depth > 0 ); // Must specify the number of processors used in depth
                assert( Nprocs_in_time * Nprocs_in_depth == wSize ); // Total number of processors does no match with specified values

                if      ( II == 0 ) { Nprocs_in_dim = Nprocs_in_time;  }
                else if ( II == 1 ) { Nprocs_in_dim = Nprocs_in_depth; }

                my_count = ( (int)count[II] ) / Nprocs_in_dim;
                overflow = (int)( count[II] - my_count * Nprocs_in_dim );

                Index1to4( wRank, Itime_proc,      Idepth_proc,     Ilat_proc, Ilon_proc,
                                  Nprocs_in_time,  Nprocs_in_depth, 1,         1          );
                if      ( II == 0 ) { Iproc_in_dim = Itime_proc;  }
                else if ( II == 1 ) { Iproc_in_dim = Idepth_proc; }

                start[II] = (size_t) (   
                          std::min(Iproc_in_dim,            overflow) * (my_count + 1)
                        + std::max(Iproc_in_dim - overflow, 0       ) *  my_count
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
    fflush(stdout);
    if (num_dims == 4) { fprintf( stdout, " MPI rank %d has sizes:a %'zu %'zu %'zu %'zu\n", wRank, count[0], count[1], count[2], count[3] ); }
    #endif

    // Now resize the vector to the appropriate size
    var.resize(num_pts);

    // Now read in the data
    if (num_pts < 4 * pow(512,3)) {
        retval = nc_get_vara_double(ncid, var_id, start, count, &var[0]);
        if (retval != NC_NOERR ) { NC_ERR(retval, __LINE__, __FILE__); }
    } else {
        #if DEBUG >= 2
        if (wRank == 0) { fprintf(stdout, "Data is large, so will read in chunks\n"); }
        #endif
        // We can't read in more than 4 * 512^3 at a time, so break into chunks 
        //   right now we assume that splitting the first dimension is sufficient
        const unsigned int num_chunks = ceil( num_pts / ( 4 * pow(512,3) ) );

        size_t target = 0;
        const size_t full_time_count = count[0];
        size_t curr_count, pts_per_slice = 1;


        for (int Idim = 1; Idim < num_dims; ++Idim) {
            pts_per_slice *= count[Idim];
        }

        for (unsigned int Ichunk = 0; Ichunk < num_chunks; ++Ichunk) {

            curr_count = floor( full_time_count / num_chunks );
            if (Ichunk < ( full_time_count % num_chunks ) ) {
                curr_count++;
            }

            count[0] = curr_count;

            retval = nc_get_vara_double(ncid, var_id, start, count, &var[target]);
            if (retval != NC_NOERR ) { NC_ERR(retval, __LINE__, __FILE__); }

            start[0] += curr_count;
            target += curr_count * pts_per_slice;

        }
    }

    // Apply scale factor if appropriate
    double scale = 1.;
    retval = nc_get_att_double(ncid, var_id, "scale_factor", &scale);
    if (retval != NC_NOERR ) { NC_ERR(retval, __LINE__, __FILE__); }
    #if DEBUG >= 2
    if (wRank == 0) { fprintf(stdout, "  scale factor = %'g\n", scale); }
    #endif
    if (scale != 1.) { for (size_t II = 0; II < num_pts; II++) { var.at(II) = var.at(II) * scale; } }

    // Apply offset if appropriate
    double offset = 0.;
    retval = nc_get_att_double(ncid, var_id, "add_offset", &offset);
    if (retval != NC_NOERR ) { NC_ERR(retval, __LINE__, __FILE__); }
    #if DEBUG >= 2
    if (wRank == 0) { fprintf(stdout, "  additive offset = %'g\n", offset); }
    #endif
    if (offset != 0.) { for (size_t II = 0; II < num_pts; II++) { var.at(II) = var.at(II) + offset; } }


    // Determine masking, if desired
    double fill_val = 1e100;  // backup value
    double var_max = -1e10, var_min = 1e10;
    size_t num_land = 0, num_water = 0, num_unmasked = 0;

    if (mask != NULL) { mask->resize(var.size()); }

    // Get the relevant fill value
    nc_get_att_double(ncid, var_id, "_FillValue", &fill_val);
    if (retval != NC_NOERR ) { NC_ERR(retval, __LINE__, __FILE__); }

    #if DEBUG >= 2
    if (wRank == 0) { fprintf(stdout, "  fill value = %'g\n", fill_val); }
    #endif

    // If we're at 99% of the fill_val, call it land
    for (size_t II = 0; II < var.size(); II++) {
        if (fabs(var.at(II)) > 0.999 * (fabs(fill_val*scale + offset))) {
            if (constants::FILTER_OVER_LAND) {
                // If requested to filter over land, then fill in the mask now
                if (mask != NULL) { mask->at(II) = true; }
                num_unmasked++;
                var.at(II) = land_fill_value;
            } else {
                if (mask != NULL) { mask->at(II) = false; }
                num_land++;
            }
        } else {
            var_max = std::max( var_max, var.at(II) );
            var_min = std::min( var_min, var.at(II) );
            if (mask != NULL) { mask->at(II) = true; }
            num_water++;
        }
    }

    #if DEBUG >= 1
    if (wRank == 0) {
        fprintf(stdout, "  Land cover = %'.4g%% (%'zu water vs %'zu land) (%'zu land converted to water) \n", 
                100 * ((double)num_land) / (num_land + num_water + num_unmasked),
                num_water + num_unmasked, num_land, num_unmasked);
    }
    #endif

    #if DEBUG >= 2
    if (wRank == 0) { 
        fprintf(stdout, "  var_max = %g\n", var_max);
        fprintf(stdout, "  var_min = %g\n", var_min);
        fprintf(stdout, "\n\n"); 
    }
    #endif

    MPI_Barrier(comm);
    retval = nc_close(ncid);
    if (retval != NC_NOERR ) { NC_ERR(retval, __LINE__, __FILE__); }
}
