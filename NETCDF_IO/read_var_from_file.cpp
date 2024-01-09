
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
    #if DEBUG >= 1
    // NC_FORMAT_CLASSIC, NC_FORMAT_64BIT_OFFSET, NC_FORMAT_CDF5, NC_FORMAT_NETCDF4, NC_FORMAT_NETCDF4_CLASSIC
    if (wRank == 0) {
        fprintf( stdout, "Detected file format (%d) is ", input_nc_format );
        if ( input_nc_format == NC_FORMAT_CLASSIC )         { fprintf( stdout, "NC_FORMAT_CLASSIC (%d)", NC_FORMAT_CLASSIC ); }
        if ( input_nc_format == NC_FORMAT_64BIT_OFFSET )    { fprintf( stdout, "NC_FORMAT_64BIT_OFFSET (%d)", NC_FORMAT_64BIT_OFFSET ); }
        if ( input_nc_format == NC_FORMAT_CDF5 )            { fprintf( stdout, "NC_FORMAT_CDF5 (%d)", NC_FORMAT_CDF5 ); }
        if ( input_nc_format == NC_FORMAT_NETCDF4 )         { fprintf( stdout, "NC_FORMAT_NETCDF4 (%d)", NC_FORMAT_NETCDF4 ); }
        if ( input_nc_format == NC_FORMAT_NETCDF4_CLASSIC ) { fprintf( stdout, "NC_FORMAT_NETCDF4_CLASSIC (%d)", NC_FORMAT_NETCDF4_CLASSIC ); }
        fprintf( stdout, " \n" );
        fflush(stdout);
    }
    #endif
    assert( input_nc_format == NC_FORMAT_NETCDF4 ); // input file must be netCDF-4 format. Use `nccopy -k netCDF-4 input.nc output.nc` to change file version

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
    #if DEBUG >= 3
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
                else                { Nprocs_in_dim = 0; assert(false); }  // II <= 1 so won't happen

                my_count = ( (int)count[II] ) / Nprocs_in_dim;
                overflow = (int)( count[II] - my_count * Nprocs_in_dim );

                Index1to4( wRank, Itime_proc,      Idepth_proc,     Ilat_proc, Ilon_proc,
                                  Nprocs_in_time,  Nprocs_in_depth, 1,         1          );
                if      ( II == 0 ) { Iproc_in_dim = Itime_proc;  }
                else if ( II == 1 ) { Iproc_in_dim = Idepth_proc; }
                else                { Iproc_in_dim = -1; assert(false); }  // II <= 1 so won't happen

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
    #endif

    // Now resize the vector to the appropriate size
    var.resize(num_pts);

    // Now read in the data
    //const size_t MAX_PTS_PER_READ = 4 * pow(512,3);
    const size_t MAX_PTS_PER_READ = 2 * pow(512,3);
    if (num_pts < MAX_PTS_PER_READ) {
        retval = nc_get_vara_double(ncid, var_id, start, count, &var[0]);
        if (retval != NC_NOERR ) { NC_ERR(retval, __LINE__, __FILE__); }
    } else {
        #if DEBUG >= 2
        if (wRank == 0) { fprintf(stdout, "Data is large, so will read in chunks\n"); }
        #endif
        // We can't read in more than 4 * 512^3 at a time, so break into chunks 
        //   right now we assume that splitting the first dimension is sufficient
        //   Surprise! Just splitting the first dimension is no longer sufficient.
        //   Time to expand this feature.
        //
        //   So, figure out what dimensions we need to split in
        unsigned int num_chunks_time, num_chunks_depth, num_chunks_dim3;
        if ( ( num_pts / count[0] ) < MAX_PTS_PER_READ  ) {
            // If splitting in just time is good enough
            num_chunks_time = std::min( count[0], 1 + (size_t)ceil( num_pts / MAX_PTS_PER_READ ) );
            num_chunks_depth = 1;
            num_chunks_dim3 = 1;
        } else if ( ( num_dims > 1 ) and ( num_pts / (count[0]*count[1]) ) < MAX_PTS_PER_READ ) {
            // If splitting in just time and depth is good enough
            num_chunks_time = count[0];
            num_chunks_depth = std::min( count[1], 1 + (size_t) ceil( (num_pts / num_chunks_time) / MAX_PTS_PER_READ ) );
            num_chunks_dim3 = 1;
        } else if ( num_dims > 2 ) {
            // Otherwise, split in 'lat' too
            num_chunks_time = count[0];
            num_chunks_depth = count[1];
            num_chunks_dim3 = std::min( count[2], 1 + (size_t) ceil( (num_pts / (num_chunks_time*num_chunks_depth)) / MAX_PTS_PER_READ ) );
        } else {
            fprintf( stderr, "This data is too large to be read with current implementations.\n" );
            assert(false);
        }

        #if DEBUG >= 2
        fprintf( stdout, "Using %d, %d, %d chunks in time, depth, and dim3.\n", num_chunks_time, num_chunks_depth, num_chunks_dim3 );
        #endif

        size_t pts_per_slice, target = 0;
        size_t start_chunk[num_dims], count_chunk[num_dims];
        size_t curr_time_count, curr_depth_count, curr_dim3_count;
        start_chunk[0] = 0;
        for (unsigned int Ichunk_time = 0; Ichunk_time < num_chunks_time; ++Ichunk_time) {
            curr_time_count = (count[0] / num_chunks_time) + ( (Ichunk_time < (count[0] % num_chunks_time)) ? 1 : 0 );
            pts_per_slice = curr_time_count;
            count_chunk[0] = curr_time_count;

            if (num_dims > 1) { start_chunk[1] = 0; }

            for (unsigned int Ichunk_depth = 0; Ichunk_depth < num_chunks_depth; ++Ichunk_depth) {
                if ( num_dims > 1) {
                    curr_depth_count = (count[1] / num_chunks_depth) + ( (Ichunk_depth < (count[1] % num_chunks_depth)) ? 1 : 0 );
                    pts_per_slice *= curr_depth_count;
                    count_chunk[1] = curr_depth_count;
                }

                if (num_dims > 2) { start_chunk[2] = 0; }
                if (num_dims > 3) { start_chunk[3] = 0; }

                for (unsigned int Ichunk_dim3 = 0; Ichunk_dim3 < num_chunks_dim3; ++Ichunk_dim3) {
                    if ( num_dims > 2 ) {
                        curr_dim3_count = (count[2] / num_chunks_dim3) + ( (Ichunk_dim3 < (count[2] % num_chunks_dim3)) ? 1 : 0 );
                        pts_per_slice *= curr_dim3_count;
                        count_chunk[2] = curr_dim3_count;
                    }

                    if (num_dims > 3) {
                        pts_per_slice *= count[3];
                        count_chunk[3] = count[3];
                    }
                    
                    #if DEBUG>=2
                    if (num_dims > 2) {
                        fprintf( stdout, "Loading chunk (%d, %d, %d) with starts (%zu, %zu, %zu) and counts (%zu, %zu, %zu) [%zu points total, with target %zu].\n", 
                            Ichunk_time, Ichunk_depth, Ichunk_dim3, 
                            start_chunk[0], start_chunk[1], start_chunk[2],
                            count_chunk[0], count_chunk[1], count_chunk[2], 
                            pts_per_slice, target );
                    } else if (num_dims > 1) {
                        fprintf( stdout, "Loading chunk (%d, %d) with starts (%zu, %zu) and counts (%zu, %zu) [%zu points total, with target %zu].\n", 
                            Ichunk_time, Ichunk_depth,
                            start_chunk[0], start_chunk[1],
                            count_chunk[0], count_chunk[1],
                            pts_per_slice, target );
                    } else {
                        fprintf( stdout, "Loading chunk %d with start %zu and counts %zu [%zu points total, with target %zu].\n", 
                            Ichunk_time,
                            start_chunk[0],
                            count_chunk[0],
                            pts_per_slice, target );
                    }
                    #endif

                    // do reading
                    retval = nc_get_vara_double(ncid, var_id, start_chunk, count_chunk, &var[target]);
                    if (retval != NC_NOERR ) { NC_ERR(retval, __LINE__, __FILE__); }
                    target += pts_per_slice;

                    
                    if (num_dims > 2) { start_chunk[2] += curr_dim3_count; }
                }

                if (num_dims > 1) { start_chunk[1] += curr_depth_count; }
            }

            start_chunk[0] += curr_time_count;
        }
    }
    
    // Determine masking, if desired
    double fill_val = 1e100;  // backup value
    size_t num_land = 0, num_water = 0, num_unmasked = 0;

    if (mask != NULL) { mask->resize(var.size()); }

    // Get the relevant fill value
    nc_get_att_double(ncid, var_id, "_FillValue", &fill_val);
    if (retval != NC_NOERR ) { NC_ERR(retval, __LINE__, __FILE__); }

    #if DEBUG >= 2
    if (wRank == 0) { fprintf(stdout, "  fill value = %'g\n", fill_val); }
    #endif
    double var_max = -std::fabs(fill_val), var_min = std::fabs(fill_val);

    // Apply scale factor if appropriate
    double scale = 1.;
    retval = nc_get_att_double(ncid, var_id, "scale_factor", &scale);
    if (retval != NC_NOERR ) { NC_ERR(retval, __LINE__, __FILE__); }
    #if DEBUG >= 2
    if (wRank == 0) { fprintf(stdout, "  scale factor = %'g\n", scale); }
    #endif

    // Apply offset if appropriate
    double offset = 0.;
    retval = nc_get_att_double(ncid, var_id, "add_offset", &offset);
    if (retval != NC_NOERR ) { NC_ERR(retval, __LINE__, __FILE__); }
    #if DEBUG >= 2
    if (wRank == 0) { fprintf(stdout, "  additive offset = %'g\n", offset); }
    #endif

    // Masked if equal to fill value
    for (size_t II = 0; II < var.size(); II++) {
        if ( var.at(II) == fill_val ) {
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

            // Apply scale factor and offset to non-masked values
            if (scale  != 1.) { var.at(II) = var.at(II) * scale; }
            if (offset != 0.) { var.at(II) = var.at(II) + offset; }
        }
    }

    var_max = var_max * scale + offset;
    var_min = var_min * scale + offset;

    #if DEBUG >= 1
    if (wRank == 0) {
        fprintf(stdout, "  Land cover = %'.4g%% (%'zu water vs %'zu land) (%'zu land converted to water) \n", 
                100 * ((double)num_land) / (num_land + num_water + num_unmasked),
                num_water + num_unmasked, num_land, num_unmasked);
    }
    #endif

    #if DEBUG >= 1
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
