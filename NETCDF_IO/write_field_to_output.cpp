#include <fenv.h>
#include <vector>
#include <string>
#include <mpi.h>
#include <math.h>
#include <cassert>
#include "../netcdf_io.hpp"
#include "../constants.hpp"

void write_field_to_output(
        const std::vector<double> & field,
        const std::string & field_name,
        const size_t * start,
        const size_t * count,
        const std::string & filename,
        const std::vector<bool> * mask,
        MPI_Comm comm
        ) {

    // During writing, ignore floating point exceptions
    std::fenv_t fe_env;
    feholdexcept( &fe_env );

    int wRank, wSize;
    MPI_Comm_rank( comm, &wRank );
    MPI_Comm_size( comm, &wSize );

    #if DEBUG >= 2
    if (wRank == 0) { fprintf(stdout, "  Preparing to write %s to %s (%'zu points).\n", field_name.c_str(), filename.c_str(), field.size()); }
    fflush(stdout);
    #endif

    // Open the NETCDF file
    int FLAG = NC_NETCDF4 | NC_WRITE | NC_MPIIO;
    int ncid=0, retval;
    char buffer [50];
    snprintf(buffer, 50, filename.c_str());
    MPI_Barrier(comm);
    retval = nc_open_par(buffer, FLAG, comm, MPI_INFO_NULL, &ncid);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

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
    
    size_t num_pts = 1;
    for (int Idim = 0; Idim < num_dims; Idim++) {
        num_pts *= count[Idim];
    }

    std::vector<signed short> reduced_field; 
    std::vector<double> output_field;
    size_t index;
    double add_offset, scale_factor;

    // This is the maximum value of the transformed variable
    const double max_val =
        constants::CAST_TO_SINGLE ?
            ( constants::fill_value < 0 ? 
              constants::fill_value + 2 : 
              constants::fill_value - 2 ) :
            ( constants::fill_value_double < 0 ? 
              constants::fill_value_double + 2 : 
              constants::fill_value_double - 2 ) ;

    if (constants::CAST_TO_INT) {
        // If we want to reduce output size, pack into short ints
        //   floats are 32bit, short ints are 16bit, so we can cut
        //   file size in half. Of course, this is at the cost
        //   of precision.

        reduced_field.resize(field.size());

        #if DEBUG >= 2
        if (wRank == 0) { fprintf(stdout, "    Preparing to package the field.\n"); }
        fflush(stdout);
        #endif
        package_field(reduced_field, scale_factor, add_offset, field, mask);

        // We need to record the scale and translation used to encode in signed shorts
        retval = nc_put_att_double( ncid, field_varid, "scale_factor", NC_DOUBLE, 1, &scale_factor );
        if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

        retval = nc_put_att_double( ncid, field_varid, "add_offset",   NC_DOUBLE, 1, &add_offset );
        if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

        retval = nc_put_vara_short( ncid, field_varid, start, count, &(reduced_field[0]) );
        if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

    } else {

        // Get the median value to let us use offsets
        //      to do this, get min and max (first, MPI_local values)
        //      initalize with first element before looping over all values
        double  fmax_loc = 0,
                fmin_loc = 0;
        #pragma omp parallel \
        default(none) shared(field, mask) private(index) \
        reduction(max : fmax_loc) reduction(min : fmin_loc)
        {
            #pragma omp for collapse(1) schedule(dynamic)
            for (index = 0; index < field.size(); index++) {
                if ( (mask == NULL) or ( mask->at(index) ) ) {
                    fmax_loc = std::max(fmax_loc, field.at(index));
                    fmin_loc = std::min(fmin_loc, field.at(index));
                }
            }
        }

        // Now communicate across MPI to get true max and min values
        double fmin, fmax;
        MPI_Allreduce(&fmax_loc, &fmax, 1, MPI_DOUBLE, MPI_MAX, comm);
        MPI_Allreduce(&fmin_loc, &fmin, 1, MPI_DOUBLE, MPI_MIN, comm);

        const double fmiddle = 0.5 * ( fmax + fmin );
        const double frange  = fmax - fmin;

        #if DEBUG >= 2
        if (wRank == 0) { 
            fprintf(stdout, "    fmin, fmax, fmiddle, frange = %'g, %'g, %'g, %'g\n", 
                    fmin, fmax, fmiddle, frange);
            fflush(stdout);
        }
        #endif

        // Get the multiplicative scale factor. If it's extreme, then truncate it.
        scale_factor = frange == 0. ? 1. : fabs( frange / max_val );

        #if DEBUG >= 2
        if (wRank == 0) { 
            fprintf(stdout, "    scale_factor, add_offtset, max_transform_val = %'g, %'g, %'g\n", scale_factor, fmiddle, max_val);
            fflush(stdout);
        }
        #endif

        retval = nc_put_att_double( ncid, field_varid, "scale_factor", NC_DOUBLE, 1, &scale_factor );
        if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

        retval = nc_put_att_double(ncid, field_varid, "add_offset", NC_DOUBLE, 1, &fmiddle );
        if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

        output_field.resize(field.size());
        #pragma omp parallel \
        default(none) \
        shared(output_field, field, mask, scale_factor) \
        private(index) \
        firstprivate( fmiddle )
        {
            #pragma omp for collapse(1) schedule(static)
            for (index = 0; index < field.size(); index++) {
                if ( (mask == NULL) or ( mask->at(index) ) ) {
                    output_field.at(index) = ( field.at(index) - fmiddle ) / scale_factor;
                } else {
                    output_field.at(index) = constants::CAST_TO_SINGLE ?
                                                constants::fill_value :
                                                constants::fill_value_double;
                }
            }
        }

        const size_t MAX_PTS_PER_READ = 2 * pow(512,3);
        if (num_pts < MAX_PTS_PER_READ) {

            retval = nc_put_vara_double(ncid, field_varid, start, count, &(output_field[0]));
            if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

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

                #if DEBUG >= 2
                fprintf( stdout, "Using %d, chunks in dim1.\n", num_chunks_time );
                #endif
            } else if ( ( num_dims > 1 ) and ( num_pts / (count[0]*count[1]) ) < MAX_PTS_PER_READ ) {
                // If splitting in just time and depth is good enough
                num_chunks_time = count[0];
                num_chunks_depth = std::min( count[1], 1 + (size_t) ceil( (num_pts / num_chunks_time) / MAX_PTS_PER_READ ) );
                num_chunks_dim3 = 1;

                #if DEBUG >= 2
                fprintf( stdout, "Using %d, %d, chunks in dim1 and dim2.\n", num_chunks_time, num_chunks_depth );
                #endif
            } else if ( num_dims > 2 ) {
                // Otherwise, split in 'lat' too
                num_chunks_time = count[0];
                num_chunks_depth = count[1];
                num_chunks_dim3 = std::min( count[2], 1 + (size_t) ceil( (num_pts / (num_chunks_time*num_chunks_depth)) / MAX_PTS_PER_READ ) );

                #if DEBUG >= 2
                fprintf( stdout, "Using %d, %d, %d, chunks in dim1, dim2, and dim3.\n", num_chunks_time, num_chunks_depth, num_chunks_dim3 );
                #endif
            } else {
                fprintf( stderr, "This data is too large to be written with current implementations.\n" );
                assert(false);
            }

            size_t pts_per_slice, target = 0;
            size_t start_chunk[num_dims], count_chunk[num_dims];
            size_t curr_time_count, curr_depth_count, curr_dim3_count;
            start_chunk[0] = start[0];
            for (unsigned int Ichunk_time = 0; Ichunk_time < num_chunks_time; ++Ichunk_time) {
                curr_time_count = (count[0] / num_chunks_time) + ( (Ichunk_time < (count[0] % num_chunks_time)) ? 1 : 0 );
                pts_per_slice = curr_time_count;
                count_chunk[0] = curr_time_count;

                if (num_dims > 1) { start_chunk[1] = start[1]; }

                for (unsigned int Ichunk_depth = 0; Ichunk_depth < num_chunks_depth; ++Ichunk_depth) {
                    if ( num_dims > 1) {
                        curr_depth_count = (count[1] / num_chunks_depth) + ( (Ichunk_depth < (count[1] % num_chunks_depth)) ? 1 : 0 );
                        pts_per_slice *= curr_depth_count;
                        count_chunk[1] = curr_depth_count;
                    }

                    if (num_dims > 2) { start_chunk[2] = start[2]; }
                    if (num_dims > 3) { start_chunk[3] = start[3]; }

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

                        // do writing
                        retval = nc_put_vara_double(ncid, field_varid, start_chunk, count_chunk, &(output_field[target]));
                        if (retval) { NC_ERR(retval, __LINE__, __FILE__); }
                        target += pts_per_slice;


                        if (num_dims > 2) { start_chunk[2] += curr_dim3_count; }
                    }

                    if (num_dims > 1) { start_chunk[1] += curr_depth_count; }
                }

                start_chunk[0] += curr_time_count;
            }
        }

    }

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
