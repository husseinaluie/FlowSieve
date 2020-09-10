#include <vector>
#include <string>
#include <mpi.h>
#include <math.h>
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

    int wRank, wSize;
    MPI_Comm_rank( comm, &wRank );
    MPI_Comm_size( comm, &wSize );

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

    std::vector<signed short> reduced_field; 
    std::vector<double> output_field;
    int index;
    double fmin, fmax, fmiddle, frange = 0, 
           fmin_loc = 1e100, fmax_loc = -1e100, add_offset, scale_factor;
    double max_val =   constants::fill_value < 0 
                     ? constants::fill_value + 2 
                     : constants::fill_value - 2;
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
                ncid, field_varid, "scale_factor", NC_DOUBLE, 1, &scale_factor);
        if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

        retval = nc_put_att_double(
                ncid, field_varid, "add_offset",   NC_DOUBLE, 1, &add_offset);
        if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

        retval = nc_put_vara_short(ncid, field_varid, start, count, &(reduced_field[0]));
        if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

    } else {

        // Get the median value to let us use offsets
        #pragma omp parallel \
        default(none) shared(field, mask) private(index) \
        reduction(max : fmax_loc) reduction(min : fmin_loc)
        {
            #pragma omp for collapse(1) schedule(guided)
            for (index = 0; index < (int) field.size(); index++) {
                if ( (mask == NULL) or ( mask->at(index) ) ) {
                    fmax_loc = std::max(fmax_loc, field.at(index));
                    fmin_loc = std::min(fmin_loc, field.at(index));
                }
            }
        }
        MPI_Allreduce(&fmax_loc, &fmax, 1, MPI_DOUBLE, MPI_MAX, comm);
        MPI_Allreduce(&fmin_loc, &fmin, 1, MPI_DOUBLE, MPI_MIN, comm);

        fmiddle = 0.5 * ( fmax + fmin );
        frange  = fmax - fmin;

        #if DEBUG >= 2
        if (wRank == 0) { 
            fprintf(stdout, "    fmin, fmax, fmiddle, frange = %'g, %'g, %'g, %'g\n", 
                    fmin, fmax, fmiddle, frange);
            fflush(stdout);
        }
        #endif

        scale_factor = std::fabs( frange / max_val );
        if (scale_factor < 1e-25) { scale_factor = 1.; }
        if (scale_factor > 1e+15) { scale_factor = 1e+15; }

        retval = nc_put_att_double(
                ncid, field_varid, "scale_factor", NC_DOUBLE, 1, &scale_factor);
        if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

        retval = nc_put_att_double(ncid, field_varid, 
                        "add_offset", NC_DOUBLE, 1, &fmiddle);
        if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

        output_field.resize(field.size());
        #pragma omp parallel \
        default(none) shared(output_field, field, mask, fmiddle, scale_factor) \
        private(index)
        {
            #pragma omp for collapse(1) schedule(static)
            for (index = 0; index < (int) field.size(); index++) {
                if ( (mask == NULL) or ( mask->at(index) ) ) {
                    output_field.at(index) = ( field.at(index) - fmiddle ) / scale_factor;
                } else {
                    output_field.at(index) = constants::fill_value;
                }
            }
        }

        // Otherwise, just write the 32-bit float field to the file
        retval = nc_put_vara_double(ncid, field_varid, start, count, &(output_field[0]));
        if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

    }

    // Close the file
    MPI_Barrier(comm);
    retval = nc_close(ncid);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

    #if DEBUG >= 1
    if (wRank == 0) { 
        fprintf(stdout, "  - wrote %s to %s -\n", field_name.c_str(), filename.c_str());
        fflush(stdout);
    }
    #endif
}
