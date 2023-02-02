#include <vector>
#include <string>
#include <mpi.h>
#include <math.h>
#include "../netcdf_io.hpp"
#include "../constants.hpp"

void write_integral_to_post(
        const std::vector<
            std::vector<double> > & field,  /**< [in] data to be written to the file*/
        std::string field_name,             /**< [in] name of the variable in the netcdf file */
        std::string field_suffix,           /**< [in] name of the variable in the netcdf file */
        size_t * start,                     /**< [in] starting indices for the write */
        size_t * count,                     /**< [in] size of the write in each dimension */
        const char * filename,              /**< [in] name of the netcdf file */
        const int region_dim,               /**< [in] integer index for region dimension */
        const MPI_Comm comm                 /**< [in] MPI Communicator */
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

    //
    //// Get the median value to let us use offsets
    //
    double fmax, fmin, fmax_loc, fmin_loc;
    size_t index, Iregion;
    const size_t num_regions = field.size(),
                 len_fields  = field.front().size();
    #pragma omp parallel \
    default(none) shared(field) private(Iregion, index) \
    reduction(max : fmax_loc) reduction(min : fmin_loc)
    {
        #pragma omp for collapse(2) schedule(static)
        for (Iregion = 0; Iregion < num_regions; Iregion++) {
            for (index = 0; index < len_fields; index++) {
                fmax_loc = std::max(fmax_loc, field.at(Iregion).at(index));
                fmin_loc = std::min(fmin_loc, field.at(Iregion).at(index));
            }
        }
    }
    MPI_Allreduce(&fmax_loc, &fmax, 1, MPI_DOUBLE, MPI_MAX, comm);
    MPI_Allreduce(&fmin_loc, &fmin, 1, MPI_DOUBLE, MPI_MIN, comm);

    const double fmiddle = 0.5 * (fmax + fmin);
    const double frange  = fmax - fmin;
    double scale_factor;

    #if DEBUG >= 2
    if (wRank == 0) { 
        fprintf(stdout, "    fmin, fmax, fmiddle, frange = %'g, %'g, %'g, %'g\n", 
                fmin, fmax, fmiddle, frange);
        fflush(stdout);
    }
    #endif

    std::vector< std::vector< signed short > > int_fields;
    std::vector< std::vector< float        > > float_fields;
    std::vector< std::vector< double       > > double_fields;

    // Get the variable ID for the field
    int field_varid;
    retval = nc_inq_varid(ncid, (field_name + field_suffix).c_str(), &field_varid );
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

    if (constants::CAST_TO_INT) {

        double local_double;
        signed short local_int;

        // Number of Discrete Representable Values
        //   (less two for numerical reasons)
        const int ndrv =   constants::fill_value_s < 0 
                         ? constants::fill_value_s + 2 
                         : constants::fill_value_s - 2;

        scale_factor = frange / ndrv;

        int_fields.resize( field.size(), std::vector<signed short>(field.front().size(), 0.) );

        #pragma omp parallel \
        default(none) shared(int_fields, field) \
        private(Iregion, index, local_double, local_int)
        {
            #pragma omp for collapse(2) schedule(static)
            for (Iregion = 0; Iregion < num_regions; Iregion++) {
                for (index = 0; index < len_fields; index++) {
                    local_double = (field.at(Iregion).at(index) - fmiddle) / frange;
                    local_int = (signed short) round(ndrv * local_double);
                    int_fields.at(Iregion).at(index) = local_int;
                }
            }
        }

        // We need to record the scale and translation used to encode in signed shorts
        retval = nc_put_att_double(ncid, field_varid, "scale_factor", NC_DOUBLE, 1, &scale_factor);
        if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

        retval = nc_put_att_double(ncid, field_varid, "add_offset",   NC_DOUBLE, 1, &fmiddle);
        if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

        // Now write the reduced fields
        for (size_t Iregion = 0; Iregion < field.size(); ++Iregion) {
            start[region_dim] = Iregion;
            retval = nc_put_vara_short(ncid, field_varid, start, count, &(int_fields.at(Iregion)[0]));
            if (retval) { NC_ERR(retval, __LINE__, __FILE__); }
        }

    } else if (constants::CAST_TO_SINGLE) {

        const double max_val =   constants::fill_value < 0 
                               ? constants::fill_value + 2 
                               : constants::fill_value - 2;

        //scale_factor = std::fabs( frange / max_val );
        scale_factor = fabs( frange / max_val );
        if (scale_factor < 1e-25) { scale_factor = 1.; }
        if (scale_factor > 1e+15) { scale_factor = 1e+15; }

        retval = nc_put_att_double(ncid, field_varid, "scale_factor", NC_DOUBLE, 1, &scale_factor);
        if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

        retval = nc_put_att_double(ncid, field_varid, "add_offset", NC_DOUBLE, 1, &fmiddle);
        if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

        float_fields.resize( field.size(), std::vector<float>(field.front().size(), 0.) );
        #pragma omp parallel \
        default(none) shared(float_fields, field, scale_factor) \
        private(Iregion, index)
        {
            #pragma omp for collapse(2) schedule(static)
            for (Iregion = 0; Iregion < num_regions; Iregion++) {
                for (index = 0; index < len_fields; index++) {
                    float_fields.at(Iregion).at(index) = ( field.at(Iregion).at(index) - fmiddle ) / scale_factor;
                }
            }
        }

        // Write the integral for each region
        for (size_t Iregion = 0; Iregion < field.size(); ++Iregion) {
            start[region_dim] = Iregion;
            retval = nc_put_vara_float( ncid, field_varid, start, count, &(float_fields.at(Iregion)[0]));
            if (retval) { NC_ERR(retval, __LINE__, __FILE__); }
        }

    } else {

        const double max_val =   constants::fill_value < 0 
                               ? constants::fill_value + 2 
                               : constants::fill_value - 2;

        //scale_factor = std::fabs( frange / max_val );
        scale_factor = fabs( frange / max_val );
        if (scale_factor < 1e-25) { scale_factor = 1.; }
        if (scale_factor > 1e+15) { scale_factor = 1e+15; }

        retval = nc_put_att_double(ncid, field_varid, "scale_factor", NC_DOUBLE, 1, &scale_factor);
        if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

        retval = nc_put_att_double(ncid, field_varid, "add_offset", NC_DOUBLE, 1, &fmiddle);
        if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

        double_fields.resize( field.size(), std::vector<double>(field.front().size(), 0.) );
        #pragma omp parallel \
        default(none) shared(double_fields, field, scale_factor) \
        private(Iregion, index)
        {
            #pragma omp for collapse(2) schedule(static)
            for (Iregion = 0; Iregion < num_regions; Iregion++) {
                for (index = 0; index < len_fields; index++) {
                    double_fields.at(Iregion).at(index) = ( field.at(Iregion).at(index) - fmiddle ) / scale_factor;
                }
            }
        }

        // Write the integral for each region
        for (size_t Iregion = 0; Iregion < field.size(); ++Iregion) {
            start[region_dim] = Iregion;
            retval = nc_put_vara_double( ncid, field_varid, start, count, &(double_fields.at(Iregion)[0]));
            if (retval) { NC_ERR(retval, __LINE__, __FILE__); }
        }

    }

    // Close the file
    MPI_Barrier(comm);
    retval = nc_close(ncid);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

    #if DEBUG >= 1
    if (wRank == 0) { fprintf(stdout, "  - wrote %s to %s -\n", 
            (field_name + field_suffix).c_str(), filename); }
    #endif
}
