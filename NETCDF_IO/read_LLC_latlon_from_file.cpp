
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
 *  @param[in]      comm                the MPI communicator world
 *
 */

void read_LLC_latlon_from_file(
        std::vector<double> &var,
        const std::string & var_name,
        const std::string & filename,
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

    for (int II = 0; II < num_dims; II++) {
        start[II] = 0;
        retval = nc_inq_dim(ncid, dim_ids[II] , NULL, &count[II]);
        if (retval != NC_NOERR ) { NC_ERR(retval, __LINE__, __FILE__); }
        #if DEBUG >= 2
        if (wRank == 0) { fprintf(stdout, "%'zu ", count[II]); }
        #endif

        num_pts *= count[II];
    }
    #if DEBUG >= 2
    if (wRank == 0) { fprintf(stdout, "\n"); }
    fflush(stdout);
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
    
    // Determine masking, if desired
    double fill_val = 1e100;  // backup value
    double var_max = -1e10, var_min = 1e10;
    
    // Get the relevant fill value
    nc_get_att_double(ncid, var_id, "_FillValue", &fill_val);
    if (retval != NC_NOERR ) { NC_ERR(retval, __LINE__, __FILE__); }

    #if DEBUG >= 2
    if (wRank == 0) { fprintf(stdout, "  fill value = %'g\n", fill_val); }
    #endif

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
        var_max = std::max( var_max, var.at(II) );
        var_min = std::min( var_min, var.at(II) );

        // Apply scale factor and offset to non-masked values
        if (scale  != 1.) { var.at(II) = var.at(II) * scale; }
        if (offset != 0.) { var.at(II) = var.at(II) + offset; }
    }

    var_max = var_max * scale + offset;
    var_min = var_min * scale + offset;

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
