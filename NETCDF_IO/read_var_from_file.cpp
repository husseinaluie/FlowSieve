
#include "../netcdf_io.hpp"
#include "../constants.hpp"
#include <vector>
#include <math.h>

// Write to netcdf file
void read_var_from_file(
        std::vector<double> &var,   /**< [in] Vector into which to store the variable */
        const char * var_name,      /**< [in] Name of the variable*/
        const char * filename,      /**< [in] Name of the file*/
        std::vector<double> *mask   /**< [in] Pointer to mask array to be created*/
        ) {

    #if DEBUG >= 1
    fprintf(stdout, "Attempting to read %s from %s\n", var_name, filename);
    #endif

    // Open the NETCDF file
    int FLAG = NC_NETCDF4 | NC_NOWRITE;
    int ncid=0, retval;
    char buffer [50];
    snprintf(buffer, 50, filename);
    if (( retval = nc_open(buffer, FLAG, &ncid) )) { NC_ERR(retval, __LINE__, __FILE__); }

    char varname [50];
    snprintf(varname, 50, var_name);

    int var_id, num_dims;
    int dim_ids[NC_MAX_VAR_DIMS];

    // Get the ID for the variable
    if ((retval = nc_inq_varid(ncid, varname, &var_id ))) { NC_ERR(retval, __LINE__, __FILE__); }

    // Get information about the variable
    if ((retval = nc_inq_var(ncid, var_id, NULL, NULL, &num_dims, dim_ids, NULL ))) { NC_ERR(retval, __LINE__, __FILE__); }
    #if DEBUG >= 2
    if (num_dims == 1) {
        fprintf(stdout, "  has %d dimension of size ", num_dims);
    } else {
        fprintf(stdout, "  has %d dimensions of size ", num_dims);
    }
    #endif

    // Get the size of each dimension
    size_t start[num_dims], count[num_dims];
    size_t num_pts = 1;
    for (int II = 0; II < num_dims; II++) {
        start[II] = 0;
        if ((retval = nc_inq_dim(ncid, dim_ids[II] , NULL, &count[II]  ))) { NC_ERR(retval, __LINE__, __FILE__); }
        #if DEBUG >= 2
        fprintf(stdout, "%zu ", count[II]);
        #endif
        num_pts *= count[II];
    }
    #if DEBUG >= 2
    fprintf(stdout, "\n");
    #endif

    // Now resize the vector to the appropriate size
    var.resize(num_pts);

    // Now read in the data
    if ((retval = nc_get_vara_double(ncid, var_id, start, count, &var[0]))) { NC_ERR(retval, __LINE__, __FILE__); }

    // Apply scale factor if appropriate
    double scale = 1.;
    if ((retval = nc_get_att_double(ncid, var_id, "scale_factor", &scale))) { NC_ERR(retval, __LINE__, __FILE__); }
    if (scale != 1.) {
        #if DEBUG >= 2
        fprintf(stdout, "  scale factor = %g\n", scale);
        #endif
        for (size_t II = 0; II < num_pts; II++) {
            var.at(II) = var.at(II) * scale;
        }
    }

    // Apply offset if appropriate
    double offset = 0.;
    if ((retval = nc_get_att_double(ncid, var_id, "add_offset", &offset))) { NC_ERR(retval, __LINE__, __FILE__); }
    if (offset != 0.) {
        #if DEBUG >= 2
        fprintf(stdout, "  additive offset = %g\n", offset);
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
        if ((retval = nc_get_att_double(ncid, var_id, "_FillValue", &fill_val))) { NC_ERR(retval, __LINE__, __FILE__); }
        #if DEBUG >= 2
        fprintf(stdout, "  fill value = %g\n", fill_val);
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
        fprintf(stdout, "  Land cover = %.4g%%\n", 100 * ((double)num_land) / (num_land + num_water));
        #endif
    }
}
