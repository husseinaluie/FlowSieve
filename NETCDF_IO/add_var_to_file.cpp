
#include <vector>
#include <string>
#include "../netcdf_io.hpp"
#include "../constants.hpp"
#include <cassert>

void add_var_to_file(
        const std::string var_name,
        const char ** dim_list,
        const int num_dims,
        const char * filename
        ) {

    static_assert( 
                   not(     (constants::CAST_TO_SINGLE) 
                        and (constants::CAST_TO_INT) 
                       ),
                   "Cannot cast to both single and int. Update cast flags in constants.hpp\n"
                 );

    int datatype;
    if      (constants::CAST_TO_SINGLE) { datatype = NC_FLOAT;  }
    else if (constants::CAST_TO_INT   ) { datatype = NC_SHORT;  }
    else                                { datatype = NC_DOUBLE; }

    // Open the NETCDF file
    int FLAG = NC_WRITE;
    int ncid=0, retval;
    char buffer [50];
    snprintf(buffer, 50, filename);
    retval = nc_open(buffer, FLAG, &ncid);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

    // Extract dimension ids sizes
    int dim_ids[num_dims];
    for (int dim_ind = 0; dim_ind < num_dims; dim_ind++) {
        retval = nc_inq_dimid(ncid, dim_list[dim_ind], &dim_ids[dim_ind] );
        if (retval) { NC_ERR(retval, __LINE__, __FILE__); }
    }

    // Declare the variable
    int var_id;
    char varname [50];
    snprintf(varname, 50, var_name.c_str());
    retval = nc_def_var(ncid, varname, datatype, num_dims, dim_ids, &var_id);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

    // Add the fill value
    const double fill_value = constants::fill_value;
    const signed short fill_value_s = constants::fill_value_s;
    if (constants::CAST_TO_INT) {
        retval = nc_put_att_short( ncid, var_id, "_FillValue", datatype, 1, &fill_value_s);
    } else {
        retval = nc_put_att_double(ncid, var_id, "_FillValue", datatype, 1, &fill_value);
    }
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

    // Close the file
    retval = nc_close(ncid);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

    #if DEBUG >= 2
    fprintf(stdout, "  - added %s to %s -\n", varname, filename);
    #endif
}
