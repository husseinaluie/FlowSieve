
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

    #if DEBUG >= 2
    fprintf(stdout, "  - added %s to %s -\n", varname, filename);
    #endif

    ////
    //    Now add various meta-data statements to the file for clarity (if available)
    ////
    
    // First, check 'type' (i.e. time average, spatial average, toroidal potential, etc)
    std::string substring = var_name;

    const size_t    time_avg_index  = var_name.find("_time_average"),
                    Okubo_avg_index = var_name.find("_OkuboWeiss_average"),
                    area_avg_index  = var_name.find("_area_average"),
                    zonal_avg_index = var_name.find("_zonal_average");
    if ( time_avg_index != std::string::npos ) {
        nc_put_att_text( ncid, var_id, "variable_type", 
                         constants::time_average_description.size(), 
                         constants::time_average_description.c_str() );
        substring = var_name.substr( 0, time_avg_index );

    } else if ( Okubo_avg_index != std::string::npos ) {
        nc_put_att_text( ncid, var_id, "variable_type", 
                         constants::OkuboWeiss_average_description.size(), 
                         constants::OkuboWeiss_average_description.c_str() );
        substring = var_name.substr( 0, Okubo_avg_index );

    } else if ( area_avg_index != std::string::npos ) {
        nc_put_att_text( ncid, var_id, "variable_type", 
                         constants::spatial_average_description.size(),    
                         constants::spatial_average_description.c_str() );
        substring = var_name.substr( 0, area_avg_index );

    } else if ( zonal_avg_index != std::string::npos ) {
        nc_put_att_text( ncid, var_id, "variable_type", 
                         constants::zonal_average_description.size(),    
                         constants::zonal_average_description.c_str() );
        substring = var_name.substr( 0, zonal_avg_index );

    }

    // Next, check if we have a long-form description and add that
    if ( constants::variable_descriptions.count( substring ) ) {
        nc_put_att_text( ncid, var_id, "long_name", 
                         constants::variable_descriptions.at( substring ).size(),
                         constants::variable_descriptions.at( substring ).c_str() );
    }

    // Next, check if we have defined units, and add them
    if ( constants::variable_units.count( substring ) ) {
        nc_put_att_text( ncid, var_id, "units",     
                         constants::variable_units.at( substring ).size(),
                         constants::variable_units.at( substring ).c_str() );
    }

    ////
    //    Close the file
    ////
    retval = nc_close(ncid);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

}
