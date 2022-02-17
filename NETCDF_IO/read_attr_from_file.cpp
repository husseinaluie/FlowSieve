
#include "../netcdf_io.hpp"
#include "../constants.hpp"
#include <string.h>
#include <cassert>
#include <math.h>

// Write to netcdf file
void read_attr_from_file(
        double &attr,
        const char * attr_name,
        const std::string filename,
        const char * var_name,
        const MPI_Comm comm
        ) {

    assert( check_file_existence( filename ) );

    int wRank, wSize;
    MPI_Comm_rank( comm, &wRank );
    MPI_Comm_size( comm, &wSize );

    #if DEBUG >= 1
    if (wRank == 0) {
        fprintf(stdout, "Attempting to read %s from %s\n", attr_name, filename.c_str());
    }
    #endif

    // Open the NETCDF file
    //int FLAG = NC_NETCDF4 | NC_NOWRITE | NC_MPIIO;
    int FLAG = NC_NETCDF4 | NC_MPIIO;
    int ncid=0, retval;
    char buffer [50];
    snprintf(buffer, 50, filename.c_str());
    retval = nc_open_par(buffer, FLAG, comm, MPI_INFO_NULL, &ncid);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

    // Get information about the variable
    int var_id = NC_GLOBAL, num_dims;
    int dim_ids[NC_MAX_VAR_DIMS];
    if (var_name == NULL) {
        // If no var_name given, then assume a global attribute
        var_id = NC_GLOBAL;
    } else {
        // Otherwise, get the appropriate variable id
        retval = nc_inq_var(ncid, var_id, NULL, NULL, &num_dims, dim_ids, NULL );
        if (retval) { NC_ERR(retval, __LINE__, __FILE__); }
    }

    char attrname [50];
    snprintf(attrname, 50, attr_name);

    nc_type att_type;
    nc_inq_atttype(ncid, var_id, attrname, &att_type);
    #if DEBUG >= 1
    if (wRank == 0) {
        switch (att_type) {
            case NC_SHORT   : fprintf( stdout, "Attribute type is NC_SHORT\n"     ); break;
            case NC_USHORT  : fprintf( stdout, "Attribute type is NC_USHORT\n"    ); break;
            case NC_INT     : fprintf( stdout, "Attribute type is NC_INT\n"       ); break;
            case NC_UINT    : fprintf( stdout, "Attribute type is NC_UINT\n"      ); break;
            case NC_INT64   : fprintf( stdout, "Attribute type is NC_INT64\n"     ); break;
            case NC_UINT64  : fprintf( stdout, "Attribute type is NC_UINT64\n"    ); break;
            case NC_FLOAT   : fprintf( stdout, "Attribute type is NC_FLOAT\n"     ); break;
            case NC_DOUBLE  : fprintf( stdout, "Attribute type is NC_DOUBLE\n"    ); break;
            case NC_STRING  : fprintf( stdout, "Attribute type is NC_STRING\n"    ); break;
            case NC_CHAR    : fprintf( stdout, "Attribute type is NC_CHAR\n"      ); break;
            case NC_BYTE    : fprintf( stdout, "Attribute type is NC_BYTE\n"      ); break;
            case NC_UBYTE   : fprintf( stdout, "Attribute type is NC_UBYTE\n"     ); break;
            default         : fprintf( stdout, "Attribute type not recognized.\n" ); break;
        }
    }
    #endif
    short att_val_short;
    unsigned short att_val_ushort;
    int att_val_int;
    unsigned int att_val_uint;
    long att_val_long;
    unsigned long att_val_ulong;
    float att_val_float;
    double att_val_double;
    switch (att_type) {
        case NC_SHORT   : nc_get_att(ncid, var_id, attrname, &att_val_short);   attr = double(att_val_short);   break;
        case NC_USHORT  : nc_get_att(ncid, var_id, attrname, &att_val_ushort);  attr = double(att_val_ushort);  break;
        case NC_INT     : nc_get_att(ncid, var_id, attrname, &att_val_int);     attr = double(att_val_int);     break;
        case NC_UINT    : nc_get_att(ncid, var_id, attrname, &att_val_uint);    attr = double(att_val_uint);    break;
        case NC_INT64   : nc_get_att(ncid, var_id, attrname, &att_val_long);    attr = double(att_val_long);    break;
        case NC_UINT64  : nc_get_att(ncid, var_id, attrname, &att_val_ulong);   attr = double(att_val_ulong);   break;
        case NC_FLOAT   : nc_get_att(ncid, var_id, attrname, &att_val_float);   attr = double(att_val_float);   break;
        case NC_DOUBLE  : nc_get_att(ncid, var_id, attrname, &att_val_double);  attr = double(att_val_double);  break;
    }

    //nc_get_att(ncid, var_id, attrname, &attr);

    #if DEBUG >= 1
    if (wRank == 0) {
        fprintf(stdout, "    %s = %g\n", attr_name, attr);
    }
    #endif

    MPI_Barrier(comm);
    retval = nc_close(ncid);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }
}
