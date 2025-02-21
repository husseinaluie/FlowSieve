
#include "../netcdf_io.hpp"
#include "../constants.hpp"
#include <string.h>
#include <cassert>
#include <math.h>
#include <fenv.h>

// Write to netcdf file
void read_attr_from_file(
        double &attr,
        const std::string & attr_name,
        const std::string & filename,
        const std::string & var_name,
        const MPI_Comm comm
        ) {

    assert( check_file_existence( filename ) );

    int wRank, wSize;
    MPI_Comm_rank( comm, &wRank );
    MPI_Comm_size( comm, &wSize );

    #if DEBUG >= 1
    if (wRank == 0) {
        fprintf(stdout, "Attempting to read %s from %s\n", attr_name.c_str(), filename.c_str());
    }
    #endif

    // Open the NETCDF file
    //int FLAG = NC_NETCDF4 | NC_NOWRITE | NC_MPIIO;
    int FLAG = NC_NETCDF4 | NC_MPIIO;
    int ncid=0, retval;

    // Some netcdf functions [in some netcdf versions] cause floating-point errors
    //  so, we need to disable floating point exceptions when we try to open files.
    fedisableexcept( FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW );
    retval = nc_open_par( filename.c_str(), FLAG, comm, MPI_INFO_NULL, &ncid);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

    // Now we can restore fp-exception handling, after clearing out any that were raised
    feclearexcept( FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW ); // erase whatever exceptions were raised
    feenableexcept( FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW ); // re-enable exceptions

    // Get information about the variable
    int var_id = NC_GLOBAL, num_dims;
    int dim_ids[NC_MAX_VAR_DIMS];
    //if (var_name == NULL) {
    if ( var_name.empty() ) {
        // If no var_name given, then assume a global attribute
        var_id = NC_GLOBAL;
    } else {
        // Otherwise, get the appropriate variable id
        retval = nc_inq_var(ncid, var_id, NULL, NULL, &num_dims, dim_ids, NULL );
        if (retval) { NC_ERR(retval, __LINE__, __FILE__); }
    }

    nc_type att_type;
    nc_inq_atttype(ncid, var_id, attr_name.c_str(), &att_type);
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
        case NC_SHORT   : nc_get_att(ncid, var_id, attr_name.c_str(), &att_val_short);   
                          attr = double(att_val_short);   break;
        case NC_USHORT  : nc_get_att(ncid, var_id, attr_name.c_str(), &att_val_ushort);  
                          attr = double(att_val_ushort);  break;
        case NC_INT     : nc_get_att(ncid, var_id, attr_name.c_str(), &att_val_int);     
                          attr = double(att_val_int);     break;
        case NC_UINT    : nc_get_att(ncid, var_id, attr_name.c_str(), &att_val_uint);    
                          attr = double(att_val_uint);    break;
        case NC_INT64   : nc_get_att(ncid, var_id, attr_name.c_str(), &att_val_long);    
                          attr = double(att_val_long);    break;
        case NC_UINT64  : nc_get_att(ncid, var_id, attr_name.c_str(), &att_val_ulong);   
                          attr = double(att_val_ulong);   break;
        case NC_FLOAT   : nc_get_att(ncid, var_id, attr_name.c_str(), &att_val_float);   
                          attr = double(att_val_float);   break;
        case NC_DOUBLE  : nc_get_att(ncid, var_id, attr_name.c_str(), &att_val_double);  
                          attr = double(att_val_double);  break;
    }

    #if DEBUG >= 1
    if (wRank == 0) {
        fprintf(stdout, "    %s = %g\n", attr_name.c_str(), attr);
    }
    #endif

    MPI_Barrier(comm);
    retval = nc_close(ncid);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }
}
