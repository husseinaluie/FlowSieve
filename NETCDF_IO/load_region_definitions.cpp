
#include "../netcdf_io.hpp"
#include "../constants.hpp"
#include "../functions.hpp"
#include <cassert>

void dataset::load_region_definitions(
        const std::string filename,
        const std::string dim_name,
        const std::string var_name,
        const MPI_Comm comm
        ) {

    assert( check_file_existence( filename.c_str() ) ); // Make sure the file exists

    int wRank, wSize;
    MPI_Comm_rank( comm, &wRank );
    MPI_Comm_size( comm, &wSize );

    assert( ( Nlat > 0 ) and ( Nlon > 0 ) ); // The lat and lon grids should already exist
    assert( Nlat * Nlon < 4 * pow(512,3) ); // Grid is too large. IO scripts will need to be modified to handle large fields.

    // Open the NETCDF file
    const int str_len = 250;
    int FLAG = NC_NETCDF4 | NC_MPIIO;
    int ncid=0, retval;
    char buffer [str_len];
    snprintf(buffer, str_len, filename.c_str());

    #if DEBUG >= 1
    if (wRank == 0) {
        fprintf(stdout, "Attempting to read region information from %s\n", buffer);
        fflush(stdout);
    }
    #endif

    retval = nc_open_par(buffer, FLAG, comm, MPI_INFO_NULL, &ncid);
    if (retval != NC_NOERR ) { NC_ERR(retval, __LINE__, __FILE__); }

    int dim_id, name_id;
    retval = nc_inq_dimid(ncid, dim_name.c_str(), &dim_id );
    if (retval != NC_NOERR ) { NC_ERR(retval, __LINE__, __FILE__); }

    retval = nc_inq_varid(ncid, dim_name.c_str(), &name_id );
    if (retval != NC_NOERR ) { NC_ERR(retval, __LINE__, __FILE__); }

    size_t num_regions;
    retval = nc_inq_dim(ncid, dim_id, NULL, &num_regions);
    if (retval != NC_NOERR ) { NC_ERR(retval, __LINE__, __FILE__); }
    region_names.resize(num_regions);
    #if DEBUG >= 1
    if (wRank == 0) {
        fprintf(stdout, "  %'zu regions to be read in\n", num_regions);
    }
    #endif

    int var_id;
    retval = nc_inq_varid(ncid, var_name.c_str(), &var_id );
    if (retval != NC_NOERR ) { NC_ERR(retval, __LINE__, __FILE__); }


    // Set up to read in each region one at a time.
    size_t start[3], count[3];
    start[1] = 0;
    start[2] = 0;
    count[0] = 1;
    count[1] = Nlat;
    count[2] = Nlon;

    const size_t name_count = 1;
    
    //  Buffer to store region info as it's read in
    //    These won't have been stored as booleans, but we will
    //    convert them now to save on memory
    std::vector<double> var_buffer( Nlat * Nlon );
    std::string region_name;
    region_name.resize(str_len);

    char* name_buffer = &region_name[0];
    
    for (size_t Iregion = 0; Iregion < num_regions; ++Iregion) {

        start[0] = Iregion;

        // Get the name of the region
        retval = nc_get_vara_string( ncid, name_id, &Iregion, &name_count, &name_buffer );
        region_name = name_buffer;

        region_names.at(Iregion) = region_name;

        retval = nc_get_vara_double(ncid, var_id, start, count, &var_buffer[0]);
        if (retval != NC_NOERR ) { NC_ERR(retval, __LINE__, __FILE__); }

        regions.insert( std::pair< std::string, std::vector<bool> >( region_name, std::vector<bool>( var_buffer.size() ) ) );
        for (size_t II = 0; II < var_buffer.size(); ++II) {
            regions.at( region_name ).at( II ) = var_buffer.at( II ) == 1 ? true : false;
        }

        #if DEBUG >= 1
        if (wRank == 0) {
            fprintf( stdout, " Read in region %s\n", region_name.c_str() );
        }
        #endif
    }

    // While we're here, might as well also compute the area of each region
    compute_region_areas();

}
