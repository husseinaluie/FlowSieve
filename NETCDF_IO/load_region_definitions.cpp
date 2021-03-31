
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

    int wRank, wSize, Nprocs_in_dim, Iproc_in_dim;
    MPI_Comm_rank( comm, &wRank );
    MPI_Comm_size( comm, &wSize );

    assert( ( Nlat > 0 ) and ( Nlon > 0 ) ); // The lat and lon grids should already exist
    assert( Nlat * Nlon < 4 * pow(512,3) ); // Grid is too large. IO scripts will need to be modified to handle large fields.

    // Open the NETCDF file
    const int str_len = 100;
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
    }

    // While we're here, might as well also compute the area of each region
    assert( mask.size() > 0 ); // must read in mask before computing region areas
    assert( areas.size() > 0 ); // must compute cell areas before computing region areas
    region_areas.resize( num_regions * Ntime * Ndepth );

    double local_area;
    size_t Ilat, Ilon, reg_index, index, area_index;

    const int chunk_size = get_omp_chunksize(Nlat, Nlon);

    for (size_t Iregion = 0; Iregion < num_regions; ++Iregion) {
        fprintf( stdout, " Read in region %s with %'zu points\n", region_names.at(Iregion).c_str(), regions.at( region_names.at(Iregion) ).size() );
        for (size_t Itime = 0; Itime < Ntime; ++Itime) {
            for (size_t Idepth = 0; Idepth < Ndepth; ++Idepth) {

                local_area = 0.;

                #pragma omp parallel default(none)\
                private( Ilat, Ilon, index, area_index )\
                shared( mask, areas, Iregion, Itime, Idepth ) \
                reduction(+ : local_area)
                { 
                    #pragma omp for collapse(2) schedule(guided, chunk_size)
                    for (Ilat = 0; Ilat < Nlat; ++Ilat) {
                        for (Ilon = 0; Ilon < Nlon; ++Ilon) {

                            index = Index(Itime, Idepth, Ilat, Ilon,
                                          Ntime, Ndepth, Nlat, Nlon);

                            if ( mask.at(index) ) { // Skip land areas

                                area_index = Index(0, 0, Ilat, Ilon,
                                                   1, 1, Nlat, Nlon);

                                if ( regions.at( region_names.at( Iregion ) ).at(area_index) ) {
                                    local_area += areas.at(area_index); 
                                }
                            }
                        }
                    }
                }
                reg_index = Index(0, Itime, Idepth, Iregion,
                                  1, Ntime, Ndepth, num_regions);

                region_areas.at(reg_index) = local_area;
                
                fprintf(stdout, "   area = %g (%g)\n", region_areas.at(reg_index), local_area);
            }
        }
    }

}
