#include <math.h>
#include <vector>
#include <mpi.h>
#include <cassert>
#include "../netcdf_io.hpp"
#include "../constants.hpp"

void initialize_adjacency_file(
        const dataset & source_data,
        const std::vector<std::string> & vars,
        const char * filename,
        const double filter_scale,
        const MPI_Comm comm
        ) {

    int wRank=-1, wSize=-1;
    MPI_Comm_rank( MPI_COMM_WORLD, &wRank );
    MPI_Comm_size( MPI_COMM_WORLD, &wSize );

    assert( constants::GRID_TYPE == constants::GridType::LLC ); // No adjacency to write otherwise

    const int NeighbourPopulation = source_data.num_neighbours;

    #if DEBUG>=1
    if (wRank == 0) { fprintf(stdout, "\nPreparing to initialize the output file.\n"); }
    #endif

    // Create some tidy names for variables
    const std::vector<double>   &latitude   = source_data.latitude,
                                &longitude  = source_data.longitude,
                                &areas      = source_data.areas;

    // Open the NETCDF file
    int FLAG = NC_NETCDF4 | NC_CLOBBER | NC_MPIIO;
    int ncid=0, retval;
    char buffer [50];
    snprintf(buffer, 50, filename);
    retval = nc_create_par(buffer, FLAG, comm, MPI_INFO_NULL, &ncid);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

    #if DEBUG>=2
    if (wRank == 0) { fprintf(stdout, "    Logging the filter scale\n"); }
    #endif
    if ( filter_scale >= 0 ) {
        retval = nc_put_att_double(ncid, NC_GLOBAL, "filter_scale", NC_DOUBLE, 1, &filter_scale);
        if (retval) { NC_ERR(retval, __LINE__, __FILE__); }
    }

    retval = nc_put_att_double(ncid, NC_GLOBAL, "rho0", NC_DOUBLE, 1, &constants::rho0);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

    // Record coordinate type
    #if DEBUG>=2
    if (wRank == 0) { fprintf(stdout, "    Logging the grid type\n"); }
    #endif
    if (constants::CARTESIAN) {
        retval = nc_put_att_text(ncid, NC_GLOBAL, "coord-type", 10, "cartesian");
    } else {
        retval = nc_put_att_text(ncid, NC_GLOBAL, "coord-type", 10, "spherical");
    }
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

    // Extract dimension sizes
    const int   Nlat    = latitude.size(),
                Nlon    = longitude.size();

    // Define the dimensions
    #if DEBUG>=2
    if (wRank == 0) { fprintf(stdout, "    Defining the dimensions\n"); }
    #endif
    int lat_dimid, lon_dimid, neighbour_dimid;
    
    retval = nc_def_dim(ncid, "latlon",  Nlat,      &lat_dimid);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }
    lon_dimid = lat_dimid;

    retval = nc_def_dim(ncid, "neighbour",  NeighbourPopulation+1, &neighbour_dimid);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

    // Define coordinate variables
    #if DEBUG>=2
    if (wRank == 0) { fprintf(stdout, "    Defining the dimension variables\n"); }
    #endif
    int lat_varid, lon_varid, neighbour_varid;
    retval = nc_def_var(ncid, "latitude",  NC_DOUBLE, 1, &lat_dimid,   &lat_varid);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }
    retval = nc_def_var(ncid, "longitude", NC_DOUBLE, 1, &lon_dimid,   &lon_varid);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }
    //retval = nc_def_var(ncid, "neighbour", NC_DOUBLE, 1, &neighbour_dimid,   &neighbour_varid);
    //if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

    if (not(constants::CARTESIAN)) {
        #if DEBUG>=2
        if (wRank == 0) { fprintf(stdout, "    Add scale factors for Rad to Degrees\n"); }
        #endif
        const double rad_to_degree = 180. / M_PI;
        retval = nc_put_att_double(ncid, lon_varid, "scale_factor", 
                NC_DOUBLE, 1, &rad_to_degree);
        if (retval) { NC_ERR(retval, __LINE__, __FILE__); }
        retval = nc_put_att_double(ncid, lat_varid, "scale_factor", 
                NC_DOUBLE, 1, &rad_to_degree);
        if (retval) { NC_ERR(retval, __LINE__, __FILE__); }
    }

    // Write the coordinate variables
    #if DEBUG>=2
    if (wRank == 0) { fprintf(stdout, "    Write the dimensions\n"); }
    #endif
    size_t start[1], count[1];
    start[0] = 0;
    count[0] = Nlat;
    retval = nc_put_vara_double(ncid, lat_varid,   start, count, &latitude[0]);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

    count[0] = Nlon;
    retval = nc_put_vara_double(ncid, lon_varid,   start, count, &longitude[0]);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

    // Close the file
    retval = nc_close(ncid);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

    #if DEBUG >= 2
    if (wRank == 0) { fprintf(stdout, "\nAdjacency file (%s) initialized.\n", buffer); }
    #endif

    if (wRank == 0) {
        #if DEBUG>=2
        if (wRank == 0) { fprintf(stdout, "    Root rank will now add each variable.\n"); }
        #endif
        // Loop through and add the desired variables
        // Dimension names (in order!)
        const char* dim_names_LLC[]      = {"latlon", "neighbour"};
        const int ndims_LLC = 2;
        for (size_t varInd = 0; varInd < vars.size(); ++varInd) {
            add_var_to_file( vars.at(varInd), dim_names_LLC, ndims_LLC, buffer );
        }
    }

    // Add some global attributes from constants.hpp
    #if DEBUG>=2
    if (wRank == 0) { fprintf(stdout, "    Now add a series of computational details\n"); }
    #endif
    add_attr_to_file("R_earth",                                      constants::R_earth,    filename);

    #if DEBUG >= 2
    if (wRank == 0) { fprintf(stdout, "\n"); }
    #endif
}
