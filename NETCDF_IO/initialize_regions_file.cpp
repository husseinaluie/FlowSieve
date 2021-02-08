#include <math.h>
#include <vector>
#include <string>
#include <mpi.h>
#include "../netcdf_io.hpp"
#include "../constants.hpp"
#include "../postprocess.hpp"

void initialize_regions_file(
        const std::vector<double> & latitude,
        const std::vector<double> & longitude,
        const char * filename,
        const MPI_Comm comm
        ) {

    int wRank=-1, wSize=-1;
    MPI_Comm_rank( MPI_COMM_WORLD, &wRank );
    MPI_Comm_size( MPI_COMM_WORLD, &wSize );

    // Open the NETCDF file
    int FLAG = NC_NETCDF4 | NC_CLOBBER | NC_MPIIO;
    int ncid=0, retval;
    char buffer [50];
    snprintf(buffer, 50, filename);
    retval = nc_create_par(buffer, FLAG, comm, MPI_INFO_NULL, &ncid);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

    // Record coordinate type
    retval = nc_put_att_text(ncid, NC_GLOBAL, "coord-type", 10, constants::CARTESIAN ? "cartesian" : "spherical");
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

    // Extract dimension sizes
    const int Nlat    = latitude.size();
    const int Nlon    = longitude.size();
    const int Nregion = RegionTest::all_regions.size();

    // Define the dimensions
    int lat_dimid, lon_dimid, reg_dimid;
    retval = nc_def_dim(ncid, "latitude",  Nlat,      &lat_dimid);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }
    retval = nc_def_dim(ncid, "longitude", Nlon,      &lon_dimid);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }
    retval = nc_def_dim(ncid, "region",    Nregion,   &reg_dimid);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

    // Define coordinate variables
    int lat_varid, lon_varid, reg_varid;
    retval = nc_def_var(ncid, "latitude",  NC_DOUBLE,  1, &lat_dimid,   &lat_varid);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }
    retval = nc_def_var(ncid, "longitude", NC_DOUBLE,  1, &lon_dimid,   &lon_varid);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }
    retval = nc_def_var(ncid, "region",    NC_STRING, 1, &reg_dimid,   &reg_varid);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

    if (not(constants::CARTESIAN)) {
        const double rad_to_degree = 180. / M_PI;
        retval = nc_put_att_double(ncid, lon_varid, "scale_factor", 
                NC_DOUBLE, 1, &rad_to_degree);
        if (retval) { NC_ERR(retval, __LINE__, __FILE__); }
        retval = nc_put_att_double(ncid, lat_varid, "scale_factor", 
                NC_DOUBLE, 1, &rad_to_degree);
        if (retval) { NC_ERR(retval, __LINE__, __FILE__); }
    }

    // Write the coordinate variables
    size_t start[1], count[1];
    start[0] = 0;

    count[0] = Nlat;
    retval = nc_put_vara_double(ncid, lat_varid,   start, count, &latitude[0]);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

    count[0] = Nlon;
    retval = nc_put_vara_double(ncid, lon_varid,   start, count, &longitude[0]);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

    // We're also going to store the region areas
    int reg_area_varid;
    retval = nc_def_var(ncid, "region_areas", NC_DOUBLE, 1, &reg_dimid, &reg_area_varid);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

    // Close the file
    retval = nc_close(ncid);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

    #if DEBUG >= 2
    if (wRank == 0) { fprintf(stdout, "\nOutput file (%s) initialized.\n", buffer); }
    #endif

    if (wRank == 0) {
        // time averages
        const char* dim_names[] = {"region", "latitude", "longitude"};
        const int ndims = 3;
        add_var_to_file("region_definitions", dim_names, ndims, buffer);
    }

    #if DEBUG >= 2
    if (wRank == 0) { fprintf(stdout, "\n"); }
    #endif
}
