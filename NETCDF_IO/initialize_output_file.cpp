#include <math.h>
#include <vector>
#include <mpi.h>
#include "../netcdf_io.hpp"
#include "../constants.hpp"

void initialize_output_file(
        const std::vector<double> & time,
        const std::vector<double> & depth,
        const std::vector<double> & longitude,
        const std::vector<double> & latitude,
        const std::vector<double> & areas,
        const std::vector<std::string> & vars,
        const char * filename,
        const double filter_scale,
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

    retval = nc_put_att_double(ncid, NC_GLOBAL, "filter_scale", NC_DOUBLE, 1, &filter_scale);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

    retval = nc_put_att_double(ncid, NC_GLOBAL, "rho0", NC_DOUBLE, 1, &constants::rho0);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

    // Record coordinate type
    if (constants::CARTESIAN) {
        retval = nc_put_att_text(ncid, NC_GLOBAL, "coord-type", 10, "cartesian");
    } else {
        retval = nc_put_att_text(ncid, NC_GLOBAL, "coord-type", 10, "spherical");
    }
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

    // Extract dimension sizes
    const int Ntime   = time.size();
    const int Ndepth  = depth.size();
    const int Nlat    = latitude.size();
    const int Nlon    = longitude.size();

    // Define the dimensions
    int time_dimid, depth_dimid, lat_dimid, lon_dimid;
    retval = nc_def_dim(ncid, "time",      Ntime,     &time_dimid);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }
    retval = nc_def_dim(ncid, "depth",     Ndepth,    &depth_dimid);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }
    retval = nc_def_dim(ncid, "latitude",  Nlat,      &lat_dimid);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }
    retval = nc_def_dim(ncid, "longitude", Nlon,      &lon_dimid);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

    // Define coordinate variables
    int time_varid, depth_varid, lat_varid, lon_varid;
    retval = nc_def_var(ncid, "time",      NC_DOUBLE, 1, &time_dimid,  &time_varid);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }
    retval = nc_def_var(ncid, "depth",     NC_DOUBLE, 1, &depth_dimid, &depth_varid);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }
    retval = nc_def_var(ncid, "latitude",  NC_DOUBLE, 1, &lat_dimid,   &lat_varid);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }
    retval = nc_def_var(ncid, "longitude", NC_DOUBLE, 1, &lon_dimid,   &lon_varid);
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
    count[0] = Ntime;
    retval = nc_put_vara_double(ncid, time_varid,  start, count, &time[0]);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

    count[0] = Ndepth;
    retval = nc_put_vara_double(ncid, depth_varid, start, count, &depth[0]);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

    count[0] = Nlat;
    retval = nc_put_vara_double(ncid, lat_varid,   start, count, &latitude[0]);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

    count[0] = Nlon;
    retval = nc_put_vara_double(ncid, lon_varid,   start, count, &longitude[0]);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

    // Write the cell areas for convenience
    int area_dimids[2];
    area_dimids[0] = lat_dimid;
    area_dimids[1] = lon_dimid;
    int area_varid;
    retval = nc_def_var(ncid, "cell_areas", NC_DOUBLE, 2, area_dimids, &area_varid);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

    size_t area_start[2], area_count[2];
    area_start[0] = 0;
    area_start[1] = 0;
    area_count[0] = Nlat;
    area_count[1] = Nlon;
    retval = nc_put_vara_double(ncid, area_varid, area_start, area_count, &areas[0]);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

    // Close the file
    retval = nc_close(ncid);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

    #if DEBUG >= 2
    if (wRank == 0) { fprintf(stdout, "\nOutput file (%s) initialized.\n", buffer); }
    #endif

    if (wRank == 0) {
        // Loop through and add the desired variables
        // Dimension names (in order!)
        const char* dim_names[] = {"time", "depth", "latitude", "longitude"};
        const int ndims = 4;
        for (size_t varInd = 0; varInd < vars.size(); ++varInd) {
            add_var_to_file(vars.at(varInd), dim_names, ndims, buffer);
        }
    }

    // Add some global attributes from constants.hpp
    add_attr_to_file("R_earth",                                      constants::R_earth,    filename);
    add_attr_to_file("rho0",                                         constants::rho0,       filename);
    add_attr_to_file("g",                                            constants::g,          filename);
    add_attr_to_file("differentiation_convergence_order",   (double) constants::DiffOrd,    filename);
    add_attr_to_file("KERNEL_OPT",                          (double) constants::KERNEL_OPT, filename);
    add_attr_to_file("KernPad",                             (double) constants::KernPad,    filename);

    #if DEBUG >= 2
    if (wRank == 0) { fprintf(stdout, "\n"); }
    #endif
}
