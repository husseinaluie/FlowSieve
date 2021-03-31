#include <math.h>
#include <vector>
#include <string>
#include <mpi.h>
#include "../netcdf_io.hpp"
#include "../constants.hpp"
#include "../postprocess.hpp"

void initialize_postprocess_file(
        const dataset & source_data,
        const std::vector<double> & OkuboWeiss_dim_vals,
        const std::vector<std::string> & int_vars,
        const char * filename,
        const double & filter_scale,
        const bool include_OkuboWeiss,
        const MPI_Comm comm
        ) {

    // Create some tidy names for variables
    const std::vector<double>   &time       = source_data.time,
                                &depth      = source_data.depth,
                                &latitude   = source_data.latitude,
                                &longitude  = source_data.longitude;

    // Get some MPI info
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
    retval = nc_put_att_text(ncid, NC_GLOBAL, "coord-type", 10, constants::CARTESIAN ? "cartesian" : "spherical");
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

    // Extract dimension sizes
    const int   Ntime   = time.size(),
                Ndepth  = depth.size(),
                Nlat    = latitude.size(),
                Nlon    = longitude.size(),
                Nregion = source_data.region_names.size(),
                Nokubo  = OkuboWeiss_dim_vals.size();

    // Define the dimensions
    int time_dimid, depth_dimid, lat_dimid, lon_dimid, reg_dimid, Okubo_dimid;
    retval = nc_def_dim(ncid, "time",      Ntime,     &time_dimid);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }
    retval = nc_def_dim(ncid, "depth",     Ndepth,    &depth_dimid);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }
    retval = nc_def_dim(ncid, "latitude",  Nlat,      &lat_dimid);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }
    retval = nc_def_dim(ncid, "longitude", Nlon,      &lon_dimid);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }
    retval = nc_def_dim(ncid, "region",    Nregion,   &reg_dimid);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }
    if (include_OkuboWeiss) {
        retval = nc_def_dim(ncid, "OkuboWeiss", Nokubo, &Okubo_dimid);
        if (retval) { NC_ERR(retval, __LINE__, __FILE__); }
    }

    // Define coordinate variables
    int time_varid, depth_varid, lat_varid, lon_varid, reg_varid, Okubo_varid;
    retval = nc_def_var(ncid, "time",      NC_DOUBLE,  1, &time_dimid,  &time_varid);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }
    retval = nc_def_var(ncid, "depth",     NC_DOUBLE,  1, &depth_dimid, &depth_varid);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }
    retval = nc_def_var(ncid, "latitude",  NC_DOUBLE,  1, &lat_dimid,   &lat_varid);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }
    retval = nc_def_var(ncid, "longitude", NC_DOUBLE,  1, &lon_dimid,   &lon_varid);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }
    retval = nc_def_var(ncid, "region",    NC_STRING, 1,  &reg_dimid,   &reg_varid);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }
    if (include_OkuboWeiss) {
        retval = nc_def_var(ncid, "OkuboWeiss",    NC_DOUBLE, 1,  &Okubo_dimid,   &Okubo_varid);
        if (retval) { NC_ERR(retval, __LINE__, __FILE__); }
    }

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

    if (include_OkuboWeiss) {
        count[0] = Nokubo;
        retval = nc_put_vara_double(ncid, Okubo_varid, start, count, &OkuboWeiss_dim_vals[0]);
        if (retval) { NC_ERR(retval, __LINE__, __FILE__); }
    }

    // We're also going to store the region areas
    int area_dims[3];
    area_dims[0] = time_dimid;
    area_dims[1] = depth_dimid;
    area_dims[2] = reg_dimid;
    int reg_area_varid;
    retval = nc_def_var(ncid, "region_areas", NC_DOUBLE, 3, area_dims, &reg_area_varid);
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

        // region averages
        const char* dim_names[] = {"time", "depth", "region"};
        const int ndims = 3;
        for (size_t varInd = 0; varInd < int_vars.size(); ++varInd) {
            add_var_to_file(int_vars.at(varInd)+"_avg", dim_names, ndims, buffer);
            //add_var_to_file(int_vars.at(varInd)+"_std", dim_names, ndims, buffer);
        }

        // time averages
        const char* dim_names_time_ave[] = {"depth", "latitude", "longitude"};
        const int ndims_time_ave = 3;
        for (size_t varInd = 0; varInd < int_vars.size(); ++varInd) {
            add_var_to_file(int_vars.at(varInd)+"_time_average", dim_names_time_ave, ndims_time_ave, buffer);
            //add_var_to_file(int_vars.at(varInd)+"_time_std_dev", dim_names_time_ave, ndims_time_ave, buffer);
        }

        // region averages : averaged over OkuboWeiss contours
        if (include_OkuboWeiss) {
            const char* dim_names[] = {"time", "depth", "OkuboWeiss", "region"};
            const int ndims = 4;
            add_var_to_file("area_OkuboWeiss", dim_names, ndims, buffer);
            for (size_t varInd = 0; varInd < int_vars.size(); ++varInd) {
                add_var_to_file(int_vars.at(varInd)+"_avg_OkuboWeiss", dim_names, ndims, buffer);
                //add_var_to_file(int_vars.at(varInd)+"_std_OkuboWeiss", dim_names, ndims, buffer);
            }
        }
    }

    // Add some global attributes from constants.hpp
    add_attr_to_file("R_earth",                                      constants::R_earth,    filename);
    add_attr_to_file("rho0",                                         constants::rho0,       filename);
    add_attr_to_file("g",                                            constants::g,          filename);
    add_attr_to_file("differentiation_convergence_order",   (double) constants::DiffOrd,    filename);
    add_attr_to_file("KERNEL_OPT",                          (double) constants::KERNEL_OPT, filename);
    add_attr_to_file("KernPad",                             (double) constants::KernPad,    filename);

    // Write region names - this has to be done separately for reasons
    write_regions_to_post( filename, source_data.region_names );

    // Write region areas
    size_t start_r[] = { source_data.myStarts.at(0), source_data.myStarts.at(1),  0       }, 
           count_r[] = { source_data.Ntime,          source_data.Ndepth,          Nregion };
    write_field_to_output( source_data.region_areas, "region_areas", start_r, count_r, filename, NULL);

    #if DEBUG >= 2
    if (wRank == 0) { fprintf(stdout, "\n"); }
    #endif
}
