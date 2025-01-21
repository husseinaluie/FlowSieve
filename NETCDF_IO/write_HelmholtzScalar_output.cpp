#include <fenv.h>
#include <vector>
#include <string>
#include <mpi.h>
#include <math.h>
#include <cassert>
#include <fenv.h>
#include "../netcdf_io.hpp"
#include "../constants.hpp"
#include "../preprocess.hpp"
#include "../functions.hpp"


void write_HelmholtzScalar_output(
        const std::string & filename,
        const HelmholtzDataClass & Helmholtz_data,
        const dataset & source_data,
        const dataset & solution_data,
        const MPI_Comm comm
        ){


    // Psi, Phi is in solution data
    // u, v, vort, div are in source_data
    // convergence metrics are in Helmholtz_data


    int wRank=-1, wSize=-1;
    MPI_Comm_rank( MPI_COMM_WORLD, &wRank );
    MPI_Comm_size( MPI_COMM_WORLD, &wSize );

    #if DEBUG>=1
    if (wRank == 0) { fprintf(stdout, "\nPreparing to initialize the output file.\n"); }
    #endif

    // Create some tidy names for variables
    const std::vector<double>   &time       = source_data.time,
                                &depth      = source_data.depth,
                                &latitude   = source_data.latitude,
                                &longitude  = source_data.longitude,
                                &areas      = source_data.areas;

    // Extract dimension sizes
    const int Ntime   = time.size(),
              Ndepth  = depth.size(),
              Nlat    = latitude.size(),
              Nlon    = longitude.size();
    const std::vector<int>  &myStarts = source_data.myStarts,
                            &myCounts = source_data.myCounts;

    const int GridType = constants::GRID_TYPE;

    // Open the NETCDF file
    int FLAG = NC_NETCDF4 | NC_CLOBBER | NC_MPIIO;
    int ncid=0, retval;
    retval = nc_create_par( filename.c_str(), FLAG, comm, MPI_INFO_NULL, &ncid);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

    // Define the dimensions
    #if DEBUG>=2
    if (wRank == 0) { fprintf(stdout, "    Defining the dimensions\n"); }
    #endif
    int time_dimid, depth_dimid, lat_dimid, lon_dimid;
    retval = nc_def_dim(ncid, "time",      Ntime,     &time_dimid);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }
    retval = nc_def_dim(ncid, "depth",     Ndepth,    &depth_dimid);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }
    
    if ( GridType == constants::GridType::MeshGrid ) {
        retval = nc_def_dim(ncid, "latitude",  Nlat,      &lat_dimid);
        if (retval) { NC_ERR(retval, __LINE__, __FILE__); }
        retval = nc_def_dim(ncid, "longitude", Nlon,      &lon_dimid);
        if (retval) { NC_ERR(retval, __LINE__, __FILE__); }
    } else if ( GridType == constants::GridType::LLC ) {
        retval = nc_def_dim(ncid, "latlon",  Nlat,      &lat_dimid);
        if (retval) { NC_ERR(retval, __LINE__, __FILE__); }
        lon_dimid = lat_dimid;
    }

    int IterCycle_dimid;
    retval = nc_def_dim(ncid, "IterCycle", Helmholtz_data.vel_2_errors.size(), &IterCycle_dimid);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

    // Define coordinate variables
    #if DEBUG>=2
    if (wRank == 0) { fprintf(stdout, "    Defining the dimension variables\n"); }
    #endif
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
    if ( GridType == constants::GridType::MeshGrid ) {
        #if DEBUG>=2
        if (wRank == 0) { fprintf(stdout, "    Write the cell areas\n"); }
        #endif
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
    } else if ( GridType == constants::GridType::LLC ) {
        #if DEBUG>=2
        if (wRank == 0) { fprintf(stdout, "    Write the cell areas\n"); }
        #endif
        int area_dimids[1];
        area_dimids[0] = lat_dimid;
        int area_varid;
        retval = nc_def_var(ncid, "cell_areas", NC_DOUBLE, 1, area_dimids, &area_varid);
        if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

        size_t area_start[1], area_count[1];
        area_start[0] = 0;
        area_count[0] = Nlat;
        retval = nc_put_vara_double(ncid, area_varid, area_start, area_count, &areas[0]);
        if (retval) { NC_ERR(retval, __LINE__, __FILE__); }
    }

    // Close the file
    retval = nc_close(ncid);
    if (retval) { NC_ERR(retval, __LINE__, __FILE__); }

    #if DEBUG >= 2
    if (wRank == 0) { fprintf(stdout, "\nFile (%s) initialized.\n", filename.c_str() ); }
    #endif

    std::vector<std::string> vars = { "scalar" };
    if (wRank == 0) {
        // Loop through and add the desired variables
        // Dimension names (in order!)
        const char* dim_names_MeshGrid[] = {"time", "depth", "latitude", "longitude"};
        const char* dim_names_LLC[]      = {"time", "depth", "latlon"};
        const int ndims_MeshGrid = 4;
        const int ndims_LLC      = 3;
        for (size_t varInd = 0; varInd < vars.size(); ++varInd) {
            if ( GridType == constants::GridType::MeshGrid ) {
                add_var_to_file( vars.at(varInd), dim_names_MeshGrid, ndims_MeshGrid, filename );
            } else if ( GridType == constants::GridType::LLC ) {
                add_var_to_file( vars.at(varInd), dim_names_LLC, ndims_LLC, filename );
            }
        }
    }

    size_t starts[ ( GridType == constants::GridType::MeshGrid ) ? 4 : 3 ],
           counts[ ( GridType == constants::GridType::MeshGrid ) ? 4 : 3 ];
    starts[0] = size_t(myStarts.at(0));
    counts[0] = size_t(myCounts.at(0));

    starts[1] = size_t(myStarts.at(1));
    counts[1] = size_t(myCounts.at(1));

    starts[2] = size_t(myStarts.at(2));
    counts[2] = size_t(myCounts.at(2));

    if ( GridType == constants::GridType::MeshGrid ) {
        starts[3] = size_t(myStarts.at(3));
        counts[3] = size_t(myCounts.at(3));
    }
    write_field_to_output( solution_data.variables.at("land_filled_scalar"), "scalar", starts, counts, filename );


    // And add some attributes for reference
    add_attr_to_file("tolerance", Helmholtz_data.tolerance, filename);
    add_attr_to_file("iterations_per_cycle", Helmholtz_data.iterations_per_cycle, filename);
    add_attr_to_file("iteration_max", Helmholtz_data.iteration_max, filename);

    #if DEBUG >= 2
    if (wRank == 0) { fprintf(stdout, "\n"); }
    #endif
}

