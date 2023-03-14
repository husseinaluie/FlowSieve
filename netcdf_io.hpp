#ifndef NETCDF_IO_HPP
#define NETCDF_IO_HPP 1

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <string>
#include <mpi.h>
#include <math.h>

#include "netcdf.h"
#include "netcdf_par.h"
#include "hdf5.h"

#include "functions.hpp" // this provides class info

/*!
 * \file
 * \brief Collection of functions for handling netcdf IO.
 */


/*! 
 * \brief Wrapper to handle netcdf errors.
 *
 * Short-hand function to handle netcdf errors
 *
 * Prints: the netcdf error
 *         the line number at which the error occurred
 *         the file in which the error occurred
 *
 * @param[in] e         netcdf error code
 * @param[in] line_num  line at which the error occurred
 * @param[in] file_name name of the file in which the error occurred
 *
 */
void NC_ERR(
        const int e, 
        const int line_num, 
        const char* file_name
        );


/*!
 * \brief check if a file exists
 *
 * @param[in] name      Filename to test
 *
 */
bool check_file_existence (
        const std::string& name
        );


/*! 
 * \brief Initialize netcdf output file for filtered fields.
 *
 *  Given the filter scale (filter_scale) creates filter_###km.nc, 
 *  a netcdf4-formatted file, (### is formatted to print the filter scale)
 *  which will contain the filtered fields following the CF-convention.
 *
 *  Dimension ordering is:
 *    time, depth, longitude, latitude
 *
 *  vars is a list of variables to initialize when creating the file
 *
 *  filter_scale and rho0 are included as global attributes
 *
 *  The output longitude and latitude fields are given a scale factor
 *    to convert from radians to degrees
 *
 * @param[in] source_data                   dataset class storing dimension information (as well as other stuff)
 * @param[in] vars                          name of variables to write
 * @param[in] filename                      name for the output file
 * @param[in] filter_scale                  lengthscale used in the filter
 * @param[in] comm                          MPI Communicator (defaults to MPI_COMM_WORLD)
 *
 */
void initialize_output_file(
        const dataset & source_data,
        const std::vector<std::string> & vars,
        const char * filename,
        const double filter_scale = -1,
        MPI_Comm = MPI_COMM_WORLD
        );

void initialize_subset_file(
        const std::vector<double> & time,
        const std::vector<double> & depth,
        const std::vector<double> & windows,
        const int & Nsamples,
        const std::vector<std::string> & vars,
        const char * filename,
        const double filter_scale,
        MPI_Comm = MPI_COMM_WORLD
        );

/*! 
 * \brief Initialize netcdf output file for the postprocessing results.
 *
 *  The post-processing file contains time- and region- averaged values.
 *
 *  Given the filter scale (filter_scale) creates postprocess_###km.nc, 
 *  a netcdf4-formatted file, (### is formatted to print the filter scale).
 *
 *  Dimension ordering is:
 *    time, depth, longitude, latitude
 *
 *  vars is a list of variables to initialize when creating the file
 *
 *  filter_scale and rho0 are included as global attributes
 *
 * @param[in] time,depth,longitude,latitude vectors giving the dimensions 
 * @param[in] regions                       names of the regions used in integration
 * @param[in] int_vars                      names of the variables to write
 * @param[in] filename                      name for the output file
 * @param[in] filter_scale                  lengthscale used in the filter
 * @param[in] comm                          MPI Communicator (defaults to MPI_COMM_WORLD)
 *
 */
void initialize_postprocess_file(
        const dataset & source_data,
        const std::vector<double> & OkuboWeiss_dim_vals,
        const std::vector<std::string> & int_vars,
        const char * filename,
        const double & filter_scale,
        const bool include_OkuboWeiss,
        const MPI_Comm comm = MPI_COMM_WORLD
        );



/*!
 * \brief Create a file that contains the region definitions.
 *
 * Automatically created when post-processing is applied.
 *
 * @param[in] latitude  longitude vector (1D)
 * @param[in] longitude latitude vector (1D)
 * @param[in] regions   name of regions
 * @param[in] filename  name of created file
 * @param[in] comm      MPI communicator (defaults to MPI_COMM_WORLD)
 *
 */
void initialize_regions_file(
        const std::vector<double> & latitude,
        const std::vector<double> & longitude,
        const char * filename,
        const MPI_Comm comm = MPI_COMM_WORLD
        );

void initialize_particle_file(
        const std::vector<double> & time,
        const std::vector<double> & trajectory,
        std::vector<std::string> & vars,
        const std::string & filename,
        const MPI_Comm comm = MPI_COMM_WORLD
        );

void initialize_projected_particle_file(
        const std::vector<double> & time,
        const std::vector<double> & trajectory,
        std::vector<std::string> & vars,
        const std::string & filename,
        const MPI_Comm comm = MPI_COMM_WORLD
        );

/*!
 * \brief Write on scales worth of a single field.
 *
 *  Write a single filter-scale worth of a single field,
 *  called field_name, to the previously initialized file,
 *  called filename.
 *
 *  start and count correspond to the netcdf put_var arguments 
 *    of the same name
 *
 * @param[in] field data    to be written to the file
 * @param[in] field_name    name of the variable in the netcdf file
 * @param[in] start         starting indices for the write
 * @param[in] count         size of the write in each dimension
 * @param[in] filename      name of the netcdf file
 * @param[in] mask          (pointer to) mask that distinguishes land/water cells (default is NULL)
 * @param[in] comm          MPI Communicator
 *
 */
void write_field_to_output(
        const std::vector<double> & field, 
        const std::string & field_name,
        const size_t * start, 
        const size_t * count,
        const std::string & filename,
        const std::vector<bool> * mask = NULL,
        MPI_Comm = MPI_COMM_WORLD
        );


void write_integral_to_post(
        const std::vector<
            std::vector<double> > & field,
        std::string field_name,
        std::string field_suffix,
        size_t * start,
        size_t * count,
        const char * filename,
        const int region_dim,
        const MPI_Comm comm = MPI_COMM_WORLD
        );


void write_time_average_to_post(
        const std::vector< double > & field,
        std::string field_name,
        std::string field_suffix,
        size_t * start,
        size_t * count,
        const char * filename,
        const std::vector<bool> * mask,
        const MPI_Comm comm = MPI_COMM_WORLD
        );


void write_regions_to_post(
        const char * filename,
        const std::vector< std::string > & region_names,
        const MPI_Comm comm = MPI_COMM_WORLD
        );

void read_var_from_file(
        std::vector<double> &var,
        const std::string & var_name,
        const std::string & filename,
        std::vector<bool> *mask = NULL,
        std::vector<int> *myCounts = NULL,
        std::vector<int> *myStarts = NULL,
        const int Nprocs_in_time = 1,
        const int Nprocs_in_depth = 1,
        const bool do_splits = true,
        const int force_split_dim = -1,
        const double land_fill_value = 0.,
        const MPI_Comm = MPI_COMM_WORLD 
        );

void read_mask_from_file(
        std::vector<bool> &mask,
        const std::string & var_name,
        const std::string & filename,
        const int Nprocs_in_time = 1,
        const int Nprocs_in_depth = 1,
        const bool do_splits = true,
        const int force_split_dim = -1,
        const double land_fill_value = 0.,
        const MPI_Comm = MPI_COMM_WORLD 
        );



/*!
 *  \brief Read an attribute from a netcdf file
 *
 *  @param[in,out]  attr        var into which to store the variable
 *  @param[in]      attr_name   Name of the attribute
 *  @param[in]      filename    Name of the file
 *  @param[in]      var_name    Name of associated variable
 *  @param[in]      comm        MPI Communicator
 */
void read_attr_from_file(
        double &attr,
        const char * attr_name,
        const std::string filename,
        const char * var_name = NULL,
        const MPI_Comm comm = MPI_COMM_WORLD 
        );


/*!
 *  \brief Add a new variable to a netcdf file
 *
 *  Note: this declares a new variable, but doesn't write to it.
 *  Adds the variable attribute '_FillValue' given as fill_value in constants.hpp
 *
 *  @param[in] var_name name of the variable to add
 *  @param[in] dim_list dimensions of the variable (in order)
 *  @param[in] num_dims number of dimensions of the variable
 *  @param[in] filename file name to add the variable
 */
void add_var_to_file(
        const std::string var_name,
        const char ** dim_list,
        const int num_dims,
        const char * filename
        );


/*!
 *  \brief Add a new (double) attribute to a netcdf file
 *
 *  This assumes that the attribute is a double
 *
 *  @param[in] varname      name of the attribute
 *  @param[in] value        value of the attribute
 *  @param[in] filename     name of the file to modify
 *  @param[in] comm         MPI communicator (defaults to MPI_COMM_WORLD)
 *
 */
void add_attr_to_file(
        const char * varname,
        const double value,
        const char * filename,
        const MPI_Comm comm = MPI_COMM_WORLD
        );


/*!
 *  \brief Convert a vector<double> to a vector<signed short> for output
 *
 *  This decreases the storage size by a factor of 2 (with, of course, a
 *  corresponding decrease in precision).
 *
 */
void package_field(
        std::vector<signed short> & packaged,
        double & scale_factor,
        double & add_offset,
        const std::vector<double> & original,
        const std::vector<bool> * mask,
        const MPI_Comm comm = MPI_COMM_WORLD
        );



// LLC utilities
void read_LLC_latlon_from_file(
        std::vector<double> &var,
        const std::string & var_name,
        const std::string & filename,
        const MPI_Comm comm = MPI_COMM_WORLD
        );

void initialize_adjacency_file(
        const dataset & source_data,
        const std::vector<std::string> & vars,
        const char * filename,
        const double filter_scale = -1,
        const MPI_Comm comm = MPI_COMM_WORLD
        );

#endif
