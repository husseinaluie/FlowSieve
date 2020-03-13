#ifndef NETCDF_IO_HPP
#define NETCDF_IO_HPP 1

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <string>
#include <mpi.h>

#include "netcdf.h"
#include "netcdf_par.h"
#include "hdf5.h"

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
 * @param[in] time,depth,longitude,latitude vectors giving the dimensions 
 * @param[in] area                          cell areas
 * @param[in] vars                          name of variables to write
 * @param[in] filename                      name for the output file
 * @param[in] filter_scale                  lengthscale used in the filter
 * @param[in] comm                          MPI Communicator (defaults to MPI_COMM_WORLD)
 *
 */
void initialize_output_file(
        const std::vector<double> & time,
        const std::vector<double> & depth,
        const std::vector<double> & longitude,
        const std::vector<double> & latitude,
        const std::vector<double> & areas,
        const std::vector<std::string> & vars,
        const char * filename,
        const double filter_scale,
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
        const std::vector<double> & time,
        const std::vector<double> & depth,
        const std::vector<double> & latitude,
        const std::vector<double> & longitude,
        const std::vector<std::string> & regions,
        const std::vector<std::string> & int_vars,
        const char * filename,
        const double & filter_scale,
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
        const char * filename,
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
        const char * field_name,
        const size_t * start, 
        const size_t * count,
        const char * filename,
        const std::vector<double> * mask = NULL,
        MPI_Comm = MPI_COMM_WORLD
        );


void write_integral_to_post(
        const std::vector<
            std::vector<double> > & field,
        //const char * field_name,
        std::string field_name,
        std::string field_suffix,
        size_t * start,
        size_t * count,
        const char * filename,
        const MPI_Comm comm = MPI_COMM_WORLD
        );


void write_time_average_to_post(
        const std::vector< double > & field,
        std::string field_name,
        std::string field_suffix,
        size_t * start,
        size_t * count,
        const char * filename,
        const MPI_Comm comm = MPI_COMM_WORLD
        );


void write_regions_to_post(
        const char * filename,
        const MPI_Comm comm = MPI_COMM_WORLD
        );

/*!
 *  \brief Read a specific variable from a specific file.
 *
 *  Accounts for variable attributes 'scale_factor' and 'add_offset'.
 *
 *  If mask != NULL, then determine the mask based on variable attribute '_FillValue'
 *
 *  @param[in,out]  var         vector into which to store the loaded variable
 *  @param[in]      var_name    name of the variable to be read
 *  @param[in]      filename    name of the file from which to load the variable
 *  @param[in,out]  mask        point to where a mask array should be stored (if not NULL)
 *  @param[in,out]  myCounts    the sizes of each dimension (on this MPI process) if not NULL
 *  @param[in,out]  myStarts    the starting index for each dimension, if not NULL
 *  @param[in]      do_splits   boolean indicating if the arrays should be split over MPI procs.
 *  @param[in]      comm        the MPI communicator world
 *
 */
void read_var_from_file(
        std::vector<double> &var,
        const char * var_name,
        const char * filename,
        std::vector<double> *mask = NULL,
        std::vector<int> *myCounts = NULL,
        std::vector<int> *myStarts = NULL,
        const bool do_splits = true,
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
        const char * filename,
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
        const std::vector<double> * mask,
        const MPI_Comm comm = MPI_COMM_WORLD
        );

#endif
