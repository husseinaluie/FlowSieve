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
 */
void NC_ERR(const int e, const int line_num, const char* file_name);


/*! 
 * \brief Initialize netcdf output file for filtered fields.
 *
 *  Given the filter scale (filter_scale) creates filter_###km.nc, 
 *  a netcdf4-formatted file, (### is formatted to print the filter scale)
 *  which will contain the filtered fields following the CF-convention.
 *
 *  Dimension ordering is:
 *    time, depth, latitude, longitude
 *
 *  vars is a list of variables to initialize when creating the file
 *
 *  filter_scale and rho0 are included as global attributes
 *
 *  The output longitude and latitude fields are given a scale factor
 *    to convert from radians to degrees
 *
 */
void initialize_output_file(
        const std::vector<double> & time,
        const std::vector<double> & depth,
        const std::vector<double> & longitude,
        const std::vector<double> & latitude,
        const std::vector<double> & mask,
        const std::vector<std::string> & vars,
        const char * filename,
        const double filter_scale);

/*!
 * \brief Write on scales worth of a single field.
 *
 *  Write a single filter-scale worth of a single field,
 *  called field_name, to the previously initialized file,
 *  called filename.
 *
 *  start and count correspond to the netcdf put_var arguments 
 *    of the same name
 */
void write_field_to_output(
        const std::vector<double> & field, 
        const char * field_name,
        const size_t * start, 
        const size_t * count,
        const char * filename,
        MPI_Comm = MPI_COMM_WORLD);

/*!
 *  \brief Read a specific variable from a specific file.
 *
 *  Accounts for variable attributes 'scale_factor' and 'add_offset'.
 *
 *  If mask != NULL, then determine the mask based on variable attribute '_FillValue'
 */
void read_var_from_file(
        std::vector<double> &var,
        const char * var_name,
        const char * filename,
        std::vector<double> *mask = NULL,
        std::vector<int> *myCounts = NULL,
        std::vector<int> *myStarts = NULL,
        const MPI_Comm = MPI_COMM_WORLD );

/*!
 *  \brief Add a new variable to a netcdf file
 *
 *  Note: this declares a new variable, but doesn't write to it.
 *  Adds the variable attribute '_FillValue' given as fill_value in constants.hpp
 */
void add_var_to_file(
        const std::string var_name,
        const char ** dim_list,
        const int num_dims,
        const char * filename);
        //const char * filename = "filter_output.nc");

#endif
