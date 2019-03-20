#ifndef NETCDF_IO_HPP
#define NETCDF_IO_HPP 1

#include <stdio.h>
#include <stdlib.h>
#include <vector>
//#include <mpi.h>

#include "netcdf.h"
//#include "netcdf_par.h"
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

/*! \brief Read necessary input data into workspace.
 *
 * Read source data (from input.nc) for coarse graining.
 *
 * Current assumptions:
 *    Variables to read: uo (as u_lon), vo (as u_lat)
 *
 *    Dimensions: time, depth, longitude, latitude (in that order)
 *
 * There's no u_r in the current data files, so we'll zero it out.
 */
void read_source(
        std::vector<double> &longitude,
        std::vector<double> &latitude,
        std::vector<double> &time,
        std::vector<double> &depth,
        std::vector<double> &u_r,
        std::vector<double> &u_lon,
        std::vector<double> &u_lat,
        std::vector<double> &rho,
        std::vector<double> &p,
        std::vector<double> &mask);

/*! 
 * \brief Initialize netcdf output file for filtered fields.
 *
 *  Creates filter_output.nc, a netcdf4-formatted file,
 *  which will contain the filtered fields in a 
 *  pseudo-CF convention.
 *
 *  Dimension ordering is:
 *    filter_scale, time, depth, latitude, longitude
 */
void initialize_output_file(
        const std::vector<double> &time,      
        const std::vector<double> &depth, 
        const std::vector<double> &longitude, 
        const std::vector<double> &latitude, 
        const std::vector<double> &scales,    
        const std::vector<double> &mask,
        const char * filename = "filter_output.nc");

/*!
 * \brief Write on scales worth of a single field.
 *
 *  Write a single filter-scale worth of a single field,
 *  called field_name, to the previously initialized netcdf
 *  file "filter_output.nc".
 */
void write_field_to_output(
        const std::vector<double> & field, 
        const char * field_name,
        const size_t * start, 
        const size_t * count,
        const char * filename = "filter_output.nc");

/*!
 *  \brief Read a specific variable from a specific file.
 */
void read_var_from_file(
        std::vector<double> &var,
        const char * var_name,
        const char * filename,
        std::vector<double> *mask );

/*!
 *  \brief Add a new variable to a netcdf file
 *
 *  Note: this declares a new variable, but doesn't write to it.
 */
void add_var_to_file(
        const char * var_name,
        const char * filename,
        const char ** dim_list,
        const int num_dims);

#endif
