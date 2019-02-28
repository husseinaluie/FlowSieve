#ifndef NETCDF_IO_HPP
#define NETCDF_IO_HPP 1

#include <stdio.h>
#include <stdlib.h>
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
        int & Nlon,          int & Nlat,
        int & Ntime,         int & Ndepth,
        double ** longitude, double ** latitude,
        double ** time,      double ** depth,
        double ** u_r,       double ** u_lon,
        double ** u_lat,     double ** mask);

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
        const int Ntime, const int Ndepth, const int Nlon, 
                         const int Nlat, const int Nscales,
        const double * time,      const double * depth, 
        const double * longitude, const double * latitude, 
        const double * scales,    const double * mask);

/*! \brief Write one scale's worth of outputs to file.
 *
 *  Write a single filter-scale worth of filtered velocities
 *  to the previously initialized netcdf file
 *  "filter_output.nc".
 */
void write_to_output(
        const double * coarse_u_r, 
        const double * coarse_u_lon, 
        const double * coarse_u_lat, 
        const int Iscale, const int Ntime, const int Ndepth,
        const int Nlat,   const int Nlon);

/*! \brief Write one scales worth of vorticity to file
 *
 *  Write a single filter-scale worth of filtered vorticity
 *  to the previously initialized netcdf file
 *  "filter_output.nc".
 */
void write_vorticity(
        const double * vort_r, 
        const double * vort_lon, 
        const double * vort_lat,
        const int Iscale, const int Ntime, const int Ndepth, 
        const int Nlat,   const int Nlon);

/*!
 * \brief Write on scales worth of energy transfers to file
 *
 *  Write a single filter-scale worth of energy transfers
 *  to the previously initialized netcdf file
 *  "filter_output.nc".
 */
void write_energy_transfer(
        const double * energy_transfer, const int Iscale, 
        const int Ntime, const int Ndepth, const int Nlat, const int Nlon);

#endif
