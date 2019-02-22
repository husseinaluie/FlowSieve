#ifndef NETCDF_IO_HPP
#define NETCDF_IO_HPP 1

#include <stdio.h>
#include <stdlib.h>
//#include <mpi.h>

#include "netcdf.h"
//#include "netcdf_par.h"
#include "hdf5.h"

void NC_ERR(const int e, const int line_num, const char* file_name);

void read_source(
        int & Nlon,          int & Nlat,
        int & Ntime,         int & Ndepth,
        double ** longitude, double ** latitude,
        double ** time,      double ** depth,
        double ** u_r,       double ** u_lon,
        double ** u_lat,     double ** mask);

void initialize_output_file(
        const int Ntime, const int Ndepth, const int Nlon, 
                         const int Nlat, const int Nscales,
        const double * time,      const double * depth, 
        const double * longitude, const double * latitude, 
        const double * scales,    const double * mask);

void write_to_output(
        const double * coarse_u_r, const double * coarse_u_lon, const double * coarse_u_lat, 
        const int Iscale, const int Ntime, const int Ndepth,
        const int Nlat,   const int Nlon);

#endif
