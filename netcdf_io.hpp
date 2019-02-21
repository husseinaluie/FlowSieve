#ifndef NETCDF_IO_HPP
#define NETCDF_IO_HPP 1

#include <stdio.h>
#include <stdlib.h>
//#include <mpi.h>

#include "netcdf.h"
//#include "netcdf_par.h"
#include "hdf5.h"

void NC_ERR(const int e, const int line_num, const char* file_name);

/*
void write_filtered_field(int iter, double curr_time, 
        double * my_x, double * my_y,
        double * my_u, double * my_v, double * my_h, double * my_vort,
        int bRank, int my_rank_x, int my_rank_y, 
        int my_Nx, int my_Ny, int full_Nx, int full_Ny);
*/

void read_source(
        int & Nlon,          int & Nlat,
        int & Ntime,         int & Ndepth,
        double ** longitude, double ** latitude,
        double ** time,      double ** depth,
        double ** u_lon,     double ** u_lat);

#endif
