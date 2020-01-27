#ifndef PREPROCESS_HPP
#define PREPROCESS_HPP 1

#include <stdio.h>
#include <stdlib.h>
#include "ALGLIB/linalg.h"
#include <mpi.h>
#include <vector>

/*!
 *  \brief Apply ALGLIB interpolation routines to fill in masked areas
 */
void interpolate_over_land(
        std::vector<double> &interp_field,
        const std::vector<double> &field,
        const std::vector<double> &time,
        const std::vector<double> &depth,
        const std::vector<double> &latitude,
        const std::vector<double> &longitude,
        const std::vector<double> &mask);

void interpolate_over_land_from_coast(
        std::vector<double> &interp_field,
        const std::vector<double> &field,
        const std::vector<double> &time,
        const std::vector<double> &depth,
        const std::vector<double> &latitude,
        const std::vector<double> &longitude,
        const std::vector<double> &mask);

void get_coast(
        std::vector<double> &lon_coast,
        std::vector<double> &lat_coast,
        std::vector<double> &field_coast,
        const std::vector<double> &lon_full,
        const std::vector<double> &lat_full,
        const std::vector<double> &field_full,
        const std::vector<double> &mask,
        const int Itime,
        const int Idepth,
        const int Ntime,
        const int Ndepth,
        const int Nlat,
        const int Nlon);

void Apply_Toroidal_Projection(
        std::vector<double> & u_lon,
        std::vector<double> & u_lat,
        const std::vector<double> & time,
        const std::vector<double> & depth,
        const std::vector<double> & latitude,
        const std::vector<double> & longitude,
        const std::vector<double> & mask,
        const std::vector<int>    & myCounts,
        const std::vector<int>    & myStarts,
        const MPI_Comm comm = MPI_COMM_WORLD
        );

void toroidal_vel_from_F(  
        std::vector<double> & vel_lon,
        std::vector<double> & vel_lat,
        const std::vector<double> & F,
        const std::vector<double> & longitude,
        const std::vector<double> & latitude,
        const int Ntime,
        const int Ndepth,
        const int Nlat,
        const int Nlon,
        const std::vector<double> & mask
    );

void toroidal_curl_u_dot_er(
        std::vector<double> & out_arr,
        const std::vector<double> & u_lon,
        const std::vector<double> & u_lat,
        const std::vector<double> & longitude,
        const std::vector<double> & latitude,
        const int Ntime,
        const int Ndepth,
        const int Nlat,
        const int Nlon,
        const std::vector<double> & mask
        );

void toroidal_sparse_Lap(
        alglib::sparsematrix & Lap,
        const std::vector<double> & latitude,
        const std::vector<double> & longitude,
        const int Nlat,
        const int Nlon,
        const std::vector<double> & mask
        );

#endif
