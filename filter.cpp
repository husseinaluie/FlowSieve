#include <fenv.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <algorithm>
#include <math.h>
//#include <mpi.h>
//#include <omp.h>

#include "netcdf_io.hpp"
#include "functions.hpp"

#ifndef DEBUG
    #define DEBUG false
#endif

int main(int argc, char *argv[]) {

    const bool debug = DEBUG;

    if (debug) {
        fprintf(stdout, "\n\n");
        fprintf(stdout, "Compiled at %s on %s.\n", __TIME__, __DATE__);
        fprintf(stdout, "\n\n");
    }

    // Enable all floating point exceptions but FE_INEXACT
    //feenableexcept(FE_ALL_EXCEPT & ~FE_INEXACT);

    //MPI_Init( &argc, &argv);

    double *longitude, *latitude, *time, *depth;
    double *u_lon, *u_lat;
    int Nlon, Nlat, Ntime, Ndepth;

    // For the time being, hard-code the filter scales
    const int Nfilt = 4;
    const double filter_scales [4] = {50e3, 150e3, 250e3, 350e3};

    // Read in source data / get size information
    read_source(
            Nlon, Nlat, Ntime, Ndepth,
            &longitude, &latitude, 
            &time, &depth, 
            &u_lon, &u_lat
            );

    // Convert coordinate to radians
    for (int ii = 0; ii < Nlon; ii++) {
        longitude[ii] = longitude[ii] * M_PI / 180;
    }
    for (int ii = 0; ii < Nlat; ii++) {
        latitude[ii] = latitude[ii] * M_PI / 180;
    }


    // Compute the area of each 'cell'
    //   which will be necessary for integration
    double *areas;
    areas = new double[Nlon * Nlat];
    compute_areas(areas, longitude, latitude, Nlon, Nlat);

    // Now convert the Spherical velocities to Cartesian
    //   (although we will still be on a spherical
    //     coordinate system)
    double *u_x, *u_y, *u_z;
    u_x = new double[Ntime * Ndepth * Nlon * Nlat];
    u_y = new double[Ntime * Ndepth * Nlon * Nlat];
    u_z = new double[Ntime * Ndepth * Nlon * Nlat];

    int index;
    for (int Itime = 0; Itime < Ntime; Itime++) {
        for (int Idepth = 0; Idepth < Ndepth; Idepth++) {
            for (int Ilat = 0; Ilat < Nlat; Ilat++) {
                for (int Ilon = 0; Ilon < Nlon; Ilon++) {

                    // Convert our four-index to a one-index
                    index = Index(Itime, Idepth, Ilat, Ilon,
                                  Ntime, Ndepth, Nlat, Nlon);

                    // Note that 0 is in place of u_r, since 
                    // the current datasets don't have u_r
                    vel_Spher_to_Cart(u_x[index], u_y[index],   u_z[index],
                                      0.,         u_lon[index], u_lat[index],
                                      longitude[Ilon], latitude[Ilat]);
                }
            }
        }
    }

    // Done!
    fprintf(stdout, "\n\n");
    fprintf(stdout, "Process completed.\n");
    return 0;

}
