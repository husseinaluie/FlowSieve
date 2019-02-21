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
    double *u_r, *u_lon, *u_lat;
    int Nlon, Nlat, Ntime, Ndepth;

    // For the time being, hard-code the filter scales
    const int Nfilt = 1;
    const double filter_scales [4] = {200e3};

    // Read in source data / get size information
    if (debug) {
        fprintf(stdout, "Reading in source data.\n");
    }
    read_source(
            Nlon, Nlat, Ntime, Ndepth,
            &longitude, &latitude, 
            &time, &depth, 
            &u_r, &u_lon, &u_lat
            );

    // Convert coordinate to radians
    for (int ii = 0; ii < Nlon; ii++) {
        longitude[ii] = longitude[ii] * M_PI / 180;
    }
    for (int ii = 0; ii < Nlat; ii++) {
        latitude[ii] = latitude[ii] * M_PI / 180;
    }
    double dlon = longitude[1] - longitude[0];
    double dlat = latitude[1]  - latitude[0];

    // Compute the area of each 'cell'
    //   which will be necessary for integration
    if (debug) {
        fprintf(stdout, "Converting to cell areas.\n");
    }
    double *areas;
    areas = new double[Nlon * Nlat];
    compute_areas(areas, longitude, latitude, Nlon, Nlat);

    // Now pass the arrays along to the filtering routines
    filtering(u_r, u_lon, u_lat,
              filter_scales, Nfilt,
              dlon, dlat,
              Ntime, Ndepth, Nlon, Nlat,
              areas, 
              time, depth,
              longitude, latitude);

    // Done!
    fprintf(stdout, "\n\n");
    fprintf(stdout, "Process completed.\n");
    return 0;

}
