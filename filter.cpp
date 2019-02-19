#include <fenv.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <algorithm>
#include <math.h>
//#include <mpi.h>
//#include <omp.h>
#include "netcdf_io.hpp"

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

    if (debug) {
        fprintf(stdout, "Hello, World!\n");
    }

    double *longitude, *latitude, *time, *depth;
    double *u_lon, *u_lat;

    read_source(&longitude, &latitude, &time, &depth, &u_lon, &u_lat);

    return 0;
}
