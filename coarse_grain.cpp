#include <fenv.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <algorithm>
#include <math.h>
#include <vector>
//#include <mpi.h>
#include <omp.h>

#include "netcdf_io.hpp"
#include "functions.hpp"

#ifndef DEBUG
    #define DEBUG 0
#endif

int main(int argc, char *argv[]) {

    #if DEBUG >= 0
    fprintf(stdout, "\n\n");
    fprintf(stdout, "Compiled at %s on %s.\n", __TIME__, __DATE__);
    fprintf(stdout, "  Version %d.%d.%d \n\n", MAJOR_VERSION, MINOR_VERSION, PATCH_VERSION);
    fprintf(stdout, "\n\n");
    #endif

    // Enable all floating point exceptions but FE_INEXACT
    //feenableexcept(FE_ALL_EXCEPT & ~FE_INEXACT);

    //MPI_Init( &argc, &argv);

    // Print processor assignments
    int tid, nthreads;
    const int max_threads = omp_get_max_threads();
    omp_set_num_threads( max_threads );
    #pragma omp parallel default(none) private(tid, nthreads) shared(stdout)
    {
        tid = omp_get_thread_num();
        nthreads = omp_get_num_threads();
        #if DEBUG >= 1
        fprintf(stdout, "Hello from thread %d of %d.\n", tid, nthreads);
        #endif
    }

    std::vector<double> longitude;
    std::vector<double> latitude;
    std::vector<double> time;
    std::vector<double> depth;

    std::vector<double> u_r;
    std::vector<double> u_lon;
    std::vector<double> u_lat;

    std::vector<double> rho;
    std::vector<double> p;

    std::vector<double> mask;

    // For the time being, hard-code the filter scales
    //   add zero as the bottom
    /*
    const int Nfilt = 15;
    const double scales [Nfilt+1] = 
            {100e3, 150e3, 200e3, 250e3, 
                300e3, 350e3, 400e3, 450e3, 
                500e3, 750e3, 1000e3,
                1250e3, 1500e3, 1750e3, 
                2000e3, 0};
    */
    /*
    const int Nfilt = 10;
    const double scales [Nfilt+1] = {100e3, 150e3, 200e3, 
        250e3, 300e3, 350e3, 400e3, 450e3, 500e3, 750e3, 0}; 
    */
    const int Nfilt = 1;
    const double scales [Nfilt+1] = 
            {100e3, 0};

    std::vector<double> filter_scales;
    filter_scales.assign(scales, scales + Nfilt);

    // Read in source data / get size information
    #if DEBUG >= 1
    fprintf(stdout, "Reading in source data.\n\n");
    #endif

    read_source(longitude, latitude, time,  depth, 
                u_r, u_lon, u_lat, 
                rho, p,
                mask );

    const int Nlon   = longitude.size();
    const int Nlat   = latitude.size();
    //const int Ntime  = time.size();
    //const int Ndetph = depth.size();

    // Convert coordinate to radians
    int ii;
    #pragma omp parallel default(none) private(ii) shared(longitude)
    { 
        #pragma omp for collapse(1) schedule(dynamic)
        for (ii = 0; ii < Nlon; ii++) {
            longitude.at(ii) = longitude.at(ii) * M_PI / 180;
        }
    }
    #pragma omp parallel default(none) private(ii) shared(latitude)
    {
        #pragma omp for collapse(1) schedule(dynamic)
        for (ii = 0; ii < Nlat; ii++) {
            latitude.at(ii) = latitude.at(ii) * M_PI / 180;
        }
    }

    // Compute the area of each 'cell'
    //   which will be necessary for integration
    #if DEBUG >= 1
    fprintf(stdout, "Computing the cell areas.\n\n");
    #endif

    std::vector<double> areas(Nlon * Nlat);
    compute_areas(areas, longitude, latitude);

    // Now pass the arrays along to the filtering routines
    filtering(u_r, u_lon, u_lat,
              rho, p,
              filter_scales,
              areas, 
              time, depth,
              longitude, latitude,
              mask);


    // Done!
    #if DEBUG >= 0
    fprintf(stdout, "\n\n");
    fprintf(stdout, "Process completed.\n");
    #endif
    return 0;

}
