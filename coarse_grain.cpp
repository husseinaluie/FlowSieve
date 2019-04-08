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
#include "constants.hpp"

int main(int argc, char *argv[]) {

    // For the time being, hard-code the filter scales
    //   add zero as the end for padding
    // Scale must be increasing (ignoring the end)
    // A zero scale will cause everything to nan out (again ignoring the end)
    const int Nfilt = 2;
    const double scales [Nfilt+1] = 
            {100e3, 200e3, 0};

    std::vector<double> filter_scales;
    filter_scales.assign(scales, scales + Nfilt);

    // Parse command-line flags
    char buffer [50];
    for(int ii = 1; ii < argc; ++ii) {  
        fprintf(stdout, "Argument %d : %s\n", ii, argv[ii]);
        snprintf(buffer, 50, "--version");
        if (strcmp(argv[ii], buffer) == 0) {
            print_compile_info(filter_scales);
            return 0;
        } else {
            fprintf(stderr, "Flag %s not recognized.\n", argv[ii]);
            return -1;
        }
    }

    #if DEBUG >= 0
    fprintf(stdout, "\n\n");
    fprintf(stdout, "Compiled at %s on %s.\n", __TIME__, __DATE__);
    fprintf(stdout, "  Version %d.%d.%d \n", MAJOR_VERSION, MINOR_VERSION, PATCH_VERSION);
    fprintf(stdout, "\n");
    #endif

    #if DEBUG >= 0
    #if not(CARTESIAN)
    fprintf(stdout, "Using spherical coordinates.\n");
    #elif CARTESIAN
    fprintf(stdout, "Using Cartesian coordinates.\n");
    #endif
    #endif

    // Enable all floating point exceptions but FE_INEXACT
    //feenableexcept(FE_ALL_EXCEPT & ~FE_INEXACT);

    //MPI_Init( &argc, &argv);

    // Print processor assignments
    int tid, nthreads;
    size_t II;
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

    // Read in source data / get size information
    #if DEBUG >= 1
    fprintf(stdout, "Reading in source data.\n\n");
    #endif

    // Read in the grid coordinates
    read_var_from_file(longitude, "longitude", "input.nc", NULL);
    read_var_from_file(latitude,  "latitude",  "input.nc", NULL);
    read_var_from_file(time,      "time",      "input.nc", NULL);
    read_var_from_file(depth,     "depth",     "input.nc", NULL);

    // Read in the velocity fields
    read_var_from_file(u_lon, "uo", "input.nc", &mask);
    read_var_from_file(u_lat, "vo", "input.nc", &mask);

    // No u_r in inputs, so initialize as zero
    u_r.resize(u_lon.size());
    #pragma omp parallel default(none) private(II) shared(u_r)
    { 
        #pragma omp for collapse(1) schedule(static)
        for (II = 0; II < u_r.size(); II++) {
            u_r.at(II) = 0.;
        }
    }


    #if COMP_BC_TRANSFERS
    // If desired, read in rho and p
    read_var_from_file(rho, "rho", "input.nc", NULL);
    read_var_from_file(p,   "p",   "input.nc", NULL);
    #endif

    const int Nlon   = longitude.size();
    const int Nlat   = latitude.size();

    #if not(CARTESIAN)
    // Convert coordinate to radians
    int ii;
    #pragma omp parallel default(none) private(ii) shared(longitude)
    { 
        #pragma omp for collapse(1) schedule(static)
        for (ii = 0; ii < Nlon; ii++) {
            longitude.at(ii) = longitude.at(ii) * M_PI / 180;
        }
    }
    #pragma omp parallel default(none) private(ii) shared(latitude)
    {
        #pragma omp for collapse(1) schedule(static)
        for (ii = 0; ii < Nlat; ii++) {
            latitude.at(ii) = latitude.at(ii) * M_PI / 180;
        }
    }
    #endif

    // Compute the area of each 'cell'
    //   which will be necessary for integration
    #if DEBUG >= 1
    fprintf(stdout, "Computing the cell areas.\n\n");
    #endif

    std::vector<double> areas(Nlon * Nlat);
    compute_areas(areas, longitude, latitude);

    /*
    #if COMP_BC_TRANSFERS
    const int Ntime  = time.size();
    const int Ndetph = depth.size();
    std::vector<double> rho_mean(Ntime * Ndepth);
    compute_mean(rho_mean, rho, areas, 
            Ntime, Ndepth, Nlat, Nlon,
            mask);
    
    // Subtract the density mean
    for (int Itime = 0; Itime < Ntime; Itime++) {
        for (int Idepth = 0; Idepth < Ndepth; Idepth++) {
            for (int Ilat = 0; Ilat < Nlat; Ilat++) {
                for (int Ilon = 0; Ilon < Nlon; Ilon++) {
                    index = Index(
                            Itime, Idepth, Ilat, Ilon,
                            Ntime, Ndepth, Nlat, Nlon)
                    rho.at(index) = rho.at(index) - rho_mean(Itime*Ndepth + Idepth);
                }
            }
        }
    }
    #endif
    */

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
