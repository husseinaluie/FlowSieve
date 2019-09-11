#include <fenv.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <algorithm>
#include <math.h>
#include <vector>
#include <mpi.h>
#include <omp.h>
#include <cassert>

#include "../netcdf_io.hpp"
#include "../functions.hpp"
#include "../constants.hpp"
#include "../postprocess.hpp"

int main(int argc, char *argv[]) {
    
    // Specify the number of OpenMP threads
    //   and initialize the MPI world
    int thread_safety_provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &thread_safety_provided);
    MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI::ERRORS_THROW_EXCEPTIONS);
    //MPI_Status status;

    int wRank=-1, wSize=-1;
    MPI_Comm_rank( MPI_COMM_WORLD, &wRank );
    MPI_Comm_size( MPI_COMM_WORLD, &wSize );

    // Parse command-line flags
    char ifile[50];
    char ofile[50];
    assert(argc == 3);
    snprintf(ifile, 50, argv[1]);
    snprintf(ofile, 50, argv[2]);
    if (wRank == 0) {
        fprintf(stdout, "Argument 1 (input file)  : %s\n", argv[1]);
        fprintf(stdout, "Argument 2 (output file) : %s\n", argv[2]);
    }

    // Set thread caps
    const int max_threads = omp_get_max_threads();
    omp_set_num_threads( max_threads );

    // Get filter scale
    double filter_scale;
    read_attr_from_file(filter_scale, "filter_scale", ifile);

    // Read in grids
    std::vector<double> longitude;
    std::vector<double> latitude;
    std::vector<double> time;
    std::vector<double> depth;

    std::vector<double> mask;
    std::vector<int> myCounts, myStarts;

    // Read in source data / get size information
    #if DEBUG >= 1
    if (wRank == 0) { fprintf(stdout, "Reading in source data.\n\n"); }
    #endif

    // Read in the grid coordinates
    read_var_from_file(longitude, "longitude", ifile);
    read_var_from_file(latitude,  "latitude",  ifile);
    read_var_from_file(time,      "time",      ifile);
    read_var_from_file(depth,     "depth",     ifile);

    const int Nlon = longitude.size();
    const int Nlat = latitude.size();
    const double D2R = M_PI / 180.;

    if (not(constants::CARTESIAN)) {
        // Convert coordinate to radians
        if (wRank == 0) { fprintf(stdout, "Converting to radians.\n\n"); }
        int ii;
        #pragma omp parallel default(none) private(ii) shared(longitude, latitude)
        { 
            #pragma omp for collapse(1) schedule(static)
            for (ii = 0; ii < Nlon; ii++) {
                longitude.at(ii) = longitude.at(ii) * D2R;
            }

            #pragma omp for collapse(1) schedule(static)
            for (ii = 0; ii < Nlat; ii++) {
                latitude.at(ii) = latitude.at(ii) * D2R;
            }
        }
    }

    // Compute the area of each 'cell'
    //   which will be necessary for integration
    #if DEBUG >= 1
    if (wRank == 0) { fprintf(stdout, "Computing the cell areas.\n\n"); }
    #endif

    std::vector<double> areas(Nlon * Nlat);
    compute_areas(areas, longitude, latitude);

    // Read in the requisit fields
    std::vector<double> u_r, u_lon, u_lat, Pi;
    read_var_from_file(u_r,   "coarse_u_r",      ifile, &mask, &myCounts, &myStarts);
    read_var_from_file(u_lon, "coarse_u_lon",    ifile, &mask, &myCounts, &myStarts);
    read_var_from_file(u_lat, "coarse_u_lat",    ifile, &mask, &myCounts, &myStarts);
    read_var_from_file(Pi,    "energy_transfer", ifile, &mask, &myCounts, &myStarts);

    // Get some size parameters
    const size_t Npts   = Pi.size();
    const int    Ntime  = myCounts.at(0);
    const int    Ndepth = myCounts.at(1);
    const int    Stime  = myStarts.at(0);
    const int    Sdepth = myStarts.at(1);

    // Load in the region functions
    const std::vector< bool(*)(double,double) > region_tests = RegionTest::all_regions;
    const int num_regions = region_tests.size();
    if (wRank == 0) {
        fprintf(stdout, "Will integrate over %d spatial regions.\n", num_regions);
    }
    std::vector<double> region_areas(num_regions);

    // Domain integrals
    double initial_value = 0;

    std::vector<std::vector<double>> Pi_integrals, KE_integrals;
    Pi_integrals.resize(num_regions, std::vector<double>(Ntime * Ndepth, initial_value));
    KE_integrals.resize(num_regions, std::vector<double>(Ntime * Ndepth, initial_value));

    fprintf(stdout, "%d, %d, %d, %d\n", Ntime, Ndepth, Nlat, Nlon);

    double curr_Pi, curr_KE, dA;
    int Ilat, Ilon, Itime, Idepth, Iregion, area_index, int_index;
    size_t index;
    #pragma omp parallel default(none)\
    private(Ilat, Ilon, Itime, Idepth, Iregion, \
            index, curr_Pi, curr_KE, dA, area_index, int_index)\
    shared(latitude, longitude, Pi, u_lon, u_lat, \
            Pi_integrals, KE_integrals, areas, mask,\
            region_areas) \
    firstprivate(wRank)
    { 
        #pragma omp for collapse(1) schedule(dynamic)
        for (index = 0; index < Npts; index++) {
            if (mask.at(index) == 1) { // Skip land areas

                Index1to4(index, Itime, Idepth, Ilat, Ilon,
                                 Ntime, Ndepth, Nlat, Nlon);

                int_index = Index(Itime, Idepth, 0, 0,
                                  Ntime, Ndepth, 1, 1);

                area_index = Index(0,     0,      Ilat, Ilon,
                                   Ntime, Ndepth, Nlat, Nlon);

                curr_Pi = Pi.at(index);
                curr_KE = (pow(u_lon.at(index), 2) + pow(u_lat.at(index), 2));
                dA      = areas.at(area_index);

                for (Iregion = 0; Iregion < num_regions; ++Iregion) {
                    if ( region_tests.at(Iregion)(latitude.at(Ilat), longitude.at(Ilon)) ) {
                        Pi_integrals.at(Iregion).at(int_index) += curr_Pi * dA;
                        KE_integrals.at(Iregion).at(int_index) += curr_KE * dA;

                        if ( (Itime == 0) and (Idepth == 0) and (wRank == 0)) {
                            region_areas.at(Iregion) += dA;
                        }
                    }
                }

            }
        }
    }

    // Variables to integrate
    std::vector<std::string> integrated_vars;
    integrated_vars.push_back("Pi");
    integrated_vars.push_back("KE");

    // Open the NETCDF file
    char filename [50];
    snprintf(filename, 50, ofile);

    initialize_postprocess_file(
            time, depth, latitude, longitude, 
            RegionTest::region_names,
            integrated_vars,
            filename,
            filter_scale
            );

    // Dimension order: time - depth - region
    size_t start[3], count[3];
    start[0] = Stime;
    count[0] = Ntime;

    start[1] = Sdepth;
    count[1] = Ndepth;

    count[2] = 1;

    write_integral_to_post(Pi_integrals, "Pi", start, count, filename);
    write_integral_to_post(KE_integrals, "KE", start, count, filename);

    // Write the region areas (needed for normalization)
    size_t start_r[1], count_r[1];
    start_r[0] = 0;
    count_r[0] = (size_t) num_regions;
    write_field_to_output(region_areas, "region_areas", start_r, count_r, filename, NULL);

    // Write region names
    //   this has to be done separately for reasons
    write_regions_to_post(filename);

    // DONE!
    fprintf(stdout, "Processor %d / %d waiting to finalize.\n", wRank + 1, wSize);
    MPI_Finalize();
    return 0;
}
