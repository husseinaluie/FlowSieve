#include <math.h>
#include <mpi.h>
#include <omp.h>
#include <vector>

#include "../constants.hpp"
#include "../functions.hpp"
#include "../postprocess.hpp"
#include "../netcdf_io.hpp"

void write_regions(
        const char * filename,
        const std::vector<double> & latitude,
        const std::vector<double> & longitude,
        const std::vector<bool>   & mask,
        const std::vector<double> & areas,
        const std::vector<int>    & myCounts,
        const std::vector<int>    & myStarts,
        const MPI_Comm comm
        ) {

    const int Nlat   = myCounts.at(2);
    const int Nlon   = myCounts.at(3);

    const int Slat   = myStarts.at(2);
    const int Slon   = myStarts.at(3);

    int wRank=-1, wSize=-1;
    MPI_Comm_rank( MPI_COMM_WORLD, &wRank );
    MPI_Comm_size( MPI_COMM_WORLD, &wSize );

    // Filename
    initialize_regions_file( latitude, longitude, filename );

    const int Nregion = RegionTest::all_regions.size();
    const int Npts = Nregion * Nlat * Nlon;

    // we need to make a new mask to have the correct shape
    std::vector<double> region_vals(Npts, 0.);
    
    #if DEBUG >= 1
    if (wRank == 0) {
        fprintf(stdout, "Preparing to compute region values\n");
        fflush(stdout);
    }
    #endif
    int index, Ilat, Ilon, Iregion;
    double curr_lat, curr_lon, reg_test;
    #pragma omp parallel \
    default(none) shared(region_vals, latitude, longitude, mask) \
    private(index, Ilat, Ilon, Iregion, curr_lat, curr_lon, reg_test)
    {
        #pragma omp for collapse(3) schedule(static)
        for (Iregion = 0; Iregion < Nregion; ++Iregion) {
            for (Ilat = 0; Ilat < Nlat; ++Ilat) {
                for (Ilon = 0; Ilon < Nlon; ++Ilon) {

                    index = Index(0, Iregion, Ilat, Ilon,
                                  1, Nregion, Nlat, Nlon);

                    curr_lat = latitude.at(Ilat);
                    curr_lon = longitude.at(Ilon);

                    reg_test = (RegionTest::all_regions.at(Iregion))( curr_lat, curr_lon );

                    region_vals.at(index) = reg_test ? 1 : 0;

                }
            }
        }
    }

    #if DEBUG >= 1
    if (wRank == 0) {
        fprintf(stdout, "Preparing to write the region values\n");
        fflush(stdout);
    }
    #endif

    // Write the regions
    //   dimension order: region - lat - lon
    size_t start[3], count[3];
    start[0] = 0;
    count[0] = Nregion;

    start[1] = Slat;
    count[1] = Nlat;

    start[2] = Slon;
    count[2] = Nlon;

    write_field_to_output(region_vals, "region_definitions", 
            start, count, filename, NULL);

}
