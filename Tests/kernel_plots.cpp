#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <algorithm>
#include <math.h>
#include <vector>
#include <assert.h>
#include "../functions.hpp"
#include "../constants.hpp"
#include "../netcdf_io.hpp"

int main(int argc, char *argv[]) {

    int thread_safety_provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &thread_safety_provided);
    MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI::ERRORS_THROW_EXCEPTIONS);

    int wRank=-1, wSize=-1;
    MPI_Comm_rank( MPI_COMM_WORLD, &wRank );
    MPI_Comm_size( MPI_COMM_WORLD, &wSize );

    assert(wSize==1);

    const int Nlat = 180 * 4,
              Nlon = 360 * 4;

    const double    lon_min = 0.,
                    lon_max = 2. * M_PI,
                    lat_min = -M_PI/2.,
                    lat_max =  M_PI/2.;

    const double    D2R = M_PI / 180.;

    // 
    const double filter_scale = 3000e3;

    // Reference points - lon reference will just be uniformly spaced
    const std::vector<double> ref_lats = { -89,  70, -50,  30, -10, 0 };
    const size_t num_kernels = ref_lats.size();

    std::vector<double> kernel_values( num_kernels * Nlat * Nlon ),
                        times( num_kernels ),
                        depth{ 0. },
                        latitude( Nlat ),
                        longitude( Nlon );

    const double    dlat = (lat_max - lat_min) / Nlat,
                    dlon = (lon_max - lon_min) / Nlon;

    for (int II = 0; II < num_kernels;  II++) { times.at(II) = II; }
    for (int II = 0; II < Nlat;         II++) { latitude.at( II) = lat_min + (II+0.5) * dlat; }
    for (int II = 0; II < Nlon;         II++) { longitude.at(II) = lon_min + (II+0.5) * dlon; }

    for (size_t Ikern = 0; Ikern < num_kernels; ++Ikern) {
        double ref_lat = D2R * ref_lats.at(Ikern);
        double ref_lon = lon_min + (lon_max - lon_min) * ( (double)Ikern / num_kernels );

        for (int Ilat = 0; Ilat < Nlat; ++Ilat) {
            for (int Ilon = 0; Ilon < Nlon; ++Ilon) {
                double dist = distance( ref_lon, ref_lat, longitude.at(Ilon), latitude.at(Ilat) );
                double kern = kernel( dist, filter_scale );

                size_t index = Index( 0, Ikern, Ilat, Ilon, 1, num_kernels, Nlat, Nlon );

                kernel_values.at(index) = kern;
            }
        }
    }


    // Now write them to a file
    dataset output_data;
    output_data.time      = times;
    output_data.depth     = depth;
    output_data.latitude  = latitude;
    output_data.longitude = longitude;

    output_data.Ntime   = num_kernels;
    output_data.Ndepth  = 1;
    output_data.Nlat    = Nlat;
    output_data.Nlon    = Nlon;
    
    std::vector<std::string> vars_to_write;
    vars_to_write.push_back("kernel_values");

    const std::string fname = "kernels.nc";
    initialize_output_file( output_data, vars_to_write, fname.c_str() );

    size_t starts[4] = {0,0,0,0};
    size_t counts[4] = {num_kernels,1,Nlat,Nlon};

    write_field_to_output( kernel_values, "kernel_values", starts, counts, fname );

    MPI_Finalize();

}
