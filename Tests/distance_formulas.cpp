#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <algorithm>
#include <math.h>
#include <vector>
#include <time.h>
#include <assert.h>
#include "../functions.hpp"
#include "../constants.hpp"
#include "../netcdf_io.hpp"

double distance_simple(
        const double lon1,
        const double lat1,
        const double lon2,
        const double lat2
        ) {

    const double    Delta_lon   = lon2 - lon1,
                    cos_lat2    = cos(lat2),
                    sin_lat2    = sin(lat2),
                    cos_lat1    = cos(lat1),
                    sin_lat1    = sin(lat1),
                    cos_Delta_lon = cos(Delta_lon);

    double Delta_sigma = acos( sin_lat1 * sin_lat2 + cos_lat1 * cos_lat2 * cos_Delta_lon );
    double distance = constants::R_earth * Delta_sigma;

    return distance;
}

double distance_high_precision(
        const double lon1,
        const double lat1,
        const double lon2,
        const double lat2
        ) {

    const double    Delta_lon   = lon2 - lon1,
                    cos_lat2    = cos(lat2),
                    sin_lat2    = sin(lat2),
                    cos_lat1    = cos(lat1),
                    sin_lat1    = sin(lat1),
                    cos_Delta_lon = cos(Delta_lon);
    double numer, denom, Delta_sigma;
    numer =   pow(                                    cos_lat2 * sin(Delta_lon), 2 ) 
            + pow( cos_lat1 * sin_lat2  -  sin_lat1 * cos_lat2 * cos_Delta_lon , 2 );
    numer =  sqrt(numer);

    denom =        sin_lat1 * sin_lat2  +  cos_lat1 * cos_lat2 * cos_Delta_lon ;

    Delta_sigma = atan2(numer, denom);

    double distance = constants::R_earth * Delta_sigma;

    return distance;
}

int main(int argc, char *argv[]) {

    static_assert( not(constants::CARTESIAN), "Test only applicable to Spherical coordinates.\n" );

    fprintf(stdout, "Beginning tests for distance calculation methods.\n");

    int thread_safety_provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &thread_safety_provided);
    MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI::ERRORS_THROW_EXCEPTIONS);

    int wRank=-1, wSize=-1;
    MPI_Comm_rank( MPI_COMM_WORLD, &wRank );
    MPI_Comm_size( MPI_COMM_WORLD, &wSize );

    assert(wSize==1);

    const int Nlat = 180 * 12,
              Nlon = 360 * 12;

    // Reference points
    const double    ref_lat = 0.,
                    ref_lon = 0.;

    std::vector<double> distances_slow( Nlat * Nlon, 0. ),
                        distances_fast( Nlat * Nlon, 0. ),
                        times{ 0. },
                        depth{ 0. },
                        latitude( Nlat ),
                        longitude( Nlon );


    const double    lon_min = 0.,
                    lon_max = 2. * M_PI,
                    lat_min = -M_PI/2.,
                    lat_max =  M_PI/2.;

    const double    dlat = (lat_max - lat_min) / Nlat,
                    dlon = (lon_max - lon_min) / Nlon;

    for (int II = 0; II < Nlat; II++) { latitude.at( II) = lat_min + (II+0.5) * dlat; }
    for (int II = 0; II < Nlon; II++) { longitude.at(II) = lon_min + (II+0.5) * dlon; }

    int index;
    const double slow_distance_start_time = MPI_Wtime();
    for (int Ilat = 0; Ilat < Nlat; ++Ilat) {
        for (int Ilon = 0; Ilon < Nlon; ++Ilon) {
            index = Ilat * Nlon + Ilon;

            distances_slow.at(index) = distance_high_precision( longitude.at(Ilon), latitude.at(Ilat), ref_lon, ref_lat );
        }
    }
    const double slow_distance_stop_time = MPI_Wtime();

    const double fast_distance_start_time = MPI_Wtime();
    for (int Ilat = 0; Ilat < Nlat; ++Ilat) {
        for (int Ilon = 0; Ilon < Nlon; ++Ilon) {
            index = Ilat * Nlon + Ilon;

            distances_fast.at(index) = distance_simple( longitude.at(Ilon), latitude.at(Ilat), ref_lon, ref_lat );
        }
    }
    const double fast_distance_stop_time = MPI_Wtime();


    fprintf( stdout, "Timing Results  (Nlat, Nlon) = (%'d,%'d)\n\n", Nlat, Nlon );
    fprintf( stdout, " Simplified  calculation: %.13g \n", fast_distance_stop_time - fast_distance_start_time );
    fprintf( stdout, " Full stable calculation: %.13g \n", slow_distance_stop_time - slow_distance_start_time );
    fprintf( stdout, "\n" );

    // Now write them to a file
    dataset output_data;
    output_data.time      = times;
    output_data.depth     = depth;
    output_data.latitude  = latitude;
    output_data.longitude = longitude;

    output_data.Ntime   = 1;
    output_data.Ndepth  = 1;
    output_data.Nlat    = Nlat;
    output_data.Nlon    = Nlon;
    
    std::vector<std::string> vars_to_write;
    vars_to_write.push_back("distance_fast");
    vars_to_write.push_back("distance_slow");

    const std::string fname = "distance.nc";
    initialize_output_file( output_data, vars_to_write, fname.c_str() );

    size_t starts[4] = {0,0,0,0};
    size_t counts[4] = {1,1,Nlat,Nlon};

    write_field_to_output( distances_slow, "distance_slow", starts, counts, fname );
    write_field_to_output( distances_fast, "distance_fast", starts, counts, fname );

    MPI_Finalize();

}
