#include <math.h>
#include <time.h>
#include <algorithm>
#include <vector>
#include <omp.h>
#include <mpi.h>
#include "../../constants.hpp"
#include "../../functions.hpp"
#include "../../particles.hpp"

void particles_initial_positions(
        std::vector<double> & starting_lat, 
        std::vector<double> & starting_lon, 
        const int Npts,
        const std::vector<double> & latitude, 
        const std::vector<double> & longitude, 
        const std::vector<bool> & mask,
        const MPI_Comm comm
        ) {

    int wRank, wSize;
    MPI_Comm_rank( comm, &wRank );
    MPI_Comm_size( comm, &wSize );

    std::vector<double> lon_rng, lat_rng, lon_mid, lat_mid;

    /*
    lon_rng.push_back(         longitude.back() - longitude.front()  );
    lat_rng.push_back( 0.9 * ( latitude.back()  - latitude.front() ) );
    lon_mid.push_back( 0.5 * ( longitude.back() + longitude.front()) );
    lat_mid.push_back( 0.5 * ( latitude.back()  + latitude.front() ) );
    */
    
    const double D2R = M_PI / 180.;
    // Agulhas
    // plt.xlim( 10,  60)
    // plt.ylim(-40,  20)
    lon_rng.push_back(  50 * D2R );
    lat_rng.push_back(  60 * D2R );
    lon_mid.push_back(  35 * D2R );
    lat_mid.push_back( -10 * D2R );
    
    // Gulf
    // plt.xlim(-100, -30 )
    // plt.ylim(  10,  60 )
    lon_rng.push_back(  70 * D2R );
    lat_rng.push_back(  50 * D2R );
    lon_mid.push_back( -65 * D2R );
    lat_mid.push_back(  35 * D2R );
    
    // ACC
    // plt.xlim(-90, -20)
    // plt.ylim(-75, -15)
    lon_rng.push_back(  70 * D2R );
    lat_rng.push_back(  60 * D2R );
    lon_mid.push_back( -55 * D2R );
    lat_mid.push_back( -45 * D2R );
    
    // Kuroshio
    // plt.xlim(110, 150)
    // plt.ylim(  0,  60)
    lon_rng.push_back(   40 * D2R );
    lat_rng.push_back(   60 * D2R );
    lon_mid.push_back(  130 * D2R );
    lat_mid.push_back(   30 * D2R );
    
    // Kuroshio (spot)
    // plt.xlim(132, 136)
    // plt.ylim( 29,  31)
    /*
    lon_rng.push_back(   4 * D2R );
    lat_rng.push_back(   2 * D2R );
    lon_mid.push_back( 134 * D2R );
    lat_mid.push_back(  30 * D2R );
    */

    int left, right, bottom, top;
    size_t BL_ind, BR_ind, TL_ind, TR_ind;

    const int Nlat = latitude.size(),
              Nlon = longitude.size(),
              Nreg = lon_rng.size();

    double part_lon, part_lat;

    //srand( wRank + time(NULL) );
    srand( wRank );

    for ( int II = 0; II < Npts; ++II ) {

        part_lon = ( ((double) rand() / (RAND_MAX)) - 0.5) * lon_rng.at(II % Nreg) + lon_mid.at(II % Nreg);
        part_lat = ( ((double) rand() / (RAND_MAX)) - 0.5) * lat_rng.at(II % Nreg) + lat_mid.at(II % Nreg);

        // Check if this particle is on land
        particles_get_edges(left, right, bottom, top, 
                part_lat, part_lon, latitude, longitude);

        BL_ind = Index(0, 0, bottom, left,
                       1, 1, Nlat,   Nlon);
        BR_ind = Index(0, 0, bottom, right,
                       1, 1, Nlat,   Nlon);
        TL_ind = Index(0, 0, top,    left,
                       1, 1, Nlat,   Nlon);
        TR_ind = Index(0, 0, top,    right,
                       1, 1, Nlat,   Nlon);

        // So long as we are on land, keep picking a new point
        while ( not( mask.at(BL_ind) and mask.at(TL_ind) and mask.at(BR_ind) and mask.at(TR_ind) )
                or ( std::isnan(bottom) )
                or ( std::isnan(top)    )
              ) {

            part_lon = ( ((double) rand() / (RAND_MAX)) - 0.5) * lon_rng.at(II % Nreg) + lon_mid.at(II % Nreg);
            part_lat = ( ((double) rand() / (RAND_MAX)) - 0.5) * lat_rng.at(II % Nreg) + lat_mid.at(II % Nreg);

            // Check if this particle is on land
            particles_get_edges(left, right, bottom, top, 
                    part_lat, part_lon, latitude, longitude);

            BL_ind = Index(0, 0, bottom, left,
                           1, 1, Nlat,   Nlon);
            BR_ind = Index(0, 0, bottom, right,
                           1, 1, Nlat,   Nlon);
            TL_ind = Index(0, 0, top,    left,
                           1, 1, Nlat,   Nlon);
            TR_ind = Index(0, 0, top,    right,
                           1, 1, Nlat,   Nlon);
        }

        // Store the position
        starting_lat.at(II) = part_lat;
        starting_lon.at(II) = part_lon;

    }
}
