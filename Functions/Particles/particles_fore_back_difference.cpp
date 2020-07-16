#include <math.h>
#include <algorithm>
#include <vector>
#include <mpi.h>
#include <omp.h>
#include "../../constants.hpp"
#include "../../functions.hpp"
#include "../../particles.hpp"

void particles_fore_back_difference(
        std::vector<double> & traj_dists,
        const std::vector<double> & fore_lon_hist,
        const std::vector<double> & fore_lat_hist,
        const std::vector<double> & back_lon_hist,
        const std::vector<double> & back_lat_hist,
        const MPI_Comm comm
        ) {

    int wRank, wSize;
    MPI_Comm_rank( comm, &wRank );
    MPI_Comm_size( comm, &wSize );

    const size_t Npts = fore_lon_hist.size();
    size_t index;

    double lon1, lat1, lon2, lat2;

    #pragma omp parallel \
    default(none) \
    shared( fore_lon_hist, fore_lat_hist, back_lon_hist, back_lat_hist, traj_dists ) \
    private( index, lon1, lat1, lon2, lat2 )
    {
        #pragma omp for collapse(1) schedule(dynamic)
        for (index = 0; index < Npts; ++index) {

            lon1 = fore_lon_hist.at(index);
            lat1 = fore_lat_hist.at(index);
            lon2 = back_lon_hist.at(index);
            lat2 = back_lat_hist.at(index);

            if (     ( lon1 != constants::fill_value )
                 and ( lat1 != constants::fill_value )  
                 and ( lon2 != constants::fill_value )  
                 and ( lat2 != constants::fill_value )  
               ) {
                traj_dists.at(index) = distance( lon1, lat1, lon2, lat2 );
            }
        } 
    }
}
