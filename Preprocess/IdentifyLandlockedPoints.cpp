#include "../constants.hpp"
#include "../functions.hpp"
#include "../preprocess.hpp"
#include <algorithm>
#include <vector>
#include <omp.h>
#include <math.h>

// Identify which points have land-only neighbours
void HelmholtzDataClass::IdentifyLandlockedPoints(
        const dataset & data,
        const bool second_order_adjacency
        ) {

    size_t Ipt, neighbour_ind;
    const size_t Npts = data.mask.size();
    all_land_neighbours.resize(Npts, false);
    num_coastal = 0;
    num_all_land = 0;

    if (constants::GRID_TYPE == constants::GridType::MeshGrid) {

        int Itime, Idepth, Ilat, Ilon;
        const int Ntime  = data.Ntime,
                  Ndepth = data.Ndepth,
                  Nlat   = data.Nlat,
                  Nlon   = data.Nlon;


        #pragma omp parallel default(none) \
        shared( data, all_land_neighbours ) \
        private( Ipt, Itime, Idepth, Ilat, Ilon, neighbour_ind ) \
        firstprivate( Npts, Ntime, Ndepth, Nlat, Nlon ) \
        reduction( +:num_coastal,num_all_land )
        {
            #pragma omp for collapse(1) schedule(static)
            for (Ipt = 0; Ipt < Npts; Ipt++) {
                all_land_neighbours[Ipt] = true;
                if (data.mask[Ipt]) { all_land_neighbours[Ipt] = false; }
                Index1to4( Ipt, Itime, Idepth, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon );

                neighbour_ind = Index( 0, 0, Ilat, ((Ilon - 1)+Nlon)%Nlon, 1, 1, Nlat, Nlon );
                if ( data.mask[neighbour_ind] ) { all_land_neighbours[Ipt] = false; }

                neighbour_ind = Index( 0, 0, Ilat, ((Ilon + 1)+Nlon)%Nlon, 1, 1, Nlat, Nlon );
                if ( data.mask[neighbour_ind] ) { all_land_neighbours[Ipt] = false; }

                if ( Ilat > 0 ) {
                    neighbour_ind = Index( 0, 0, Ilat-1, Ilon, 1, 1, Nlat, Nlon );
                    if ( data.mask[neighbour_ind] ) { all_land_neighbours[Ipt] = false; }
                }

                if ( Ilat < Nlat - 1 ) {
                    neighbour_ind = Index( 0, 0, Ilat+1, Ilon, 1, 1, Nlat, Nlon );
                    if ( data.mask[neighbour_ind] ) { all_land_neighbours[Ipt] = false; }
                }

                if ( ( all_land_neighbours[Ipt] == 0 ) and ( not(data.mask[Ipt]) ) ) {
                    num_coastal++;
                }
                if ( all_land_neighbours[Ipt] == 1 ) {
                    num_all_land++;
                }
            }
        }
    } else if (constants::GRID_TYPE == constants::GridType::LLC) {
        const size_t num_neighbours = data.num_neighbours;

        for (Ipt = 0; Ipt < Npts; Ipt++) {
            bool is_water = data.mask[Ipt];

            all_land_neighbours[Ipt] = true;
            if ( is_water ) { all_land_neighbours[Ipt] = false; }
            for ( size_t n_I = 0; n_I < num_neighbours; n_I++ ) {
                neighbour_ind = data.adjacency_indices.at(Ipt)[n_I];
                if ( second_order_adjacency ) {
                    // If we're using second order adjacency, then we need to look
                    // at our neighbours neighbours
                    for ( size_t n_J = 0; n_J < num_neighbours; n_J++ ) {
                        size_t neigh_neigh = data.adjacency_indices.at(neighbour_ind)[n_J];
                        if ( data.mask[ neigh_neigh ] ) {
                            all_land_neighbours[Ipt] = false;
                        }
                    }
                } else {
                    // Otherwise, just look at neighbours
                    if ( data.mask[ neighbour_ind ] ) {
                        all_land_neighbours[Ipt] = false;
                    }
                }
            }

            // Next, check if we're coastal land (i.e. we're land, but not all neighbours are)
            if ( ( all_land_neighbours[Ipt] == false ) and not(is_water) ) {
                num_coastal++;
            }
        }
    } else {
        // raise an error
        assert(false);
    }

}
