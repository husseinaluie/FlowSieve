#include "../constants.hpp"
#include "../functions.hpp"
#include "../preprocess.hpp"
#include "../differentiation_tools.hpp"
#include <algorithm>
#include <deque>
#include <vector>
#include <omp.h>
#include <math.h>

// Create land-collapsing map [i.e. map continguous land to single point]
//void HelmholtzDataClass::CreateLandCollapsingMap(
void CreateLandCollapsingMap(
        dataset & data
        ) {

    
    const int Ntime  = data.Ntime,
              Ndepth = data.Ndepth,
              Nlat   = data.Nlat,
              Nlon   = data.Nlon;
    const size_t Npts = data.mask.size();

    size_t num_mapped_onto = 0, num_mapped_points = 0, Ipt, Jpt, Ineighbour, Itest;
    std::deque<size_t> points_to_test;

    data.pt_maps_to.resize(Npts, Npts);
    size_t num_land = 0;

    int Itime, Idepth, Ilat, Ilon, Ineigh;
    std::vector<double> diff_vec;


    // Flag points as coastal land
    std::vector<short int> pt_is_coastal(Npts, false);

    for (Ipt = 0; Ipt < Npts; Ipt++) {

        if (data.mask[Ipt]) { 
            // if water, not coastal land
            continue;
        } else if ( pt_coastal[Ipt] ) {
            // already flagged
            continue;
        } else {
            // otherwise, check all neighbours and see if any are water
            bool is_coastal = false;
            if (constants::GRID_TYPE == constants::GridType::MeshGrid) {
                // Get coordinates
                Index1to4( Ipt, Itime, Idepth, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon );

                // check each neighbour
                Ineigh = ( (Ilon+1) % Nlon + Nlon ) % Nlon;
                Jpt = Index(0, 0, Ilat, Ineigh, 1, 1, Nlat, Nlon);
                if ( data.mask[Jpt] ) { is_coastal = true; }

                Ineigh = ( (Ilon-1) % Nlon + Nlon ) % Nlon;
                Jpt = Index(0, 0, Ilat, Ineigh, 1, 1, Nlat, Nlon);
                if ( data.mask[Jpt] ) { is_coastal = true; }

                if (Ilat > 0) {
                    Jpt = Index(0, 0, Ilat-1, Ilon, 1, 1, Nlat, Nlon);
                    if ( data.mask[Jpt] ) { is_coastal = true; }
                }
                if (Ilat < Nlat-1) {
                    Jpt = Index(0, 0, Ilat+1, Ilon, 1, 1, Nlat, Nlon);
                    if ( data.mask[Jpt] ) { is_coastal = true; }
                }
            } else {
                for (size_t neighbour_ind = 0; neighbour_ind < data.num_neighbours; neighbour_ind++ ) {
                    Jpt = data.adjacency_indices.at(Ipt)[neighbour_ind];
                    if ( data.mask[Jpt] ) { is_coastal = true; }
                }
            }
            pt_is_coastal[Jpt] = is_coastal;
        }
    }


    // Now create the land collapse map
    for (Ipt = 0; Ipt < Npts; Ipt++) {

        if (data.mask[Ipt]) { 
            // Water points map to themselves
            data.pt_maps_to[Ipt] = Ipt; 
            continue;
        } else if ( coastal_test[Ipt] ) { 
            // Coastal points map to themselves
            data.pt_maps_to[Ipt] = Ipt; 
            continue; 
        } else {

            num_land++;

            // Check if we've already mapped this point
            if (data.pt_maps_to[Ipt] < Npts) { continue; }

            // Begin a depth-first search through the adjacency
            // matrix, mapping all land-only neighbours to this one
            // This point will become the representative for this 'island'
            points_to_test.clear();
            data.pt_maps_to[Ipt] = Ipt; 
            num_mapped_onto++;


            if (constants::GRID_TYPE == constants::GridType::MeshGrid) {
                // Get coordinates
                Index1to4( Ipt, Itime, Idepth, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon );

                // Add neighbours to to-test list
                Ineigh = ( (Ilon+1) % Nlon + Nlon ) % Nlon;
                points_to_test.push_back( Index(0, 0, Ilat, Ineigh, 1, 1, Nlat, Nlon) );

                Ineigh = ( (Ilon-1) % Nlon + Nlon ) % Nlon;
                points_to_test.push_back( Index(0, 0, Ilat, Ineigh, 1, 1, Nlat, Nlon) );

                if (Ilat > 0) {
                    points_to_test.push_back( Index(0, 0, Ilat-1, Ilon, 1, 1, Nlat, Nlon) );
                }
                if (Ilat < Nlat-1) {
                    points_to_test.push_back( Index(0, 0, Ilat+1, Ilon, 1, 1, Nlat, Nlon) );
                }
            } else {
                for (size_t neighbour_ind = 0; neighbour_ind < data.num_neighbours; neighbour_ind++ ) {
                    points_to_test.push_back( data.adjacency_indices.at(Ipt)[neighbour_ind] );
                }
            }

            // So long as we still have points to test, keep testing!
            while ( points_to_test.size() > 0 ) {
                Itest = points_to_test.front();
                points_to_test.pop_front();
                if ( data.pt_maps_to[Itest] < Npts ) { continue; } // already mapped, skip
                else if (data.mask[Itest]) { continue; } // water, so skip
                else {
                    // Otherwise, map it to Ipt, and add its neighbours to the test list
                    data.pt_maps_to[Itest] = Ipt;
                    num_mapped_points++;

                    if (constants::GRID_TYPE == constants::GridType::MeshGrid) {
                        // Get coordinates
                        Index1to4( Itest, Itime, Idepth, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon );

                        Ineigh = ( (Ilon+1) % Nlon + Nlon ) % Nlon;
                        Ineighbour = Index(0, 0, Ilat, Ineigh, 1, 1, Nlat, Nlon);
                        if ( data.pt_maps_to[Ineighbour] == Npts ) {
                            // Only add the neighbour if we haven't already looked at it
                            points_to_test.push_back( Ineighbour );
                        }

                        Ineigh = ( (Ilon-1) % Nlon + Nlon ) % Nlon;
                        Ineighbour = Index(0, 0, Ilat, Ineigh, 1, 1, Nlat, Nlon);
                        if ( data.pt_maps_to[Ineighbour] == Npts ) {
                            // Only add the neighbour if we haven't already looked at it
                            points_to_test.push_back( Ineighbour );
                        }

                        if (Ilat > 0) {
                            Ineighbour = Index(0, 0, Ilat-1, Ilon, 1, 1, Nlat, Nlon);
                            if ( data.pt_maps_to[Ineighbour] == Npts ) {
                                // Only add the neighbour if we haven't already looked at it
                                points_to_test.push_back( Ineighbour );
                            }
                        }
                        if (Ilat < Nlat-1) {
                            Ineighbour = Index(0, 0, Ilat+1, Ilon, 1, 1, Nlat, Nlon);
                            if ( data.pt_maps_to[Ineighbour] == Npts ) {
                                // Only add the neighbour if we haven't already looked at it
                                points_to_test.push_back( Ineighbour );
                            }
                        }
                    } else {
                        for (size_t neighbour_ind = 0; neighbour_ind < data.num_neighbours; neighbour_ind++ ) {
                            size_t Jpt = data.adjacency_indices.at(Itest)[neighbour_ind];
                            if ( data.pt_maps_to[Jpt] == Npts ) {
                                // Only add the neighbour if we haven't already looked at it
                                points_to_test.push_back( data.adjacency_indices.at(Itest)[neighbour_ind] );
                            }
                        }
                    }
                }
            }
        }
    }
    #if DEBUG >= 0
    fprintf( stdout, "Mapping %'zu land points onto %'zu 'islands' (contiguous land masses).\n",
            num_mapped_points+num_mapped_onto, num_mapped_onto );
    #endif
}
