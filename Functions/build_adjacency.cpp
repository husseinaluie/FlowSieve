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
#include <bitset>

#include "../ALGLIB/stdafx.h"
#include "../ALGLIB/interpolation.h"

#include "../functions.hpp"
#include "../constants.hpp"

// Method of the dataset class
// i.e. source_data.build_adjacency();
void dataset::build_adjacency(
        const MPI_Comm comm
        ) {

    #if DEBUG >= 1
    int wRank, wSize;
    MPI_Comm_rank( comm, &wRank );
    MPI_Comm_size( comm, &wSize );

    if (wRank == 0) {
        fprintf( stdout, "Building adjacency matrix.\n" );
        fflush(stdout);
    }
    #endif

    size_t pt_index, search_index, Ineighbour;
    const size_t num_pts = longitude.size();

    adjacency_projected_x.resize( num_pts );
    adjacency_projected_y.resize( num_pts );
    adjacency_distances.resize( num_pts );

    adjacency_ddlon_weights.resize( num_pts );
    adjacency_ddlat_weights.resize( num_pts );

    adjacency_d2dlon2_weights.resize( num_pts );
    adjacency_d2dlat2_weights.resize( num_pts );


    std::vector<double> neighbour_x, neighbour_y, neighbour_dist;
    std::vector<size_t> neighbour_ind;
    std::bitset<16> stencil_set;

    double pt_lon, pt_lat, search_lon, search_lat, poleward_lat, del_lon_lim,
           proj_x, proj_y, proj_theta, local_dist, denom, denom_sign, stencil_count,
           stencil_min_x, stencil_min_y, stencil_max_x, stencil_max_y;
    bool near_pole, similar_angle;
    int Istencil, II, JJ, angle_ind;


    std::vector< double > xi, yi, ddx_coeffs, ddy_coeffs;

    std::vector< double > LHS_vec;
    alglib::real_2d_array LHS;
    alglib::matinvreport report;
    double *LHS_list;
    alglib::ae_int_t alglib_info;
    bool solve;

    #pragma omp parallel \
    default(none) \
    shared( longitude, latitude, adjacency_indices, adjacency_projected_x, adjacency_projected_y, \
            adjacency_distances, adjacency_ddlon_weights, adjacency_ddlat_weights ) \
    private( pt_index, pt_lon, pt_lat, search_lon, search_lat, \
             proj_x, proj_y, proj_theta, near_pole, poleward_lat, del_lon_lim, local_dist, \
             Ineighbour, neighbour_x, neighbour_y, neighbour_dist, neighbour_ind, stencil_set, \
             Istencil, II, JJ, xi, yi, denom, denom_sign, stencil_count, ddx_coeffs, ddy_coeffs, \
             stencil_min_x, stencil_min_y, stencil_max_x, stencil_max_y, \
             LHS, report, alglib_info, LHS_vec, solve, LHS_list, search_index \
           ) \
    firstprivate( num_pts, stdout, alglib::xdefault )
    {
        LHS_vec.resize( (num_neighbours+1)*(num_neighbours+1), 0. );
        LHS.setlength( num_neighbours+1, num_neighbours+1 );
        LHS.attach_to_ptr( num_neighbours+1, num_neighbours+1, &LHS_vec[0] );

        #pragma omp for collapse(1) schedule(guided)
        for ( pt_index = 0; pt_index < num_pts; pt_index++ ) {

            neighbour_x.clear();
            neighbour_y.clear();
            neighbour_dist.clear();
            neighbour_ind.clear();

            neighbour_x.resize( num_neighbours, 0. );
            neighbour_y.resize( num_neighbours, 0. );
            neighbour_dist.resize( num_neighbours, 0 );
            neighbour_ind.resize( num_neighbours, 0 );

            #if DEBUG > 0
            pt_lon = longitude.at( pt_index );
            pt_lat = latitude.at(  pt_index );
            #else
            pt_lon = longitude[ pt_index ];
            pt_lat = latitude[  pt_index ];
            #endif

            for ( Ineighbour = 0; Ineighbour < num_neighbours+1; Ineighbour++ ) {

                search_index = adjacency_indices.at(pt_index)[Ineighbour];

                #if DEBUG > 0
                search_lon = longitude.at( search_index );
                search_lat = latitude.at(  search_index );
                #else
                search_lon = longitude[ search_index ];
                search_lat = latitude[  search_index ];
                #endif

                local_dist = distance( search_lon, search_lat, pt_lon, pt_lat );
                near_sided_project( proj_x, proj_y, search_lon, search_lat, pt_lon, pt_lat, 2e5 );

                neighbour_x[   Ineighbour] = proj_x;
                neighbour_y[   Ineighbour] = proj_y;
                neighbour_ind[ Ineighbour] = search_index;
                neighbour_dist[Ineighbour] = local_dist;
            }

            // Now that we've searched all of the points to the adjacent ones, 
            // store the adjacent values and move on.
            //adjacency_indices.at(    pt_index).resize(num_neighbours+1, 0);
            adjacency_projected_x.at(pt_index).resize(num_neighbours+1, 0.);
            adjacency_projected_y.at(pt_index).resize(num_neighbours+1, 0.);
            adjacency_distances.at(  pt_index).resize(num_neighbours+1, 0.);

            // the last 'neighbour' is the point itself
            //adjacency_indices.at(    pt_index)[num_neighbours] = pt_index;
            adjacency_projected_x.at(pt_index)[num_neighbours] = 0.;
            adjacency_projected_y.at(pt_index)[num_neighbours] = 0.;
            adjacency_distances.at(  pt_index)[num_neighbours] = 0.;

            LHS_vec.resize( (num_neighbours+1) * (num_neighbours+1), 0. );
            for (Ineighbour = 0; Ineighbour < num_neighbours+1; Ineighbour++) {

                if (Ineighbour < num_neighbours) {
                    //adjacency_indices.at(pt_index)[Ineighbour] = neighbour_ind[Ineighbour];
                    adjacency_projected_x.at(pt_index)[Ineighbour] = neighbour_x[Ineighbour];
                    adjacency_projected_y.at(pt_index)[Ineighbour] = neighbour_y[Ineighbour];
                    adjacency_distances.at(pt_index)[Ineighbour] = neighbour_dist[Ineighbour];
                }
            }
            /*
            for (Ineighbour = 0; Ineighbour < num_neighbours+1; Ineighbour++) {

                // 1
                LHS_vec.at( Ineighbour * (num_neighbours+1) + 0 ) = 1.;

                // x
                LHS_vec.at( Ineighbour * (num_neighbours+1) + 1 ) = 
                    neighbour_x[Ineighbour];

                // y
                LHS_vec.at( Ineighbour * (num_neighbours+1) + 2 ) = 
                    neighbour_y[Ineighbour];

                // x y 
                LHS_vec.at( Ineighbour * (num_neighbours+1) + 3 ) = 
                    neighbour_x[Ineighbour] * neighbour_y[Ineighbour];

                // x^2
                if (num_neighbours >= 4) {
                    LHS_vec.at( Ineighbour * (num_neighbours+1) + 4 ) = 
                        neighbour_x[Ineighbour] * neighbour_x[Ineighbour];
                }
                
                // y^2
                if (num_neighbours >= 5) {
                    LHS_vec.at( Ineighbour * (num_neighbours+1) + 5 ) = 
                        neighbour_y[Ineighbour] * neighbour_y[Ineighbour];
                }
                
                // x * y^2
                if (num_neighbours >= 6) {
                    LHS_vec.at( Ineighbour * (num_neighbours+1) + 6 ) = 
                        neighbour_x[Ineighbour] * neighbour_y[Ineighbour] * neighbour_y[Ineighbour];
                }
                
                // x^2 * y
                if (num_neighbours >= 7) {
                    LHS_vec.at( Ineighbour * (num_neighbours+1) + 7 ) = 
                        neighbour_x[Ineighbour] * neighbour_x[Ineighbour] * neighbour_y[Ineighbour];
                }
                
                // x^2 * y^2
                if (num_neighbours >= 8) {
                    LHS_vec.at( Ineighbour * (num_neighbours+1) + 8 ) = 
                        neighbour_x[Ineighbour] * neighbour_x[Ineighbour] * neighbour_y[Ineighbour] * neighbour_y[Ineighbour];
                }

            }
            */
            for (Ineighbour = 0; Ineighbour < num_neighbours+1; Ineighbour++) {
                double dx = neighbour_x[Ineighbour], dy = neighbour_y[Ineighbour];

                LHS( Ineighbour, 0 ) = 1;
                LHS( Ineighbour, 1 ) = dx;
                LHS( Ineighbour, 2 ) = dy;
                LHS( Ineighbour, 3 ) = (1./2) * dx*dx;
                LHS( Ineighbour, 4 ) = (1./2) * 2. * dx*dy;
                LHS( Ineighbour, 5 ) = (1./2) * dy*dy;
                LHS( Ineighbour, 6 ) = (1./6) * dx*dx*dx;
                LHS( Ineighbour, 7 ) = (1./6) * 3. * dx*dx*dy;
                LHS( Ineighbour, 8 ) = (1./6) * dy*dy*dy;
            }


            alglib::rmatrixinverse( LHS, alglib_info, report );

            adjacency_ddlon_weights.at(pt_index).resize(num_neighbours+1, 0.);
            adjacency_ddlat_weights.at(pt_index).resize(num_neighbours+1, 0.);

            if (num_neighbours >= 4) {
                adjacency_d2dlon2_weights.at(pt_index).resize(num_neighbours+1, 0.);
                adjacency_d2dlat2_weights.at(pt_index).resize(num_neighbours+1, 0.);
            }
            for ( JJ = 0; JJ < num_neighbours+1; JJ++ ) {

                // factor R*cos(lat) comes from converting proj-x deriv to lon deriv
                adjacency_ddlon_weights.at(pt_index)[JJ] = 
                    LHS( 1, JJ ) * constants::R_earth * cos(pt_lat);

                // factor R comes from converting proj-y deriv to lat deriv
                adjacency_ddlat_weights.at(pt_index)[JJ] = 
                    LHS( 2, JJ ) * constants::R_earth;

                if (num_neighbours >= 4) {
                    // second lon derivative
                    //adjacency_d2dlon2_weights.at(pt_index)[JJ] = 2 * LHS( 4, JJ ) * 
                    //    pow( constants::R_earth * cos(pt_lat), 2);
                    adjacency_d2dlon2_weights.at(pt_index)[JJ] = 
                          LHS( 3, JJ ) * pow( constants::R_earth * cos(pt_lat), 2)
                        + LHS( 2, JJ ) * constants::R_earth * sin(pt_lat) * cos(pt_lat);

                    // second lat derivative
                    //adjacency_d2dlat2_weights.at(pt_index)[JJ] = 2 * LHS( 5, JJ ) * 
                    //    pow( constants::R_earth, 2);
                    adjacency_d2dlat2_weights.at(pt_index)[JJ] = 
                        LHS( 5, JJ ) * pow( constants::R_earth, 2);
                }

            }
        }
    }
}
