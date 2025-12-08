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
    int Istencil, II, JJ, KK, angle_ind;


    std::vector< double > xi, yi, ddx_coeffs, ddy_coeffs;

    alglib::real_2d_array LHS, svd_U, svd_VT;
    alglib::real_1d_array svd_S;
    alglib::matinvreport report;
    double *LHS_list;
    alglib::ae_int_t alglib_info;
    bool solve;

    const int num_rows = 6; // number of terms used in SVD decomp
    //const int num_rows = 8; // number of terms used in SVD decomp

    #pragma omp parallel \
    default(none) \
    shared( longitude, latitude, adjacency_indices, adjacency_projected_x, adjacency_projected_y, \
            adjacency_distances, adjacency_ddlon_weights, adjacency_ddlat_weights ) \
    private( pt_index, pt_lon, pt_lat, search_lon, search_lat, \
             proj_x, proj_y, proj_theta, near_pole, poleward_lat, del_lon_lim, local_dist, \
             Ineighbour, neighbour_x, neighbour_y, neighbour_dist, neighbour_ind, stencil_set, \
             Istencil, II, JJ, KK, xi, yi, denom, denom_sign, stencil_count, ddx_coeffs, ddy_coeffs, \
             stencil_min_x, stencil_min_y, stencil_max_x, stencil_max_y, \
             LHS, report, alglib_info, solve, LHS_list, search_index, \
             svd_U, svd_S, svd_VT \
           ) \
    firstprivate( num_pts, stdout, alglib::xdefault, num_rows )
    {
        LHS.setlength( 9, num_rows );

        svd_U.setlength(9,9);
        svd_S.setlength(9);
        svd_VT.setlength(num_rows,num_rows);

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
                //near_sided_project( proj_x, proj_y, search_lon, search_lat, pt_lon, pt_lat, 2e5 );
                //xs = np.cos( lats ) * np.sin( lons - lon )
                //ys = ( np.cos( lat ) * np.sin( lats ) - np.sin( lat ) * np.cos( lats ) * np.cos( lons - lon ) )

                proj_x = cos( search_lat ) * sin( search_lon - pt_lon );
                proj_y = cos( pt_lat ) * sin( search_lat ) - sin( pt_lat ) * cos( search_lat ) * cos( search_lon - pt_lon );

                neighbour_x[   Ineighbour] = proj_x;
                neighbour_y[   Ineighbour] = proj_y;
                neighbour_ind[ Ineighbour] = search_index;
                neighbour_dist[Ineighbour] = local_dist;
            }

            // Now that we've searched all of the points to the adjacent ones, 
            // store the adjacent values and move on.
            adjacency_projected_x.at(pt_index).resize(num_neighbours+1, 0.);
            adjacency_projected_y.at(pt_index).resize(num_neighbours+1, 0.);
            adjacency_distances.at(  pt_index).resize(num_neighbours+1, 0.);

            // the last 'neighbour' is the point itself
            adjacency_projected_x.at(pt_index)[num_neighbours] = 0.;
            adjacency_projected_y.at(pt_index)[num_neighbours] = 0.;
            adjacency_distances.at(  pt_index)[num_neighbours] = 0.;

            for (Ineighbour = 0; Ineighbour < num_neighbours+1; Ineighbour++) {
                if (Ineighbour < num_neighbours) {
                    adjacency_projected_x.at(pt_index)[Ineighbour] = neighbour_x[Ineighbour];
                    adjacency_projected_y.at(pt_index)[Ineighbour] = neighbour_y[Ineighbour];
                    adjacency_distances.at(pt_index)[Ineighbour] = neighbour_dist[Ineighbour];
                }
            }

            for (Ineighbour = 0; Ineighbour < num_neighbours+1; Ineighbour++) {
                double dx = neighbour_x[Ineighbour], dy = neighbour_y[Ineighbour];

                LHS( Ineighbour, 0 ) = 1;
                LHS( Ineighbour, 1 ) = dx;
                LHS( Ineighbour, 2 ) = dy;

                LHS( Ineighbour, 3 ) = dx*dx / 2.;
                LHS( Ineighbour, 4 ) = dx*dy;
                LHS( Ineighbour, 5 ) = dy*dy / 2.;

                //LHS( Ineighbour, 6 ) = dx*dx*dy / 2.;
                //LHS( Ineighbour, 7 ) = dx*dy*dy / 2.;
            }

            alglib::rmatrixsvd( LHS, 9, num_rows, 2, 2, 2, svd_S, svd_U, svd_VT );
            // LHS = U * S * VT
            // LHS is m x n
            // U is m x m
            // S is m x n diagonal (just the diagonal entries are stored in svd_S)
            // V is n x n
            //
            // inv(LHS) = V * S^-1 * UT

            // First, scale the top six rows of UT by the elements of S^-1
            //  this is equivalent to scaling the first six columns of U
            for ( JJ = 0; JJ < 9; JJ++ ) {
                for ( II = 0; II < num_rows; II++ ) {
                    svd_U(JJ,II) = svd_U(JJ,II) / svd_S(II);
                }
            }


            adjacency_ddlon_weights.at(pt_index).resize(num_neighbours+1, 0.);
            adjacency_ddlat_weights.at(pt_index).resize(num_neighbours+1, 0.);

            for ( JJ = 0; JJ < num_neighbours+1; JJ++ ) {

                // lon deriv: factor cos(lat) comes from converting proj-x deriv to lon deriv
                adjacency_ddlon_weights.at(pt_index)[JJ] = 0;
                adjacency_ddlat_weights.at(pt_index)[JJ] = 0;
                for ( KK = 0; KK < num_rows; KK++ ) {
                    adjacency_ddlon_weights.at(pt_index)[JJ] += svd_VT(KK,1) * svd_U(JJ,KK) * cos(pt_lat);
                    adjacency_ddlat_weights.at(pt_index)[JJ] += svd_VT(KK,2) * svd_U(JJ,KK);
                }
            }

            // These are just hold-overs from an older time
            adjacency_d2dlon2_weights.at(pt_index).resize(num_neighbours+1, 0.);
            adjacency_d2dlat2_weights.at(pt_index).resize(num_neighbours+1, 0.);
        }
    }
}
