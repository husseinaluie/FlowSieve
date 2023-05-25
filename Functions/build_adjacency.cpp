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

    adjacency_indices.resize( num_pts );
    adjacency_projected_x.resize( num_pts );
    adjacency_projected_y.resize( num_pts );
    adjacency_distances.resize( num_pts );

    adjacency_ddlon_weights.resize( num_pts );
    adjacency_ddlat_weights.resize( num_pts );

    adjacency_d2dlon2_weights.resize( num_pts );
    adjacency_d2dlat2_weights.resize( num_pts );

    // Be a litlle sneaky. We want to make this search faster if possible.
    // To do that, it would be helpful to have a coarse estimate of how far apart points
    // are. 
    // Note: surface area of a sphere is 4 * pi * R^2
    //       surface area / number of points ~ area per grid point
    //       sqrt( area per grid point ) ~ side length of grid, i.e. spacing between points
    // Is this perfect? Nope. But it'll get us in the ballpark
    const double    typical_spacing = sqrt( 4 * M_PI / num_pts ) * constants::R_earth,
                    max_distance = 10 * typical_spacing,
                    angle_per_neighbour = 360. / num_neighbours;

    std::vector<double> neighbour_x, neighbour_y, neighbour_dist;
    std::vector<size_t> neighbour_ind;
    std::bitset<16> stencil_set;

    double pt_lon, pt_lat, search_lon, search_lat, poleward_lat, del_lon_lim,
           proj_x, proj_y, proj_theta, local_dist, denom, denom_sign, stencil_count,
           stencil_min_x, stencil_min_y, stencil_max_x, stencil_max_y;
    bool near_pole, similar_angle;
    int Istencil, II, JJ, angle_ind;

    double max_neighbour_dist = -1, ref_angle, curr_angle;
    size_t furthest_neighbour_ind = 0 ;

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
    private( pt_index, search_index, pt_lon, pt_lat, search_lon, search_lat, angle_ind, \
             proj_x, proj_y, proj_theta, near_pole, poleward_lat, del_lon_lim, local_dist, \
             Ineighbour, neighbour_x, neighbour_y, neighbour_dist, neighbour_ind, stencil_set, \
             Istencil, II, JJ, xi, yi, denom, denom_sign, stencil_count, ddx_coeffs, ddy_coeffs, \
             stencil_min_x, stencil_min_y, stencil_max_x, stencil_max_y, \
             LHS, report, alglib_info, LHS_vec, solve, LHS_list, \
             max_neighbour_dist, furthest_neighbour_ind, similar_angle, ref_angle, curr_angle \
           ) \
    firstprivate( num_pts, max_distance, stdout, alglib::xdefault, angle_per_neighbour )
    {
        LHS_vec.resize( (num_neighbours+1)*(num_neighbours+1), 0. );
        LHS.setlength( num_neighbours+1, num_neighbours+1 );
        LHS.attach_to_ptr( num_neighbours+1, num_neighbours+1, &LHS_vec[0] );

        #pragma omp for collapse(1) schedule(guided)
        for ( pt_index = 0; pt_index < num_pts; pt_index++ ) {

            //fprintf( stdout, "Starting point %'zu\n", pt_index );

            search_index = 0;

            neighbour_x.clear();
            neighbour_y.clear();
            neighbour_dist.clear();
            neighbour_ind.clear();

            neighbour_x.resize( num_neighbours, 0. );
            neighbour_y.resize( num_neighbours, 0. );
            neighbour_dist.resize( num_neighbours, 1e100 );
            neighbour_ind.resize( num_neighbours, 0 );

            max_neighbour_dist = 1e100;
            furthest_neighbour_ind = 0;

            #if DEBUG > 0
            pt_lon = longitude.at( pt_index );
            pt_lat = latitude.at(  pt_index );
            #else
            pt_lon = longitude[ pt_index ];
            pt_lat = latitude[  pt_index ];
            #endif

            for ( search_index = 0; search_index < num_pts; search_index++ ) {

                if ( pt_index == search_index ) { continue; }

                #if DEBUG > 0
                search_lon = longitude.at( search_index );
                search_lat = latitude.at(  search_index );
                #else
                search_lon = longitude[ search_index ];
                search_lat = latitude[  search_index ];
                #endif

                // First, a lazy pruning by latitude to avoid excessive distance calculations
                if ( fabs( pt_lat - search_lat ) >= (max_distance / constants::R_earth) ) {
                    continue;
                }

                // Next, a slightly less lazy pruning by longitude
                poleward_lat = fabs(pt_lat) + (max_distance / constants::R_earth);
                near_pole = poleward_lat >= (M_PI * 89.8 / 180);
                del_lon_lim = std::fmin( M_PI, max_distance / ( cos(poleward_lat) * constants::R_earth ) );
                if (not( (near_pole) or ( cos( pt_lon - search_lon ) > cos(del_lon_lim) ) ) ) {
                    continue;
                }

                near_sided_project( proj_x, proj_y, search_lon, search_lat, pt_lon, pt_lat, 2e5 );
                local_dist = distance( search_lon, search_lat, pt_lon, pt_lat );
                if (local_dist < 1e-4) { continue; } // discard points that are erroneously close
                                                     // I don't really know why this happens
                                                     // I think it's rounding errors in the 
                                                     // map projection operator?

                /*
                // We're dividing angles up into 360 / num_neighbours = 45 degree segments
                proj_theta = std::fmax( 0, std::fmin( 360-1e-10, (M_PI + atan2( proj_y, proj_x )) * 180. / M_PI ));
                if      ( proj_theta <=  15. ) { angle_ind = 0; }
                else if ( proj_theta <=  75. ) { angle_ind = 1; }
                else if ( proj_theta <= 105. ) { angle_ind = 2; }
                else if ( proj_theta <= 165. ) { angle_ind = 3; }
                else if ( proj_theta <= 195. ) { angle_ind = 4; }
                else if ( proj_theta <= 255. ) { angle_ind = 5; }
                else if ( proj_theta <= 285. ) { angle_ind = 6; }
                else if ( proj_theta <= 345. ) { angle_ind = 7; }
                else                           { angle_ind = 0; }
                //proj_theta = (proj_theta - 0.5 * angle_per_neighbour);
                //proj_theta += (proj_theta < 0) ? 360. : 0;
                //angle_ind = (int) std::fmax( 0, std::fmin( num_neighbours-1, round( proj_theta / angle_per_neighbour ) ) );
                //angle_ind = proj_theta / (int)( angle_per_neighbour );
                #if DEBUG >= 1
                if ( angle_ind < 0 ) { fprintf(stdout, "%d\n", angle_ind); assert(false); }
                if ( angle_ind >= num_neighbours ) { fprintf(stdout, "%d\n", angle_ind); assert(false); }
                #endif

                if ( local_dist < neighbour_dist.at(angle_ind) ) {
                    neighbour_x[   angle_ind] = proj_x;
                    neighbour_y[   angle_ind] = proj_y;
                    neighbour_ind[ angle_ind] = search_index;
                    neighbour_dist[angle_ind] = local_dist;
                }
                */

                if ( local_dist < neighbour_dist[furthest_neighbour_ind] ) {

                    // Also check that there aren't other points with similar angles
                    similar_angle = false;
                    for (II = 0; II < num_neighbours; II++) {
                        if (II == furthest_neighbour_ind) { continue; }
                        ref_angle = atan2( neighbour_y[II], neighbour_x[II] );
                        curr_angle = atan2(proj_y, proj_x);
                        if ( ( fabs( curr_angle - ref_angle ) < 5.*M_PI/180. )
                             or ( fabs( curr_angle - ref_angle ) > 2*M_PI - 5.*M_PI/180. ) )
                                {
                            //similar_angle = true;
                            // there's another point at this 'angle', so we can only replace this one
                            // so call that the 'furthest' to be replaced
                            if ( local_dist < neighbour_dist[II] ) {
                                furthest_neighbour_ind = II;
                            } else {
                                similar_angle = true;
                            }
                        }
                    }

                    if ( not(similar_angle) ) {
                        neighbour_x[   furthest_neighbour_ind] = proj_x;
                        neighbour_y[   furthest_neighbour_ind] = proj_y;
                        neighbour_ind[ furthest_neighbour_ind] = search_index;
                        neighbour_dist[furthest_neighbour_ind] = local_dist;
                        max_neighbour_dist = neighbour_dist[0];
                        furthest_neighbour_ind = 0;
                        for (II = 1; II < num_neighbours; II++) {
                            if (neighbour_dist[II] > max_neighbour_dist) {
                                furthest_neighbour_ind = II;
                                max_neighbour_dist = neighbour_dist[II];
                            }
                        }
                    }
                }
            }

            // Now that we've searched all of the points to the adjacent ones, 
            // store the adjacent values and move on.
            adjacency_indices.at(    pt_index).resize(num_neighbours+1, 0);
            adjacency_projected_x.at(pt_index).resize(num_neighbours+1, 0.);
            adjacency_projected_y.at(pt_index).resize(num_neighbours+1, 0.);
            adjacency_distances.at(  pt_index).resize(num_neighbours+1, 0.);

            // the last 'neighbour' is the point itself
            adjacency_indices.at(    pt_index)[num_neighbours] = pt_index;
            adjacency_projected_x.at(pt_index)[num_neighbours] = 0.;
            adjacency_projected_y.at(pt_index)[num_neighbours] = 0.;
            adjacency_distances.at(  pt_index)[num_neighbours] = 0.;

            LHS_vec.resize( (num_neighbours+1) * (num_neighbours+1), 0. );
            LHS_vec.at(      num_neighbours    * (num_neighbours+1) + 0 ) = 1.;
            LHS_vec.at(      num_neighbours    * (num_neighbours+1) + 1 ) = 0.;
            LHS_vec.at(      num_neighbours    * (num_neighbours+1) + 2 ) = 0.;
            LHS_vec.at(      num_neighbours    * (num_neighbours+1) + 3 ) = 0.;
            LHS_vec.at(      num_neighbours    * (num_neighbours+1) + 4 ) = 0.;
            LHS_vec.at(      num_neighbours    * (num_neighbours+1) + 5 ) = 0.;
            if (num_neighbours >= 6) {
                LHS_vec.at(      num_neighbours    * (num_neighbours+1) + 6 ) = 0.;
            }
            if (num_neighbours >= 7) {
                LHS_vec.at(      num_neighbours    * (num_neighbours+1) + 7 ) = 0.;
            }
            if (num_neighbours >= 8) {
                LHS_vec.at(      num_neighbours    * (num_neighbours+1) + 8 ) = 0.;
            }
            for (Ineighbour = 0; Ineighbour < num_neighbours; Ineighbour++) {

                adjacency_indices.at(pt_index)[Ineighbour] = neighbour_ind[Ineighbour];
                adjacency_projected_x.at(pt_index)[Ineighbour] = neighbour_x[Ineighbour];
                adjacency_projected_y.at(pt_index)[Ineighbour] = neighbour_y[Ineighbour];
                adjacency_distances.at(pt_index)[Ineighbour] = neighbour_dist[Ineighbour];

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
                LHS_vec.at( Ineighbour * (num_neighbours+1) + 4 ) = 
                    neighbour_x[Ineighbour] * neighbour_x[Ineighbour];
                
                // y^2
                LHS_vec.at( Ineighbour * (num_neighbours+1) + 5 ) = 
                    neighbour_y[Ineighbour] * neighbour_y[Ineighbour];
                
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


            alglib::rmatrixinverse( LHS, alglib_info, report );

            adjacency_ddlon_weights.at(pt_index).resize(num_neighbours+1, 0.);
            adjacency_ddlat_weights.at(pt_index).resize(num_neighbours+1, 0.);

            adjacency_d2dlon2_weights.at(pt_index).resize(num_neighbours+1, 0.);
            adjacency_d2dlat2_weights.at(pt_index).resize(num_neighbours+1, 0.);
            for ( JJ = 0; JJ < num_neighbours+1; JJ++ ) {

                // factor R*cos(lat) comes from converting proj-x deriv to lat deriv
                adjacency_ddlon_weights.at(pt_index)[JJ] = LHS( 1, JJ ) * 
                    constants::R_earth * cos(pt_lat);

                // factor R comes from converting proj-y deriv to lat deriv
                adjacency_ddlat_weights.at(pt_index)[JJ] = LHS( 2, JJ ) *
                    constants::R_earth;

                // second lon derivative
                adjacency_d2dlon2_weights.at(pt_index)[JJ] = 0.5 * LHS( 4, JJ ) * 
                    pow( constants::R_earth * cos(pt_lat), 2);

                // second lat derivative
                adjacency_d2dlat2_weights.at(pt_index)[JJ] = 0.5 * LHS( 5, JJ ) * 
                    pow( constants::R_earth, 2);

            }
        }
    }
}
