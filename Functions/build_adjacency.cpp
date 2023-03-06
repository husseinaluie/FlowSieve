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
    const size_t num_pts = longitude.size(),
                 num_neighbours = 8;

    adjacency_indices.resize( num_pts );
    adjacency_projected_x.resize( num_pts );
    adjacency_projected_y.resize( num_pts );
    adjacency_distances.resize( num_pts );

    adjacency_ddlon_weights.resize( num_pts );
    adjacency_ddlat_weights.resize( num_pts );

    // Be a litlle sneaky. We want to make this search faster if possible.
    // To do that, it would be helpful to have a coarse estimate of how far apart points
    // are. 
    // Note: surface area of a sphere is 4 * pi * R^2
    //       surface area / number of points ~ area per grid point
    //       sqrt( area per grid point ) ~ side length of grid, i.e. spacing between points
    // Is this perfect? Nope. But it'll get us in the ballpark
    const double    typical_spacing = sqrt( 4 * M_PI / num_pts ) * constants::R_earth,
                    max_distance = 10 * typical_spacing;

    std::vector<double> neighbour_x, neighbour_y, neighbour_dist;
    std::vector<size_t> neighbour_ind;
    std::bitset<16> stencil_set;

    double pt_lon, pt_lat, search_lon, search_lat, poleward_lat, del_lon_lim,
           proj_x, proj_y, proj_theta, local_dist, denom, denom_sign, stencil_count,
           stencil_min_x, stencil_min_y, stencil_max_x, stencil_max_y;
    bool near_pole;
    int Istencil, II, JJ, angle_ind;

    std::vector< double > xi, yi, ddx_coeffs, ddy_coeffs;

    #pragma omp parallel \
    default(none) \
    shared( longitude, latitude, adjacency_indices, adjacency_projected_x, adjacency_projected_y, \
            adjacency_distances, adjacency_ddlon_weights, adjacency_ddlat_weights ) \
    private( pt_index, search_index, pt_lon, pt_lat, search_lon, search_lat, angle_ind, \
             proj_x, proj_y, proj_theta, near_pole, poleward_lat, del_lon_lim, local_dist, \
             Ineighbour, neighbour_x, neighbour_y, neighbour_dist, neighbour_ind, stencil_set, \
             Istencil, II, JJ, xi, yi, denom, denom_sign, stencil_count, ddx_coeffs, ddy_coeffs, \
             stencil_min_x, stencil_min_y, stencil_max_x, stencil_max_y ) \
    firstprivate( num_pts, max_distance, stdout )
    {
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

                // Next, a slightly less laxy pruning by longitude
                poleward_lat = fabs(pt_lat) + (max_distance / constants::R_earth);
                near_pole = poleward_lat >= (M_PI * 89.8 / 180);
                del_lon_lim = std::fmin( M_PI, max_distance / ( cos(poleward_lat) * constants::R_earth ) );
                if (not( (near_pole) or ( cos( pt_lon - search_lon ) > cos(del_lon_lim) ) ) ) {
                    continue;
                }

                near_sided_project( proj_x, proj_y, search_lon, search_lat, pt_lon, pt_lat, 2e5 );
                local_dist = distance( search_lon, search_lat, pt_lon, pt_lat );
                if (local_dist < 1e-0) { continue; } // discard points that are erroneously close
                                                     // I don't really know why this happens
                                                     // I think it's rounding errors in the 
                                                     // map projection operator?

                // We're dividing angles up into 360 / num_neighbours = 45 degree segments
                proj_theta = std::fmax( 0, std::fmin( 360-1e-10, (M_PI + atan2( proj_y, proj_x )) * 180. / M_PI ));
                angle_ind = proj_theta / (int)( 360. / num_neighbours );
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

            }

            // Now that we've searched all of the points to the adjacent ones, 
            // store the adjacent values and move on.
            adjacency_indices.at(    pt_index).resize(num_neighbours, 0);
            adjacency_projected_x.at(pt_index).resize(num_neighbours, 0.);
            adjacency_projected_y.at(pt_index).resize(num_neighbours, 0.);
            adjacency_distances.at(  pt_index).resize(num_neighbours, 0.);

            for (Ineighbour = 0; Ineighbour < num_neighbours; Ineighbour++) {
                //fprintf( stdout, "%'zu: %'zu, %g, %g, %g\n", 
                //        pt_index, neighbour_ind[Ineighbour], neighbour_x[Ineighbour], neighbour_y[Ineighbour],
                //        neighbour_dist[Ineighbour] );
                adjacency_indices.at(pt_index)[Ineighbour] = neighbour_ind[Ineighbour];
                adjacency_projected_x.at(pt_index)[Ineighbour] = neighbour_x[Ineighbour];
                adjacency_projected_y.at(pt_index)[Ineighbour] = neighbour_y[Ineighbour];
                adjacency_distances.at(pt_index)[Ineighbour] = neighbour_dist[Ineighbour];
            }

            // Since we're here, let's also build the differentiation weights, since we'll be using
            // the adjacency points as our differentiation stencil.
            // 8 neighbours + centre point = 9 point stencil
            // But a bi-linear spline only requires 4 points
            // So we'll just take the average over the different available
            // differentiation stencils
            // Some will be ill-conditioned, so ignore those ones
            xi.resize(4);
            yi.resize(4);
            adjacency_ddlon_weights.at(pt_index).resize(num_neighbours+1, 0.);
            adjacency_ddlat_weights.at(pt_index).resize(num_neighbours+1, 0.);
            stencil_count = 0;

            for ( Istencil = 0; Istencil < pow(2,num_neighbours+1); Istencil++ ) {

                // First, reject any stencil subsets that don't contain exactly four points
                stencil_set = std::bitset<16>(Istencil);
                if (stencil_set.count() != 4) { continue; }

                // Pull out the four points
                // Get the xi and yi for the stencil
                // the 'last neighbour' is the point itself
                II = 0;
                stencil_min_x = 1e100;
                stencil_min_y = 1e100;
                stencil_max_x = -1e100;
                stencil_max_y = -1e100;
                for ( JJ = 0; JJ < num_neighbours+1; JJ++ ) {
                    // 
                    if ( stencil_set[JJ] == 0 ) { continue; }
                    xi[II] = (JJ == num_neighbours) ? 0. : adjacency_projected_x.at(pt_index)[JJ];
                    yi[II] = (JJ == num_neighbours) ? 0. : adjacency_projected_y.at(pt_index)[JJ];
                    II++;

                    stencil_min_x = fmin( stencil_min_x, xi[II] );
                    stencil_max_x = fmax( stencil_max_x, xi[II] );

                    stencil_min_y = fmin( stencil_min_y, yi[II] );
                    stencil_max_y = fmax( stencil_max_y, yi[II] );
                }
                if ( stencil_min_x * stencil_max_x > 1e-10 ) { continue; }
                if ( stencil_min_y * stencil_max_y > 1e-10 ) { continue; }

                // Denominator / determinant of the spline matrix
                denom = - xi[0] * (   xi[1] * (yi[0] - yi[1]) * (yi[2] - yi[3]) 
                                    - xi[2] * (yi[0] - yi[2]) * (yi[1] - yi[3]) 
                                    + xi[3] * (yi[0] - yi[3]) * (yi[1] - yi[2])
                                  )
                        - xi[1] * xi[2] * (yi[0] - yi[3]) * (yi[1] - yi[2]) 
                        + xi[1] * xi[3] * (yi[0] - yi[2]) * (yi[1] - yi[3]) 
                        - xi[2] * xi[3] * (yi[0] - yi[1]) * (yi[2] - yi[3]);
               
                denom_sign = (denom > 0) ? 1 : (denom < 0) ? -1 : 0;

                // If the denominator is too small in magnitude, discard this stencil.
                if ( fabs(denom) < 1e-0 ) { 
                    //fprintf(stdout, "Pt %'zu Rejecting stencil %d (denom = %g) (%.6e, %.6e, %.6e, %.6e) (%.6e, %.6e, %.6e, %.6e)\n", 
                    //        pt_index, Istencil, denom, xi[0], xi[1], xi[2], xi[3], yi[0], yi[1], yi[2], yi[3] ); 
                    continue; 
                }
                //stencil_count++;
                stencil_count += fabs(denom);

                ddx_coeffs.clear();
                // x1 y1 y2 - x2 y1 y2 - x1 y1 y3 + x3 y1 y3 + x2 y2 y3 -  x3 y2 y3,
                ddx_coeffs.push_back( (
                            + xi[1] * yi[1] * yi[2] 
                            - xi[2] * yi[1] * yi[2]
                            - xi[1] * yi[1] * yi[3]
                            + xi[3] * yi[1] * yi[3]
                            + xi[2] * yi[2] * yi[3]
                            - xi[3] * yi[2] * yi[3]
                            ) * denom_sign );
                            //) / denom );

                // -x0 y0 y2 + x2 y0 y2 + x0 y0 y3 - x3 y0 y3 - x2 y2 y3 +  x3 y2 y3,
                ddx_coeffs.push_back( (
                            - xi[0] * yi[0] * yi[2]
                            + xi[2] * yi[0] * yi[2]
                            + xi[0] * yi[0] * yi[3]
                            - xi[3] * yi[0] * yi[3]
                            - xi[2] * yi[2] * yi[3]
                            + xi[3] * yi[2] * yi[3]
                            ) * denom_sign );
                            //) / denom );

                //  x0 y0 y1 - x1 y0 y1 - x0 y0 y3 + x3 y0 y3 + x1 y1 y3 -  x3 y1 y3,
                ddx_coeffs.push_back( (
                            + xi[0] * yi[0] * yi[1]
                            - xi[1] * yi[0] * yi[1]
                            - xi[0] * yi[0] * yi[3]
                            + xi[3] * yi[0] * yi[3]
                            + xi[1] * yi[1] * yi[3]
                            - xi[3] * yi[1] * yi[3]
                            ) * denom_sign );
                            //) / denom );

                // -x0 y0 y1 + x1 y0 y1 + x0 y0 y2 - x2 y0 y2 - x1 y1 y2 +  x2 y1 y2
                ddx_coeffs.push_back( (
                            - xi[0] * yi[0] * yi[1]
                            + xi[1] * yi[0] * yi[1]
                            + xi[0] * yi[0] * yi[2]
                            - xi[2] * yi[0] * yi[2]
                            - xi[1] * yi[1] * yi[2]
                            + xi[2] * yi[1] * yi[2]
                            ) * denom_sign );
                            //) / denom );

                II = 0;
                for ( JJ = 0; JJ < num_neighbours+1; JJ++ ) {
                    if ( stencil_set[JJ] == 0 ) { continue; }
                    adjacency_ddlon_weights.at(pt_index)[JJ] += ddx_coeffs[II];
                    II++;
                }


                ddy_coeffs.clear();
                // -x1 x2 y1 + x1 x3 y1 + x1 x2 y2 - x2 x3 y2 - x1 x3 y3 + x2 x3 y3,
                // x1 x2 (-y1 + y2) + x1 x3 (y1 - y3) + x2 x3 (-y2 + y3),
                ddy_coeffs.push_back( (
                              xi[1] * xi[2] * (-yi[1] + yi[2])
                            + xi[1] * xi[3] * ( yi[1] - yi[3])
                            + xi[2] * xi[3] * (-yi[2] + yi[3])
                            ) * denom_sign );
                            //) / denom );

                // x0 x2 y0 - x0 x3 y0 - x0 x2 y2 + x2 x3 y2 + x0 x3 y3 -  x2 x3 y3, 
                // x0 x2 ( y0 - y2) + x2 x3 (y2 - y3) + x0 x3 (-y0 + y3),
                ddy_coeffs.push_back( (
                            + xi[0] * xi[2] * ( yi[0] - yi[2])
                            + xi[2] * xi[3] * ( yi[2] - yi[3])
                            + xi[0] * xi[3] * (-yi[0] + yi[3])
                            ) * denom_sign );
                            //) / denom );

                // -x0 x1 y0 + x0 x3 y0 + x0 x1 y1 - x1 x3 y1 - x0 x3 y3 +  x1 x3 y3,
                // x0 x1 (-y0 + y1) + x0 x3 (y0 - y3) + x1 x3 (-y1 + y3),
                ddy_coeffs.push_back( (
                            - xi[0] * xi[1] * yi[0]
                            + xi[0] * xi[3] * yi[0]
                            + xi[0] * xi[1] * yi[1]
                            - xi[1] * xi[3] * yi[1]
                            - xi[0] * xi[3] * yi[3]
                            + xi[1] * xi[3] * yi[3]
                            ) * denom_sign );
                            //) / denom );

                // x0 x1 y0 - x0 x2 y0 - x0 x1 y1 + x1 x2 y1 + x0 x2 y2 - x1 x2 y2
                // x0 x1 ( y0 - y1) + x1 x2 (y1 - y2) + x0 x2 (-y0 + y2)
                ddy_coeffs.push_back( (
                            + xi[0] * xi[1] * yi[0]
                            - xi[0] * xi[2] * yi[0]
                            - xi[0] * xi[1] * yi[1]
                            + xi[1] * xi[2] * yi[1]
                            + xi[0] * xi[2] * yi[2]
                            - xi[1] * xi[2] * yi[2]
                            ) * denom_sign );
                            //) / denom );

                II = 0;
                for ( JJ = 0; JJ < num_neighbours+1; JJ++ ) {
                    if ( stencil_set[JJ] == 0 ) { continue; }
                    adjacency_ddlat_weights.at(pt_index)[JJ] += ddy_coeffs[II];
                    II++;
                }
            }
            assert( stencil_count > 0 );
            for ( JJ = 0; JJ < num_neighbours+1; JJ++ ) {
                // factor R*cos(lat) comes from converting y deriv to lat deriv
                adjacency_ddlon_weights.at(pt_index)[JJ] *= 
                    constants::R_earth * cos(pt_lat) / stencil_count;
                // factor R comes from converting y deriv to lat deriv
                adjacency_ddlat_weights.at(pt_index)[JJ] *= constants::R_earth / stencil_count;
            }
        }
    }
}
