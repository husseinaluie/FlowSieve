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

    size_t pt_index, search_index;
    const size_t num_pts = longitude.size();

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

    double UR_x, UR_y, UR_dist, UL_x, UL_y, UL_dist, LR_x, LR_y, LR_dist, LL_x, LL_y, LL_dist;
    size_t UR_ind, UL_ind, LR_ind, LL_ind;

    double pt_lon, pt_lat, search_lon, search_lat, poleward_lat, del_lon_lim,
           proj_x, proj_y, local_dist, denom, stencil_count;
    bool near_pole;
    int Istencil, II, JJ;

    std::vector< double > xi, yi, ddx_coeffs, ddy_coeffs;

    #pragma omp parallel \
    default(none) \
    shared( longitude, latitude, adjacency_indices, adjacency_projected_x, adjacency_projected_y, \
            adjacency_distances, adjacency_ddlon_weights, adjacency_ddlat_weights ) \
    private( pt_index, search_index, pt_lon, pt_lat, search_lon, search_lat, \
             proj_x, proj_y, near_pole, poleward_lat, del_lon_lim, local_dist, \
             LL_ind, LL_x, LL_y, LL_dist, LR_ind, LR_x, LR_y, LR_dist, \
             UL_ind, UL_x, UL_y, UL_dist, UR_ind, UR_x, UR_y, UR_dist, \
             Istencil, II, JJ, xi, yi, denom, stencil_count, ddx_coeffs, ddy_coeffs ) \
    firstprivate( num_pts, max_distance )
    {
        #pragma omp for collapse(1) schedule(guided)
        for ( pt_index = 0; pt_index < num_pts; pt_index++ ) {

            search_index = 0;

            UR_dist = 1e100;
            UL_dist = 1e100;
            LR_dist = 1e100;
            LL_dist = 1e100;

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

                //if ( (proj_x >= 0) and (proj_y > 1e-10) and ( local_dist < UR_dist ) ) {
                if ( (proj_y > fabs(proj_x)) and ( local_dist < UR_dist ) ) {
                    UR_x = proj_x;
                    UR_y = proj_y;
                    UR_ind = search_index;
                    UR_dist = local_dist;
                //} else if ( (proj_x > 1e-10) and (proj_y <= 0) and ( local_dist < LR_dist ) ) {
                } else if ( (proj_x > fabs(proj_y)) and ( local_dist < LR_dist ) ) {
                    LR_x = proj_x;
                    LR_y = proj_y;
                    LR_ind = search_index;
                    LR_dist = local_dist;
                //} else if ( (proj_x < -1e-10) and (proj_y >= 0) and ( local_dist < UL_dist ) ) {
                } else if ( (proj_y < -fabs(proj_x)) and ( local_dist < UL_dist ) ) {
                    UL_x = proj_x;
                    UL_y = proj_y;
                    UL_ind = search_index;
                    UL_dist = local_dist;
                //} else if ( (proj_x <= 0) and (proj_y < -1e-10) and ( local_dist < LL_dist ) ) {
                } else if ( (proj_x < -fabs(proj_y)) and ( local_dist < LL_dist ) ) {
                    LL_x = proj_x;
                    LL_y = proj_y;
                    LL_ind = search_index;
                    LL_dist = local_dist;
                }
            }

            assert( (UR_dist <= max_distance) and 
                    (UL_dist <= max_distance) and
                    (LR_dist <= max_distance) and
                    (LL_dist <= max_distance) 
                  ); // otherwise failed to build adjacency

            // Now that we've searched all of the points to the adjacent ones, 
            // store the adjacent values and move on.
            adjacency_indices.at(pt_index).resize(4, 0.);
            adjacency_projected_x.at(pt_index).resize(4, 0.);
            adjacency_projected_y.at(pt_index).resize(4, 0.);
            adjacency_distances.at(pt_index).resize(4, 0.);

            adjacency_indices.at(pt_index)[0] = LL_ind;
            adjacency_indices.at(pt_index)[1] = LR_ind;
            adjacency_indices.at(pt_index)[2] = UL_ind;
            adjacency_indices.at(pt_index)[3] = UR_ind;

            adjacency_projected_x.at(pt_index)[0] = LL_x;
            adjacency_projected_x.at(pt_index)[1] = LR_x;
            adjacency_projected_x.at(pt_index)[2] = UL_x;
            adjacency_projected_x.at(pt_index)[3] = UR_x;

            adjacency_projected_y.at(pt_index)[0] = LL_y;
            adjacency_projected_y.at(pt_index)[1] = LR_y;
            adjacency_projected_y.at(pt_index)[2] = UL_y;
            adjacency_projected_y.at(pt_index)[3] = UR_y;

            adjacency_distances.at(pt_index)[0] = LL_dist;
            adjacency_distances.at(pt_index)[1] = LR_dist;
            adjacency_distances.at(pt_index)[2] = UL_dist;
            adjacency_distances.at(pt_index)[3] = UR_dist;


            // Since we're here, let's also build the differentiation weights, since we'll be using
            // the adjacency points as our differentiation stencil.
            // 4 neighbours + centre point = 5 point stencil
            // But a bi-linear spline only requires 4 points
            // So we'll just take the average over the five different
            // differentiation stencils
            xi.resize(4);
            yi.resize(4);
            adjacency_ddlon_weights.at(pt_index).resize(5, 0.);
            adjacency_ddlat_weights.at(pt_index).resize(5, 0.);
            stencil_count = 0;
            for ( Istencil = 0; Istencil < 5; Istencil++ ) {

                // Get the xi and yi for the stencil
                //   [ 4 -> centre point ]
                II = 0;
                for ( JJ = 0; JJ < 5; JJ++ ) {
                    // JJ indicates which point is *excluded* from the stencil
                    // with 4 -> centre point and 0-3 mapping to the neighbours
                    if (JJ == Istencil) { continue; }
                    xi[II] = (JJ == 4) ? 0. : adjacency_projected_x.at(pt_index)[JJ];
                    yi[II] = (JJ == 4) ? 0. : adjacency_projected_y.at(pt_index)[JJ];
                    II++;
                }

                // Denominator / determinant of the spline matrix
                denom = - xi[0] * (   xi[1] * (yi[0] - yi[1]) * (yi[2] - yi[3]) 
                                    - xi[2] * (yi[0] - yi[2]) * (yi[1] - yi[3]) 
                                    + xi[3] * (yi[0] - yi[3]) * (yi[1] - yi[2])
                                  )
                        - xi[1] * xi[2] * (yi[0] - yi[3]) * (yi[1] - yi[2]) 
                        + xi[1] * xi[3] * (yi[0] - yi[2]) * (yi[1] - yi[3]) 
                        - xi[2] * xi[3] * (yi[0] - yi[1]) * (yi[2] - yi[3]);

                // If the denominator is too small in magnitude, discard this stencil.
                if ( fabs(denom) < 1e-10 ) { continue; }
                stencil_count++;

                ddx_coeffs.clear();
                // x1 y1 y2 - x2 y1 y2 - x1 y1 y3 + x3 y1 y3 + x2 y2 y3 -  x3 y2 y3,
                ddx_coeffs.push_back( (
                            + xi[1] * yi[1] * yi[2] 
                            - xi[2] * yi[1] * yi[2]
                            - xi[1] * yi[1] * yi[3]
                            + xi[3] * yi[1] * yi[3]
                            + xi[2] * yi[2] * yi[3]
                            - xi[3] * yi[2] * yi[3]
                            ) / denom );

                // -x0 y0 y2 + x2 y0 y2 + x0 y0 y3 - x3 y0 y3 - x2 y2 y3 +  x3 y2 y3,
                ddx_coeffs.push_back( (
                            - xi[0] * yi[0] * yi[2]
                            + xi[2] * yi[0] * yi[2]
                            + xi[0] * yi[0] * yi[3]
                            - xi[3] * yi[0] * yi[3]
                            - xi[2] * yi[2] * yi[3]
                            + xi[3] * yi[2] * yi[3]
                            ) / denom );

                //  x0 y0 y1 - x1 y0 y1 - x0 y0 y3 + x3 y0 y3 + x1 y1 y3 -  x3 y1 y3,
                ddx_coeffs.push_back( (
                            + xi[0] * yi[0] * yi[1]
                            - xi[1] * yi[0] * yi[1]
                            - xi[0] * yi[0] * yi[3]
                            + xi[3] * yi[0] * yi[3]
                            + xi[1] * yi[1] * yi[3]
                            - xi[3] * yi[1] * yi[3]
                            ) / denom );

                // -x0 y0 y1 + x1 y0 y1 + x0 y0 y2 - x2 y0 y2 - x1 y1 y2 +  x2 y1 y2
                ddx_coeffs.push_back( (
                            - xi[0] * yi[0] * yi[1]
                            + xi[1] * yi[0] * yi[1]
                            + xi[0] * yi[0] * yi[2]
                            - xi[2] * yi[0] * yi[2]
                            - xi[1] * yi[1] * yi[2]
                            + xi[2] * yi[1] * yi[2]
                            ) / denom );

                II = 0;
                for ( JJ = 0; JJ < 5; JJ++ ) {
                    if (JJ == Istencil) { continue; }
                    // Factor 5 is from having 5 stencils
                    // factor R*cos(lat) comes from converting y deriv to lat deriv
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
                            ) / denom );

                // x0 x2 y0 - x0 x3 y0 - x0 x2 y2 + x2 x3 y2 + x0 x3 y3 -  x2 x3 y3, 
                // x0 x2 ( y0 - y2) + x2 x3 (y2 - y3) + x0 x3 (-y0 + y3),
                ddy_coeffs.push_back( (
                            + xi[0] * xi[2] * ( yi[0] - yi[2])
                            + xi[2] * xi[3] * ( yi[2] - yi[3])
                            + xi[0] * xi[3] * (-yi[0] + yi[3])
                            ) / denom );

                // -x0 x1 y0 + x0 x3 y0 + x0 x1 y1 - x1 x3 y1 - x0 x3 y3 +  x1 x3 y3,
                // x0 x1 (-y0 + y1) + x0 x3 (y0 - y3) + x1 x3 (-y1 + y3),
                ddy_coeffs.push_back( (
                            - xi[0] * xi[1] * yi[0]
                            + xi[0] * xi[3] * yi[0]
                            + xi[0] * xi[1] * yi[1]
                            - xi[1] * xi[3] * yi[1]
                            - xi[0] * xi[3] * yi[3]
                            + xi[1] * xi[3] * yi[3]
                            ) / denom );

                // x0 x1 y0 - x0 x2 y0 - x0 x1 y1 + x1 x2 y1 + x0 x2 y2 - x1 x2 y2
                // x0 x1 ( y0 - y1) + x1 x2 (y1 - y2) + x0 x2 (-y0 + y2)
                ddy_coeffs.push_back( (
                            + xi[0] * xi[1] * yi[0]
                            - xi[0] * xi[2] * yi[0]
                            - xi[0] * xi[1] * yi[1]
                            + xi[1] * xi[2] * yi[1]
                            + xi[0] * xi[2] * yi[2]
                            - xi[1] * xi[2] * yi[2]
                            ) / denom );

                II = 0;
                for ( JJ = 0; JJ < 5; JJ++ ) {
                    if (JJ == Istencil) { continue; }
                    adjacency_ddlat_weights.at(pt_index)[JJ] += ddy_coeffs[II];
                    II++;
                }
            }
            for ( JJ = 0; JJ < 5; JJ++ ) {
                // Factor 5 is from having 5 stencils
                // factor R*cos(lat) comes from converting y deriv to lat deriv
                adjacency_ddlon_weights.at(pt_index)[JJ] *= 
                    constants::R_earth * cos(pt_lat) / stencil_count;
                // Factor 5 is from having 5 stencils
                // factor R comes from converting y deriv to lat deriv
                adjacency_ddlat_weights.at(pt_index)[JJ] *= constants::R_earth / stencil_count;
            }
        }
    }
}
