#include "../constants.hpp"
#include "../functions.hpp"
#include "../preprocess.hpp"
#include "../differentiation_tools.hpp"
#include <algorithm>
#include <vector>
#include <omp.h>
#include <math.h>
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>

// Build the LHS of the problem
void HelmholtzDataClass::Build_LHS(
        const dataset & source_data
        ) {

    typedef Eigen::Triplet<double> T;
    const bool USE_TRUE_2ND_DERIV = true;
    if (USE_TRUE_2ND_DERIV) {
        fprintf( stdout, "Using explicit second derivatives.\n" );
    } else {
        fprintf( stdout, "Using repeated first derivatives for second derivatives.\n" );
    }
    const int pts_per_1st_deriv = (constants::GRID_TYPE == constants::GridType::MeshGrid)
        ? constants::DiffOrd + 1 : constants::ADJACENCY_SIZE + 1;
    const int pts_per_2nd_deriv = (constants::GRID_TYPE == constants::GridType::MeshGrid)
        ? constants::DiffOrd + 2 : pow(constants::ADJACENCY_SIZE + 1, USE_TRUE_2ND_DERIV ? 1 : 2);
    #if DEBUG >= 1
    fprintf( stdout, "Build_LHS: Building with %d points for 1st derivs and %d for second derivs.\n",
            pts_per_1st_deriv, pts_per_2nd_deriv );
    #endif
    std::vector<T> Aij_triplets( 
            (use_vort_div and use_vel) ? 
                Nrow * ( 6 * pts_per_1st_deriv + 4 * pts_per_2nd_deriv )
                :
            ( use_vel ) ? 
                Nrow *   4 * pts_per_1st_deriv
                :
                Nrow * ( 2 * pts_per_1st_deriv + 4 * pts_per_2nd_deriv )
            ,
            T(0,0,0) );
    size_t Itriplet, Ipt, Ipt_mapped, neighbour_mapped, neighbour_ind, row_skip, column_skip,
           Ineighbour;
    int Itime, Idepth, Ilat, Ilon, counter;
    double weight_val, cos_lat_inv, R_inv, R2_inv, cos2_lat_inv, val;
    bool is_pole, neighbour_is_zero;

    const int Ntime  = source_data.Ntime,
              Ndepth = source_data.Ndepth,
              Nlat   = source_data.Nlat,
              Nlon   = source_data.Nlon;

    const size_t Npts = source_data.mask.size();
    const std::vector<short int> unmask(Npts, true);
    const size_t num_neighbours = source_data.num_neighbours;

    
    //
    //// If MeshGrid
    //
    if (constants::GRID_TYPE == constants::GridType::MeshGrid) {

        int LB, Idiff, IDIFF, Ndiff;
        std::vector<double> diff_vec;

        #pragma omp parallel default(none) \
        shared( source_data, Aij_triplets, unmask ) \
        private( Ipt, neighbour_ind, Itriplet, row_skip, column_skip, is_pole, val, \
                weight_val, cos_lat_inv, R_inv, cos2_lat_inv, R2_inv, counter, \
                Ipt_mapped, neighbour_mapped, \
                Itime, Idepth, Ilat, Ilon, LB, diff_vec, Ndiff, IDIFF, Idiff ) \
        firstprivate( Npts, Nrow, Ncol, weight_err, \
                stdout, pts_per_1st_deriv, \
                Ntime, Ndepth, Nlon, Nlat )
        { 
            #pragma omp for collapse(1) schedule(static)
            for ( Ipt = 0; Ipt < Npts; Ipt++ ) {

                if ( Ipt == 0 ) { continue; } // Force to zero at corner
                if ( all_land_neighbours[Ipt] == 1 ) { continue; } // Skip points that were mapped

                Index1to4( Ipt, Itime, Idepth, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon );
                is_pole = std::fabs( std::fabs( source_data.latitude.at(Ilat) * 180.0 / M_PI ) - 90 ) < 1e-6;
                if ( is_pole ) { fprintf(stdout, "Build_LHS: SKIPPING POLE POINT!\n"); continue; }

                weight_val = weight_err ? source_data.areas.at(Ipt) : 1.;
                cos_lat_inv = 1. / cos(source_data.latitude.at(Ilat));
                R_inv = 1. / constants::R_earth;

                cos2_lat_inv = pow(1. / cos(source_data.latitude.at(Ilat)), 2);
                R2_inv = pow(1. / constants::R_earth, 2);

                // get row index under land mapping
                Ipt_mapped = Ipt - num_mapped_before_row[Ipt]; // coast included in rows
                #if DEBUG >= 1
                if ( Ipt_mapped == 0 ) {
                    fprintf( stdout, "Build_LHS: Point %zu mapped to row zero.\n", Ipt );
                }
                #endif

                // First lon deriv contributions
                LB = - 2 * Nlon;
                get_diff_vector( diff_vec, LB, source_data.longitude, "lon", 
                        Itime, Idepth, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon, unmask, 1);
                assert( LB != -2*Nlon );
                Ndiff = ( LB == - 2 * Nlon ) ? 0 : diff_vec.size();
                assert( Ndiff == pts_per_1st_deriv );
                for ( IDIFF = LB; IDIFF < LB + Ndiff; IDIFF++ ) {

                    if (constants::PERIODIC_X) { Idiff = ( IDIFF % Nlon + Nlon ) % Nlon; }
                    else                       { Idiff = IDIFF;                          }

                    neighbour_ind = Index(0, 0, Ilat, Idiff, 1, 1, Nlat, Nlon);

                    neighbour_mapped = pt_maps_to[neighbour_ind]; // convert to mapped coordinated
                    if ( neighbour_mapped == 0 ) { continue; }
                    neighbour_mapped = neighbour_mapped - num_mapped_before_col[neighbour_mapped]; 

                    val  = diff_vec.at(IDIFF-LB) * cos_lat_inv * R_inv;
                    val *= weight_val;

                    // Psi part (of u_lat)
                    column_skip = 0 * Ncol + neighbour_mapped;
                    row_skip    = 1 * Nrow + Ipt_mapped;
                    Itriplet = Ipt_mapped * pts_per_1st_deriv + (IDIFF-LB);
                    Aij_triplets[Itriplet] = T( row_skip, column_skip, val );

                    // Phi part (of u_lon)
                    column_skip = 1 * Ncol + neighbour_mapped;
                    row_skip    = 0 * Nrow + Ipt_mapped;
                    Itriplet = (2*Nrow + Ipt_mapped) * pts_per_1st_deriv + (IDIFF-LB);
                    Aij_triplets[Itriplet] = T( row_skip, column_skip, val );
                }

                // First lat deriv contributions
                LB = - 2 * Nlat;
                get_diff_vector( diff_vec, LB, source_data.latitude, "lat", 
                        Itime, Idepth, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon, unmask, 1);
                assert( LB != -2*Nlat );
                Ndiff = ( LB == - 2 * Nlat ) ? 0 : diff_vec.size();
                assert( Ndiff == pts_per_1st_deriv );
                for ( IDIFF = LB; IDIFF < LB + Ndiff; IDIFF++ ) {

                    if (constants::PERIODIC_Y) { Idiff = ( IDIFF % Nlat + Nlat ) % Nlat; }
                    else                       { Idiff = IDIFF;                          }

                    neighbour_ind = Index(0, 0, Idiff, Ilon, 1, 1, Nlat, Nlon);

                    neighbour_mapped = pt_maps_to[neighbour_ind]; // convert to mapped coordinated
                    if ( neighbour_mapped == 0 ) { continue; }
                    neighbour_mapped = neighbour_mapped - num_mapped_before_col[neighbour_mapped]; 

                    val  = diff_vec.at(IDIFF-LB) * R_inv;
                    val *= weight_val;

                    // Psi part (of u_lon)
                    column_skip = 0 * Ncol + neighbour_mapped;
                    row_skip    = 0 * Nrow + Ipt_mapped;
                    Itriplet = (Nrow + Ipt_mapped) * pts_per_1st_deriv + (IDIFF-LB);
                    Aij_triplets[Itriplet] = T( row_skip, column_skip, -val );

                    // Phi part (of u_lat)
                    column_skip = 1 * Ncol + neighbour_mapped;
                    row_skip    = 1 * Nrow + Ipt_mapped;
                    Itriplet = (3*Nrow + Ipt_mapped) * pts_per_1st_deriv + (IDIFF-LB);
                    Aij_triplets[Itriplet] = T( row_skip, column_skip, val );

                    // Also the Laplacian contribution
                    for ( counter = 0; counter < 2; counter++ ) {
                        val = - diff_vec.at(IDIFF-LB) * tan( source_data.latitude.at(Ilat) );
                        val *= weight_val * pow(R_inv, 2.);
                        val *= Tikhov;

                        column_skip = counter * Ncol + neighbour_mapped;
                        row_skip    = (2+counter) * Nrow + Ipt_mapped;
                        Itriplet  = (counter*Nrow + Ipt_mapped) * pts_per_1st_deriv + (IDIFF-LB);
                        Itriplet += 4 * Nrow * pts_per_2nd_deriv;
                        Itriplet += 4 * Nrow * pts_per_1st_deriv;
                        Aij_triplets[Itriplet] = T( row_skip, column_skip, val );
                    }
                }

                // If we're not including vorticity and divergence in the solver, we're done
                if (not(use_vort_div)) { continue; }

                // Second lon deriv contributions
                LB = - 2 * Nlon;
                get_diff_vector( diff_vec, LB, source_data.longitude, "lon", 
                        Itime, Idepth, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon, unmask, 2);
                assert( LB != -2*Nlon );
                Ndiff = ( LB == - 2 * Nlon ) ? 0 : diff_vec.size();
                assert( Ndiff == pts_per_2nd_deriv );
                for ( IDIFF = LB; IDIFF < LB + Ndiff; IDIFF++ ) {

                    if (constants::PERIODIC_X) { Idiff = ( IDIFF % Nlon + Nlon ) % Nlon; }
                    else                       { Idiff = IDIFF;                          }

                    neighbour_ind = Index(0, 0, Ilat, Idiff, 1, 1, Nlat, Nlon);

                    neighbour_mapped = pt_maps_to[neighbour_ind]; // convert to mapped coordinated
                    if ( neighbour_mapped == 0 ) { continue; }
                    neighbour_mapped = neighbour_mapped - num_mapped_before_col[neighbour_mapped]; 

                    val  = diff_vec.at(IDIFF-LB) * cos2_lat_inv * R2_inv;
                    val *= weight_val * Tikhov;

                    for ( counter = 0; counter < 2; counter++ ) {
                        column_skip = counter * Ncol + neighbour_mapped;
                        row_skip    = (2+counter) * Nrow + Ipt_mapped;
                        Itriplet  = (counter*Nrow + Ipt_mapped) * pts_per_2nd_deriv + (IDIFF-LB);
                        Itriplet += 4 * Nrow * pts_per_1st_deriv;
                        Aij_triplets[Itriplet] = T( row_skip, column_skip, val );
                    }
                }

                // Second lat deriv contributions
                LB = - 2 * Nlat;
                get_diff_vector( diff_vec, LB, source_data.latitude, "lat", 
                        Itime, Idepth, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon, unmask, 2);
                assert( LB != -2*Nlat );
                Ndiff = ( LB == - 2 * Nlat ) ? 0 : diff_vec.size();
                assert( Ndiff == pts_per_2nd_deriv );
                for ( IDIFF = LB; IDIFF < LB + Ndiff; IDIFF++ ) {

                    if (constants::PERIODIC_Y) { Idiff = ( IDIFF % Nlat + Nlat ) % Nlat; }
                    else                       { Idiff = IDIFF;                          }

                    neighbour_ind = Index(0, 0, Idiff, Ilon, 1, 1, Nlat, Nlon);

                    neighbour_mapped = pt_maps_to[neighbour_ind]; // convert to mapped coordinated
                    if ( neighbour_mapped == 0 ) { continue; }
                    neighbour_mapped = neighbour_mapped - num_mapped_before_col[neighbour_mapped]; 

                    val  = diff_vec.at(IDIFF-LB) * R2_inv;
                    val *= weight_val * Tikhov;

                    for ( counter = 0; counter < 2; counter++ ) {
                        column_skip = counter * Ncol + neighbour_mapped;
                        row_skip    = (2+counter) * Nrow + Ipt_mapped;
                        Itriplet  = (counter*Nrow + Ipt_mapped) * pts_per_2nd_deriv + (IDIFF-LB);
                        Itriplet += 4 * Nrow * pts_per_1st_deriv;
                        Itriplet += 2 * Nrow * pts_per_2nd_deriv;
                        Aij_triplets[Itriplet] = T( row_skip, column_skip, val );
                    }
                }
            }
        }

    //
    //// If LLC-type grid
    //
    } else if (constants::GRID_TYPE == constants::GridType::LLC) {

        #pragma omp parallel default(none) \
        shared( source_data, Aij_triplets, pt_maps_to, num_mapped_before_row, num_mapped_before_col ) \
        private( Ipt, Ineighbour, neighbour_ind, Itriplet, row_skip, column_skip, is_pole, val, \
                counter, Ipt_mapped, neighbour_mapped, neighbour_is_zero ) \
        firstprivate( Npts, Nrow, Ncol, num_neighbours, weight_err, Tikhov, \
                stdout, pts_per_1st_deriv, Npts_mapped )
        { 
            #pragma omp for collapse(1) schedule(static)
            for ( Ipt = 0; Ipt < Npts; Ipt++ ) {

                if ( collapse_land ) { 
                    if ( Ipt == 0 ) { continue; } // Force to zero at first point, this fixes the '+C' unknown constant
                    if ( all_land_neighbours[Ipt] == 1 ) { continue; } // Skip land-locked points
                }

                double weight_val = weight_err ? source_data.areas.at(Ipt) : 1.;
                double cos_lat_inv = 1. / cos(source_data.latitude.at(Ipt));
                double R_inv = 1. / constants::R_earth;

                if ( pt_maps_to[Ipt] != Ipt ) { throw std::runtime_error("Attempting to build row for a mapped point."); }

                // get row index under land mapping
                Ipt_mapped = Ipt - num_mapped_before_row[Ipt]; // coast included in rows
                #if DEBUG >= 1
                if ( Ipt_mapped == 0 ) {
                    fprintf( stdout, "Build_LHS: Point %zu mapped to row zero.\n", Ipt );
                }
                #endif
                if ( ( Ipt_mapped < 0 ) or (Ipt_mapped >= Nrow) ) { 
                    fprintf( stdout, "Build_LHS: BAD POINT! %'zu - %'zu |-> %'zu\n", Ipt, num_mapped_before_row[Ipt], Ipt_mapped );
                    assert(false);
                }

                for ( Ineighbour = 0; Ineighbour < num_neighbours + 1; Ineighbour++ ) {

                    neighbour_ind = (Ineighbour < num_neighbours) ? 
                        source_data.adjacency_indices.at(Ipt).at(Ineighbour) :
                        Ipt;
                    neighbour_mapped = pt_maps_to[neighbour_ind]; // convert to mapped coordinated
                    if ( neighbour_mapped != 0 ) { 
                        neighbour_mapped = neighbour_mapped - num_mapped_before_col[neighbour_mapped]; // all land removed from columns
                        neighbour_is_zero = false;
                    } else {
                        neighbour_is_zero = true;
                    }
                    if ( ( neighbour_mapped < 0 ) or (neighbour_mapped >= Ncol) ) { 
                        fprintf( stdout, "Build_LHS: BAD Neighbour! %'zu - %'zu - 1 |-> %'zu\n", neighbour_ind, num_mapped_before_col[neighbour_ind], neighbour_mapped );
                        assert(false);
                    }

                    is_pole = std::fabs( std::fabs( source_data.latitude.at(Ipt) * 180.0 / M_PI ) - 90 ) < 1e-6;
                    if ( is_pole ) { fprintf(stdout, "Build_LHS: SKIPPING POLE POINT!\n"); continue; }

                    // Hard-setting Psi[Ipt = 0] = 0 = Phi[Ipt = 0], so it has no contribution to the LHS
                    if ( use_vel and not(neighbour_is_zero) ) {
                        // LON first derivative
                        val  = source_data.adjacency_ddlon_weights.at(Ipt).at(Ineighbour);
                        val *= weight_val * cos_lat_inv * R_inv;

                        // Psi part (of u_lat)
                        column_skip = 0 * Ncol + neighbour_mapped;
                        row_skip    = 1 * Nrow + Ipt_mapped;
                        Itriplet = Ipt_mapped * pts_per_1st_deriv + Ineighbour;
                        Aij_triplets[Itriplet] = T( row_skip, column_skip, val );

                        // Phi part (of u_lon)
                        column_skip = 1 * Ncol + neighbour_mapped;
                        row_skip    = 0 * Nrow + Ipt_mapped;
                        Itriplet = (2*Nrow + Ipt_mapped) * pts_per_1st_deriv + Ineighbour;
                        Aij_triplets[Itriplet] = T( row_skip, column_skip, val );


                        // LAT first derivative
                        val  = source_data.adjacency_ddlat_weights.at(Ipt).at(Ineighbour);
                        val *= weight_val * R_inv;

                        // Psi part (of u_lon)
                        column_skip = 0 * Ncol + neighbour_mapped;
                        row_skip    = 0 * Nrow + Ipt_mapped;
                        Itriplet = (Nrow + Ipt_mapped) * pts_per_1st_deriv + Ineighbour;
                        Aij_triplets[Itriplet] = T( row_skip, column_skip, -val );

                        // Phi part (of u_lat)
                        column_skip = 1 * Ncol + neighbour_mapped;
                        row_skip    = 1 * Nrow + Ipt_mapped;
                        Itriplet = (3*Nrow + Ipt_mapped) * pts_per_1st_deriv + Ineighbour;
                        Aij_triplets[Itriplet] = T( row_skip, column_skip, val );
                    }


                    // If we're not including vorticity and divergence in the solver, we're done
                    if (not(use_vort_div)) { continue; }

                    // And do Lap for both Phi and Psi
                    for ( counter = 0; counter < 2; counter++ ) {

                        // Second LON derivative
                        if ( USE_TRUE_2ND_DERIV ) {
                            if ( not(neighbour_is_zero) ) {
                                val  = source_data.adjacency_d2dlon2_weights.at(Ipt).at(Ineighbour);
                                val *= weight_val * pow(cos_lat_inv * R_inv, 2.);
                                val *= Tikhov;

                                column_skip = neighbour_mapped;
                                column_skip += counter * Ncol;
                                row_skip    = Ipt_mapped;
                                if ( use_vel ) { row_skip += 2 * Nrow; }
                                row_skip    += counter * Nrow;
                                Itriplet = (counter*Nrow + Ipt_mapped) * pts_per_2nd_deriv + Ineighbour;
                                if ( use_vel ) { Itriplet += 4 * Nrow * pts_per_1st_deriv; }
                                Aij_triplets[Itriplet] = T( row_skip, column_skip, val );
                            }
                        } else {
                            for ( size_t D2_ind = 0; D2_ind < num_neighbours+1; D2_ind++ ) {
                                val  =   source_data.adjacency_ddlon_weights.at(Ipt).at(Ineighbour)
                                    * source_data.adjacency_ddlon_weights.at(neighbour_ind).at(D2_ind);
                                val *= weight_val * pow(R_inv, 2.) * cos_lat_inv / cos(source_data.latitude.at(neighbour_ind));
                                val *= Tikhov;

                                column_skip = source_data.adjacency_indices.at(neighbour_ind).at(D2_ind);
                                column_skip = pt_maps_to[column_skip];
                                if (column_skip == 0) {continue;}
                                column_skip = column_skip - num_mapped_before_col[column_skip];
                                column_skip += counter * Ncol;
                                row_skip    = Ipt_mapped;
                                if ( use_vel ) { row_skip += 2 * Nrow; }
                                row_skip    += counter * Nrow;
                                Itriplet = (counter*Nrow + Ipt_mapped) * pts_per_2nd_deriv + Ineighbour*(num_neighbours+1)  + D2_ind;
                                if ( use_vel ) { Itriplet += 4 * Nrow * pts_per_1st_deriv; }
                                Aij_triplets[Itriplet] = T( row_skip, column_skip, val );
                            }
                        }

                        // Second LAT derivative
                        if ( USE_TRUE_2ND_DERIV ) {
                            if ( not(neighbour_is_zero) ) {
                                val = source_data.adjacency_d2dlat2_weights.at(Ipt).at(Ineighbour);
                                val *= weight_val * pow(R_inv, 2.);
                                val *= Tikhov;

                                column_skip = neighbour_mapped;
                                column_skip += counter * Ncol;
                                row_skip    = Ipt_mapped;
                                if ( use_vel ) { row_skip += 2 * Nrow; }
                                row_skip    += counter * Nrow;
                                Itriplet    = (counter * Nrow + Ipt_mapped) * pts_per_2nd_deriv + Ineighbour;
                                if ( use_vel ) { Itriplet += 4 * Nrow * pts_per_1st_deriv; }
                                Itriplet += 2 * Nrow * pts_per_2nd_deriv;
                                Aij_triplets[Itriplet] = T( row_skip, column_skip, val );
                            }
                        } else {
                            for ( size_t D2_ind = 0; D2_ind < num_neighbours+1; D2_ind++ ) {
                                val  =   source_data.adjacency_ddlat_weights.at(Ipt).at(Ineighbour)
                                    * source_data.adjacency_ddlat_weights.at(neighbour_ind).at(D2_ind);
                                val *= weight_val * pow(R_inv, 2.);
                                val *= Tikhov;

                                column_skip = source_data.adjacency_indices.at(neighbour_ind).at(D2_ind);
                                column_skip = pt_maps_to[column_skip];
                                if (column_skip == 0) {continue;}
                                column_skip = column_skip - num_mapped_before_col[column_skip];
                                column_skip += counter * Ncol;
                                row_skip    = Ipt_mapped;
                                if ( use_vel ) { row_skip += 2 * Nrow; }
                                row_skip    += counter * Nrow;
                                Itriplet = (counter*Nrow + Ipt_mapped) * pts_per_2nd_deriv + Ineighbour*(num_neighbours+1)  + D2_ind;
                                if ( use_vel ) { Itriplet += 4 * Nrow * pts_per_1st_deriv; }
                                Itriplet += 2 * Nrow * pts_per_2nd_deriv;
                                Aij_triplets[Itriplet] = T( row_skip, column_skip, val );
                            }
                        }

                        // First LAT derivative
                        if ( not(neighbour_is_zero) ) {
                            val = - source_data.adjacency_ddlat_weights.at(Ipt).at(Ineighbour) * tan( source_data.latitude.at(Ipt) );
                            val *= weight_val * pow(R_inv, 2.);
                            val *= Tikhov;

                            column_skip = neighbour_mapped;
                            column_skip += counter * Ncol;
                            row_skip    = Ipt_mapped;
                            if ( use_vel ) { row_skip += 2 * Nrow; }
                            row_skip    += counter * Nrow;
                            Itriplet    = (counter * Nrow + Ipt_mapped) * pts_per_1st_deriv + Ineighbour;
                            if ( use_vel ) { Itriplet += 4 * Nrow * pts_per_1st_deriv; }
                            Itriplet += 4 * Nrow * pts_per_2nd_deriv;
                            Aij_triplets[Itriplet] = T( row_skip, column_skip, val );
                        }
                    }
                }
            }
        }




    }


    // Finally, however we built the [row,col,val] triplets, 
    // use them to assemble our matrix and convert to compressed row form
    #if DEBUG >= 1
    fprintf( stdout, "  Build_LHS: LHS is %'zu x %'zu\n", (use_vort_div and use_vel) ? 4*Nrow : 2*Nrow, 2*Ncol );
    #endif
    LHS.resize( (use_vort_div and use_vel) ? 4*Nrow : 2*Nrow, 2*Ncol );
    LHS.setFromTriplets( Aij_triplets.begin(), Aij_triplets.end() );
    LHS.makeCompressed();

}
