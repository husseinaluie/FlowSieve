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

// Build the LHS and RHS of the problem
void HelmholtzDataClass::Build_Scalar_LHS_and_RHS(
        const dataset & source_data
        ) {

    typedef Eigen::Triplet<double> T;
    //const bool USE_TRUE_2ND_DERIV = true;
    const int pts_per_1st_deriv = (constants::GRID_TYPE == constants::GridType::MeshGrid)
        ? constants::DiffOrd + 1 : constants::ADJACENCY_SIZE + 1;
    const int pts_per_2nd_deriv = (constants::GRID_TYPE == constants::GridType::MeshGrid)
        ? constants::DiffOrd + 2 : constants::ADJACENCY_SIZE + 1;
    const int pts_per_deriv = std::max( pts_per_1st_deriv, pts_per_2nd_deriv );

    std::vector<T> Aij_triplets( 3 * Nrow * 3 * pts_per_deriv, T(0,0,0) );
    size_t Itriplet, Ipt, Ipt_mapped, neighbour_mapped, neighbour_ind, row_skip, column_skip,
           I_neighbour;
    int Itime, Idepth, Ilat, Ilon, counter;
    double weight_val, cos_lat_inv, R_inv, R2_inv, cos2_lat_inv, val, cos_lat, tan_lat;
    bool is_pole, neighbour_is_zero;

    const size_t Npts = source_data.mask.size();
    const std::vector<short int> unmask(Npts, true);
    const size_t num_neighbours = source_data.num_neighbours;

    //
    RHS.resize( 3 * Nrow );

    
    //
    //// If MeshGrid
    //
    if (constants::GRID_TYPE == constants::GridType::MeshGrid) {

        const int Ntime  = source_data.Ntime,
                  Ndepth = source_data.Ndepth,
                  Nlat   = source_data.Nlat,
                  Nlon   = source_data.Nlon;

        size_t Ipt_S, Ipt_N, Ipt_E, Ipt_W;
        bool mask_S, mask_N, mask_E, mask_W;


        int LB, Idiff, IDIFF, Ndiff;
        std::vector<double> diff_vec;

        #pragma omp parallel default(none) \
        shared( source_data, Aij_triplets, unmask ) \
        private( Ipt, I_neighbour, neighbour_ind, Itriplet, row_skip, column_skip, is_pole, val, \
                weight_val, cos_lat_inv, R_inv, cos2_lat_inv, R2_inv, counter, \
                Ipt_mapped, neighbour_mapped, cos_lat, tan_lat, \
                Ipt_S, Ipt_N, Ipt_E, Ipt_W, mask_S, mask_N, mask_E, mask_W, \
                Itime, Idepth, Ilat, Ilon, LB, diff_vec, Ndiff, IDIFF, Idiff ) \
        firstprivate( Npts, Nrow, Ncol, weight_err, \
                stdout, pts_per_deriv, \
                Ntime, Ndepth, Nlon, Nlat )
        { 
            #pragma omp for collapse(1) schedule(static)
            for ( Ipt = 0; Ipt < Npts; Ipt++ ) {

                if ( source_data.mask[Ipt] ) { continue; } // skip water points

                Index1to4( Ipt, Itime, Idepth, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon );

                if ( Ilat > 0 ) {
                    Ipt_S = Index( 0, 0, Ilat - 1, Ilon, Ntime, Ndepth, Nlat, Nlon );
                } else {
                    Ipt_S = Ipt;
                }

                if ( Ilat < Nlat - 1 ) {
                    Ipt_N = Index( 0, 0, Ilat + 1, Ilon, Ntime, Ndepth, Nlat, Nlon );
                } else {
                    Ipt_N = Ipt;
                }

                Ipt_W = Index( 0, 0, Ilat, (Ilon-1+Nlon)%Nlon, Ntime, Ndepth, Nlat, Nlon );
                Ipt_E = Index( 0, 0, Ilat, (Ilon+1+Nlon)%Nlon, Ntime, Ndepth, Nlat, Nlon );

                mask_S = source_data.mask[Ipt_S];
                mask_N = source_data.mask[Ipt_N];
                mask_E = source_data.mask[Ipt_E];
                mask_W = source_data.mask[Ipt_W];

                cos_lat = cos( source_data.latitude.at(Ilat) );
                cos_lat_inv = 1./cos_lat;
                cos2_lat_inv = pow(1. / cos(source_data.latitude.at(Ilat)), 2);

                tan_lat = tan( source_data.latitude.at(Ilat) );

                R_inv = 1. / constants::R_earth;
                R2_inv = pow(1. / constants::R_earth, 2);

                weight_val = 1.;

                bool do_ddlat = false, do_ddlon = false, do_Lap = false;

                if ( mask_S and mask_N and mask_E and mask_W ) // surrounded by water
                {
                    // In this case want to force land equals mean of neighbours.
                    // just use Lap for that for now
                    do_Lap = true;
                } else if ( ( mask_S and mask_N and (mask_E or mask_W) ) // one land neighbour (E/W)
                        or ( (mask_S or mask_N) and mask_E and mask_W ) // one land neighbour (N/S)
                 ) {
                    // Set both components of gradient to zero
                    do_ddlat = true;
                    do_ddlon = true;
                } else if ( ( mask_S and mask_W ) or ( mask_N and mask_E) ) {
                    // SW/NE corner
                    do_ddlat = true;
                    do_ddlon = true;
                    //ddlat_weight = -cos_lat / sqrt( 1 + pow(cos_lat, 2) );
                    //ddlon_weight = -1.      / sqrt( 1 + pow(cos_lat, 2) );
                } else if ( ( mask_N and mask_W ) or ( mask_S or mask_E ) ) {
                    // NW/SE corner, so angle the gradient
                    do_ddlat = true;
                    do_ddlon = true;
                    //ddlat_weight =  cos_lat / sqrt( 1 + pow(cos_lat, 2) );
                    //ddlon_weight = -1.      / sqrt( 1 + pow(cos_lat, 2) );
                } else if ( mask_E or mask_W ) {
                    // water to E/W, so set E-W gradient zero
                    do_ddlon = true;
                } else if ( mask_N or mask_S ) {
                    // water to N/S, so set N-S gradient zero
                    do_ddlat = true;
                } else {
                    // All neighbours are land, so set Laplacian to zero
                    do_Lap = true;
                }

                if ( not(do_ddlat) and not(do_ddlon) and not(do_Lap) ) {
                    throw std::runtime_error("All derivative weights are zero!");
                }

                bool did_something = false;

                row_skip = num_land_before.at(Ipt);
                if ( row_skip >= Nrow ) { throw std::runtime_error("Row index too large."); }

                // First lon deriv contributions
                if ( do_ddlon ) {
                    LB = - 2 * Nlon;
                    get_diff_vector( diff_vec, LB, source_data.longitude, "lon", 
                            Itime, Idepth, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon, unmask, 1);
                    assert( LB != -2*Nlon );
                    Ndiff = ( LB == - 2 * Nlon ) ? 0 : diff_vec.size();
                    assert( Ndiff == pts_per_1st_deriv );
                    for ( IDIFF = LB; IDIFF < LB + Ndiff; IDIFF++ ) {

                        if (constants::PERIODIC_X) { Idiff = ( IDIFF % Nlon + Nlon ) % Nlon; }
                        else                       { Idiff = IDIFF;                          }

                        I_neighbour = Index( 0, 0, Ilat, Idiff, 1, 1, Nlat, Nlon );

                        val  = diff_vec.at(IDIFF-LB) * cos_lat_inv * R_inv;
                        val *= weight_val;

                        if ( source_data.mask[I_neighbour] ) {
                            // If neighbour is water, it goes on RHS
                            RHS[row_skip] = 0.;
                            RHS[row_skip] += - val * source_data.variables.at("scalar").at(I_neighbour);
                        } else {
                            column_skip = num_land_before.at(I_neighbour);
                            if ( column_skip >= Ncol ) { throw std::runtime_error("Col index too large."); }
                            Itriplet = row_skip * pts_per_deriv + (IDIFF-LB);
                            Aij_triplets[Itriplet] = T( row_skip, column_skip, val );
                        }
                    }
                    did_something = true;
                }

                // First lat deriv contributions
                row_skip += Nrow; // increment to next block
                if ( do_ddlat ) {
                    LB = - 2 * Nlat;
                    get_diff_vector( diff_vec, LB, source_data.latitude, "lat", 
                            Itime, Idepth, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon, unmask, 1);
                    assert( LB != -2*Nlat );
                    Ndiff = ( LB == - 2 * Nlat ) ? 0 : diff_vec.size();
                    assert( Ndiff == pts_per_1st_deriv );
                    for ( IDIFF = LB; IDIFF < LB + Ndiff; IDIFF++ ) {

                        Idiff = IDIFF;
                        I_neighbour = Index( 0, 0, Idiff, Ilon, 1, 1, Nlat, Nlon );

                        val  = diff_vec.at(IDIFF-LB) * R_inv;
                        val *= weight_val;

                        if ( source_data.mask[I_neighbour] ) {
                            // If neighbour is water, it goes on RHS
                            RHS[row_skip] += - val * source_data.variables.at("scalar").at(I_neighbour);
                        } else {
                            column_skip = num_land_before.at(I_neighbour);
                            if ( column_skip >= Ncol ) { throw std::runtime_error("Col index too large."); }
                            Itriplet = (Nrow + row_skip) * pts_per_deriv + (IDIFF-LB);
                            Aij_triplets[Itriplet] = T( row_skip, column_skip, val );
                        }
                    }
                    did_something = true;
                }

                // Laplacian
                row_skip += Nrow; // increment to next block
                if ( do_Lap ) {

                    // First lat component
                    LB = - 2 * Nlat;
                    get_diff_vector( diff_vec, LB, source_data.latitude, "lat", 
                            Itime, Idepth, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon, unmask, 1);
                    assert( LB != -2*Nlat );
                    Ndiff = ( LB == - 2 * Nlat ) ? 0 : diff_vec.size();
                    assert( Ndiff == pts_per_1st_deriv );
                    for ( IDIFF = LB; IDIFF < LB + Ndiff; IDIFF++ ) {

                        Idiff = IDIFF;
                        I_neighbour = Index( 0, 0, Idiff, Ilon, 1, 1, Nlat, Nlon );

                        val  = -diff_vec.at(IDIFF-LB) * tan_lat;
                        val *= weight_val * R2_inv;

                        if ( source_data.mask[I_neighbour] ) {
                            // If neighbour is water, it goes on RHS
                            RHS[row_skip] += - val * source_data.variables.at("scalar").at(I_neighbour);
                        } else {
                            column_skip = num_land_before.at(I_neighbour);
                            if ( column_skip >= Ncol ) { throw std::runtime_error("Col index too large."); }
                            Itriplet = row_skip * pts_per_deriv + (IDIFF-LB);
                            Aij_triplets[Itriplet] = T( row_skip, column_skip, val );
                        }
                    }
                    
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

                        I_neighbour = Index( 0, 0, Ilat, Idiff, 1, 1, Nlat, Nlon );

                        val  = diff_vec.at(IDIFF-LB) * cos2_lat_inv * R2_inv;
                        val *= weight_val;

                        if ( source_data.mask[I_neighbour] ) {
                            // If neighbour is water, it goes on RHS
                            RHS[row_skip] += - val * source_data.variables.at("scalar").at(I_neighbour);
                        } else {
                            column_skip = num_land_before.at(I_neighbour);
                            if ( column_skip >= Ncol ) { throw std::runtime_error("Col index too large."); }
                            Itriplet = (Nrow + row_skip) * pts_per_deriv + (IDIFF-LB);
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

                        Idiff = IDIFF;
                        I_neighbour = Index( 0, 0, Idiff, Ilon, 1, 1, Nlat, Nlon );

                        val  = diff_vec.at(IDIFF-LB) * R2_inv;
                        val *= weight_val;

                        if ( source_data.mask[I_neighbour] ) {
                            // If neighbour is water, it goes on RHS
                            RHS[row_skip] += - val * source_data.variables.at("scalar").at(I_neighbour);
                        } else {
                            column_skip = num_land_before.at(I_neighbour);
                            if ( column_skip >= Ncol ) { throw std::runtime_error("Col index too large."); }
                            Itriplet = (2*Nrow + row_skip) * pts_per_deriv + (IDIFF-LB);
                            Aij_triplets[Itriplet] = T( row_skip, column_skip, val );
                        }
                    }
                    did_something = true;
                }

                if (not(did_something)) {
                    throw std::runtime_error("Did nothing!");
                }
            }
        }

    //
    //// If LLC-type grid
    //
    } else if (constants::GRID_TYPE == constants::GridType::LLC) {

        throw std::runtime_error("In wrong block! (LLC)");

        #pragma omp parallel default(none) \
        shared( source_data, Aij_triplets, pt_maps_to, num_mapped_before_row, num_mapped_before_col ) \
        private( Ipt, I_neighbour, neighbour_ind, Itriplet, row_skip, column_skip, is_pole, val, \
                weight_val, cos_lat_inv, R_inv, counter, \
                Ipt_mapped, neighbour_mapped, neighbour_is_zero ) \
        firstprivate( Npts, Nrow, Ncol, num_neighbours, weight_err, Tikhov, \
                stdout, pts_per_1st_deriv, Npts_mapped )
        { 
            #pragma omp for collapse(1) schedule(static)
            for ( Ipt = 0; Ipt < Npts; Ipt++ ) {

                if ( source_data.mask[Ipt] ) { continue; } // Skip water points


                weight_val = weight_err ? source_data.areas.at(Ipt) : 1.;
                cos_lat_inv = 1. / cos(source_data.latitude.at(Ipt));
                R_inv = 1. / constants::R_earth;

                // If interior of land, then just set Laplacian to zero
                if ( all_land_neighbours[Ipt] == 1 ) {

                } else {
                    // Otherwise, we're coastal. Here we set the normal gradient to zero
                }

            }
        }


    }


    // Finally, however we built the [row,col,val] triplets, 
    // use them to assemble our matrix and convert to compressed row form
    fprintf( stdout, "  LHS is %'zu x %'zu\n", 3 * Nrow, Ncol );
    LHS.resize( 3 * Nrow, Ncol );
    LHS.setFromTriplets( Aij_triplets.begin(), Aij_triplets.end() );
    LHS.makeCompressed();

    bool all_row_zero = true;
    bool all_col_zero = true;
    bool all_val_zero = true;
    for ( Itriplet = 0; Itriplet < Aij_triplets.size(); Itriplet++ ) {

        if ( Aij_triplets[Itriplet].row() != 0 ) {
            all_row_zero = false;
        }
        if ( Aij_triplets[Itriplet].col() != 0 ) {
            all_col_zero = false;
        }
        if ( Aij_triplets[Itriplet].value() != 0 ) {
            all_val_zero = false;
        }
    }
    if ( all_row_zero ) { fprintf( stdout, "All row indices are zero!\n" ); }
    if ( all_col_zero ) { fprintf( stdout, "All column indices are zero!\n" ); }
    if ( all_val_zero ) { fprintf( stdout, "All entry values are zero!\n" ); }


    // Compute row and column norms
    std::vector<double> col_norms(Ncol,0), row_norms(3*Nrow,0);
    size_t col = 0, row = 0;
    for ( int k = 0; k < LHS.outerSize(); ++k ) {
        for ( Eigen::SparseMatrix<double>::InnerIterator it(LHS,k); it; ++it ) {
            val = it.value();

            row = it.row();   // row index
            row_norms[row] += pow(val,2) / Nrow;

            col = it.col();   // col index (here it is equal to k)
            col_norms[col] += pow(val,2) / Nrow;
        }
    }
    double min_row_norm = sqrt(row_norms[0]), max_row_norm = 0,
           min_col_norm = sqrt(col_norms[0]), max_col_norm = 0;
    size_t num_zero_rows = 0;
    for ( row = 0; row < 3*Nrow; row++ ) {
        row_norms[row] = sqrt(row_norms[row]);
        if ( row_norms[row] == 0 ) {
            num_zero_rows++;
            //fprintf( stdout, "!! Row %zu has zero norm!!\n", row );
        }
        min_row_norm = std::fmin( min_row_norm, row_norms[row] );
        max_row_norm = std::fmax( max_row_norm, row_norms[row] );
    }
    for ( col = 0; col < Ncol; col++ ) {
        col_norms[col] = sqrt(col_norms[col]);
        min_col_norm = std::fmin( min_col_norm, col_norms[col] );
        max_col_norm = std::fmax( max_col_norm, col_norms[col] );
    }

    if ( num_zero_rows > 0 ) {
        fprintf(stdout, "%'zu of %'zu rows are all zeros.\n", num_zero_rows, Nrow );
        //throw std::runtime_error("LHS has all-zero rows!");
    }

    fprintf(stdout, "Column norms were bounded between %e and %e.\n", min_col_norm, max_col_norm);
    fprintf(stdout, "Row norms were bounded between %e and %e.\n", min_row_norm, max_row_norm);


}
