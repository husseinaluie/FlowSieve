#include "../constants.hpp"
#include "../functions.hpp"
#include "../differentiation_tools.hpp"
#include "../preprocess.hpp"
#include <algorithm>
#include <vector>
#include <omp.h>
#include <math.h>
#include <cassert>
#include "../ALGLIB/stdafx.h"
#include "../ALGLIB/linalg.h"
#include "../ALGLIB/solvers.h"

/*!
 * \brief Builds the (sparse) matrix to extract u,v from Psi, Phi
 *
 * Is essentially a sparse differentiation matrix with appropriate kronecker products
 *
 * @param[in,out]   LHS_matr        Where to store the (sparse) matrix
 * @param[in]       source_data
 * @param[in]       Itime,Idepth    Indices to denote which time/depth index we're doing
 * @param[in]       mask            array to distinguish land/water
 * @param[in]       area_weight     Boolean indicating if rows in least-squares problem should be area-weighted
 *
 */

void sparse_vel_from_PsiPhi(
        alglib::sparsematrix & LHS_matr,
        const dataset & source_data,
        const int Itime,
        const int Idepth,
        const std::vector<bool> & mask,
        const bool area_weight
        ) {

    const std::vector<double>   &latitude   = source_data.latitude,
                                &longitude  = source_data.longitude,
                                &areas      = source_data.areas;

    const std::vector<int>  &myCounts = source_data.myCounts;

    const int   Ntime   = myCounts.at(0),
                Ndepth  = myCounts.at(1),
                Nlat    = myCounts.at(2),
                Nlon    = myCounts.at(3);

    int Ilat, Ilon, IDIFF, Idiff, Ndiff, LB;
    size_t index, index_sub, diff_index;
    const size_t Npts = Nlat * Nlon;
    double old_val, tmp, cos_lat_inv, tan_lat;
    std::vector<double> diff_vec;
    bool is_pole;

    const double R_inv = 1. / constants::R_earth;

    for ( Ilat = 0; Ilat < Nlat; Ilat++ ) {
        for ( Ilon = 0; Ilon < Nlon; Ilon++ ) {
            
            // If we're too close to the pole (less than 0.01 degrees), bad things happen
            is_pole = std::fabs( std::fabs( latitude.at(Ilat) * 180.0 / M_PI ) - 90 ) < 0.01;

            index = Index(Itime, Idepth, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon);
            index_sub = Index(0, 0, Ilat, Ilon, 1, 1, Nlat, Nlon);

            cos_lat_inv = 1. / cos(latitude.at(Ilat));

            //if ( not(is_pole) ) { // Skip land areas and poles

                //
                //// LON first derivative part
                //

                LB = - 2 * Nlon;
                get_diff_vector(diff_vec, LB, longitude, "lon",
                                Itime, Idepth, Ilat, Ilon,
                                Ntime, Ndepth, Nlat, Nlon,
                                mask, 1, constants::DiffOrd);

                Ndiff = diff_vec.size();

                // If LB is unchanged, then we failed to build a stencil
                if (LB != - 2 * Nlon) {
                    for ( IDIFF = LB; IDIFF < LB + Ndiff; IDIFF++ ) {

                        if (constants::PERIODIC_X) { Idiff = ( IDIFF % Nlon + Nlon ) % Nlon; }
                        else                       { Idiff = IDIFF;                          }

                        diff_index = Index(0, 0, Ilat, Idiff, 1, 1, Nlat, Nlon);

                        tmp = diff_vec.at(IDIFF-LB) * cos_lat_inv * R_inv;
                        if (area_weight) { tmp *= areas.at(index_sub); }

                        // Psi part
                        size_t  column_skip = 0,
                                row_skip    = Npts;
                        old_val = sparseget(LHS_matr, row_skip + index_sub, column_skip + diff_index);
                        alglib::sparseset(  LHS_matr, row_skip + index_sub, column_skip + diff_index, old_val + tmp);

                        // Phi part
                        column_skip = Npts;
                        row_skip    = 0;
                        old_val = sparseget(LHS_matr, row_skip + index_sub, column_skip + diff_index);
                        alglib::sparseset(  LHS_matr, row_skip + index_sub, column_skip + diff_index, old_val + tmp);
                    }
                }


                //
                //// LAT first derivative part
                //

                LB = - 2 * Nlat;
                get_diff_vector(diff_vec, LB, latitude, "lat",
                                Itime, Idepth, Ilat, Ilon,
                                Ntime, Ndepth, Nlat, Nlon,
                                mask, 1, constants::DiffOrd);

                Ndiff = diff_vec.size();

                // If LB is unchanged, then we failed to build a stencil
                if (LB != - 2 * Nlat) {
                    for ( IDIFF = LB; IDIFF < LB + Ndiff; IDIFF++ ) {

                        if (constants::PERIODIC_Y) { Idiff = ( IDIFF % Nlat + Nlat ) % Nlat; }
                        else                       { Idiff = IDIFF;                          }

                        diff_index = Index(0, 0, Idiff, Ilon, 1, 1, Nlat,  Nlon);

                        tmp = diff_vec.at(IDIFF-LB) * R_inv;
                        if (area_weight) { tmp *= areas.at(index_sub); }

                        // Psi part
                        size_t  column_skip = 0,
                                row_skip    = 0;
                        old_val = sparseget(LHS_matr, row_skip + index_sub, column_skip + diff_index);
                        alglib::sparseset(  LHS_matr, row_skip + index_sub, column_skip + diff_index, old_val - tmp);

                        // Phi part
                        column_skip = Npts;
                        row_skip    = Npts;
                        old_val = sparseget(LHS_matr, row_skip + index_sub, column_skip + diff_index);
                        alglib::sparseset(  LHS_matr, row_skip + index_sub, column_skip + diff_index, old_val + tmp);
                    }
                }
            //}
        }
    }
}
