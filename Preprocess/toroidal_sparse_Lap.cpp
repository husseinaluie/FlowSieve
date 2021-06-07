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
 * \brief Builds the (sparse!) Laplacian differentiation operator.
 * @ingroup ToroidalProjection
 *
 * Uses the differentation order specified in contants.hpp
 *
 * Currently only handles spherical coordinates.
 *
 * @param[in,out]   Lap                     Where to store the (sparse) differentiation matrix
 * @param[in]       source_data             dataset class storing various fields (longitude, latitude, etc)
 * @param[in]       Itime,Idepth            Current time-depth iteration
 * @param[in]       mask                    Array to distinguish land/water
 * @param[in]       area_weight             Bool indicating if the Laplacian should be weighted by cell-size (i.e. weight error by cell size). Default is false.
 *
 */
void toroidal_sparse_Lap(
        alglib::sparsematrix & Lap,
        const dataset & source_data,
        const int Itime,
        const int Idepth,
        const std::vector<bool>   & mask,
        const bool area_weight,
        const size_t row_skip,
        const size_t column_skip
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
    double old_val, tmp, cos2_lat_inv, tan_lat;
    std::vector<double> diff_vec;
    bool is_pole;

    const double R2_inv = 1. / pow(constants::R_earth, 2);

    for ( Ilat = 0; Ilat < Nlat; Ilat++ ) {
        for ( Ilon = 0; Ilon < Nlon; Ilon++ ) {
            
            // If we're too close to the pole (less than 0.01 degrees), bad things happen
            is_pole = std::fabs( std::fabs( latitude.at(Ilat) * 180.0 / M_PI ) - 90 ) < 0.01;

            cos2_lat_inv = 1. / pow( cos(latitude.at(Ilat)), 2 );
            tan_lat = tan(latitude.at(Ilat));

            index = Index(Itime, Idepth, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon);
            index_sub = Index(0, 0, Ilat, Ilon, 1, 1, Nlat, Nlon);

            if ( (mask.at(index)) and not(is_pole) ) { // Skip land areas and poles

                //
                //// LON second derivative part
                //

                LB = - 2 * Nlon;
                get_diff_vector(diff_vec, LB, longitude, "lon",
                                Itime, Idepth, Ilat, Ilon,
                                Ntime, Ndepth, Nlat, Nlon,
                                mask, 2, constants::DiffOrd);

                Ndiff = diff_vec.size();

                // If LB is unchanged, then we failed to build a stencil
                if (LB != - 2 * Nlon) {
                    for ( IDIFF = LB; IDIFF < LB + Ndiff; IDIFF++ ) {

                        if (constants::PERIODIC_X) { Idiff = ( IDIFF % Nlon + Nlon ) % Nlon; }
                        else                       { Idiff = IDIFF;                          }

                        diff_index = Index(0, 0, Ilat, Idiff, 1, 1, Nlat, Nlon);

                        tmp = diff_vec.at(IDIFF-LB) * cos2_lat_inv * R2_inv;
                        if (area_weight) { tmp *= areas.at(index_sub); }

                        old_val = sparseget(Lap, row_skip + index_sub, column_skip + diff_index);
                        alglib::sparseset(  Lap, row_skip + index_sub, column_skip + diff_index, old_val + tmp);
                    }
                }


                //
                //// LAT second derivative part
                //

                LB = -2 * Nlat;
                get_diff_vector(diff_vec, LB, latitude, "lat",
                                Itime, Idepth, Ilat, Ilon,
                                Ntime, Ndepth, Nlat, Nlon,
                                mask, 2, constants::DiffOrd);

                Ndiff = diff_vec.size();

                // If LB is unchanged, then we failed to build a stencil
                if (LB != - 2 * Nlat) {
                    for ( IDIFF = LB; IDIFF < LB + Ndiff; IDIFF++ ) {

                        if (constants::PERIODIC_Y) { Idiff = ( IDIFF % Nlat + Nlat ) % Nlat; }
                        else                       { Idiff = IDIFF;                          }

                        diff_index = Index(0, 0, Idiff, Ilon, 1, 1, Nlat,  Nlon);

                        tmp = diff_vec.at(IDIFF-LB) * R2_inv;
                        if (area_weight) { tmp *= areas.at(index_sub); }

                        old_val = sparseget(Lap, row_skip + index_sub, column_skip + diff_index);
                        alglib::sparseset(  Lap, row_skip + index_sub, column_skip + diff_index, old_val + tmp);
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

                        tmp = - diff_vec.at(IDIFF-LB) * tan_lat * R2_inv;
                        if (area_weight) { tmp *= areas.at(index_sub); }

                        old_val = sparseget(Lap, row_skip + index_sub, column_skip + diff_index);
                        alglib::sparseset(  Lap, row_skip + index_sub, column_skip + diff_index, old_val + tmp);
                    }
                }
            } else { // end mask if
                // If this spot is masked, then set the value to 1
                //   if we correspondingly set the RHS value to 0,
                //   then this should force a zero value over land
                alglib::sparseset(Lap, row_skip + index_sub, column_skip + index_sub, 1.);
            }
        }
    }
}
