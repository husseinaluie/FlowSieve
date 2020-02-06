#include "../constants.hpp"
#include "../functions.hpp"
#include "../differentiation_tools.hpp"
#include "../preprocess.hpp"
#include <algorithm>
#include <vector>
#include <omp.h>
#include <math.h>
#include "../ALGLIB/stdafx.h"
#include "../ALGLIB/linalg.h"
#include "../ALGLIB/solvers.h"

void toroidal_sparse_Lap(
        alglib::sparsematrix & Lap,
        const std::vector<double> & latitude,
        const std::vector<double> & longitude,
        const int Nlat,
        const int Nlon,
        const std::vector<double> & mask
        ) {

    int Ilat, Ilon, IDIFF, Idiff, Ndiff,
        diff_index, index, LB;
    double old_val, tmp, cos2_lat_inv, tan_lat;
    std::vector<double> diff_vec;

    const double R2_inv = 1. / pow(constants::R_earth, 2);

    for ( Ilat = 0; Ilat < Nlat; Ilat++ ) {
        for ( Ilon = 0; Ilon < Nlon; Ilon++ ) {

            cos2_lat_inv = 1. / pow( cos(latitude.at(Ilat)), 2 );
            tan_lat = tan(latitude.at(Ilat));

            index = Index(0, 0, Ilat, Ilon,
                          1, 1, Nlat, Nlon);

            //
            //// LON second derivative part
            //

            get_diff_vector(diff_vec, LB, longitude, "lon",
                    0, 0, Ilat, Ilon,
                    1, 1, Nlat, Nlon,
                    mask, 2, constants::DiffOrd);

            Ndiff = diff_vec.size();

            for ( IDIFF = LB; IDIFF < LB + Ndiff; IDIFF++ ) {

                if (constants::PERIODIC_X) {
                    if      (IDIFF < 0    ) { Idiff = IDIFF + Nlon; }
                    else if (IDIFF >= Nlon) { Idiff = IDIFF - Nlon; }
                    else                    { Idiff = IDIFF; }
                } else {
                    Idiff = IDIFF;
                }

                diff_index = Index(0, 0, Ilat, Idiff,
                                   1, 1, Nlat, Nlon);

                if ( (diff_index >= Nlat*Nlon) or (diff_index < 0) ) {
                    fprintf(stdout, 
                            "  LON(2) index out of bounds"
                            " : (lat,lon) = (%d, %d)"
                            " : (index, diff_index) = (%d, %d)"
                            " : IDIFF, Idiff = (%d, %d)"
                            " : LB = %d\n",
                            Ilat, Ilon, index, diff_index, IDIFF, Idiff, LB);
                    fflush(stdout);
                }

                old_val = sparseget(Lap, index, diff_index);

                tmp = diff_vec.at(IDIFF-LB) * cos2_lat_inv * R2_inv;

                alglib::sparseset(Lap, index, diff_index, old_val + tmp);
            }


            //
            //// LAT second derivative part
            //

            get_diff_vector(diff_vec, LB, latitude, "lat",
                    0, 0, Ilat, Ilon,
                    1, 1, Nlat, Nlon,
                    mask, 2, constants::DiffOrd);

            Ndiff = diff_vec.size();

            for ( IDIFF = LB; IDIFF < LB + Ndiff; IDIFF++ ) {

                if (constants::PERIODIC_Y) {
                    if      (IDIFF < 0    ) { Idiff = IDIFF + Nlat; }
                    else if (IDIFF >= Nlat) { Idiff = IDIFF - Nlat; }
                    else                    { Idiff = IDIFF; }
                } else {
                    Idiff = IDIFF;
                }

                diff_index = Index(0, 0, Idiff, Ilon,
                                   1, 1, Nlat,  Nlon);

                if ( (diff_index >= Nlat*Nlon) or (diff_index < 0) ) {
                    fprintf(stdout, 
                            "  LAT(2) index out of bounds"
                            " : (lat,lon) = (%d, %d)"
                            " : (index, diff_index) = (%d, %d)"
                            " : IDIFF, Idiff = (%d, %d)"
                            " : LB = %d\n",
                            Ilat, Ilon, index, diff_index, IDIFF, Idiff, LB);
                    fflush(stdout);
                }

                old_val = sparseget(Lap, index, diff_index);

                tmp = diff_vec.at(IDIFF-LB) * R2_inv;

                alglib::sparseset(Lap, index, diff_index, old_val + tmp);
            }


            //
            //// LAT first derivative part
            //

            get_diff_vector(diff_vec, LB, latitude, "lat",
                    0, 0, Ilat, Ilon,
                    1, 1, Nlat, Nlon,
                    mask, 1, constants::DiffOrd);

            Ndiff = diff_vec.size();

            for ( IDIFF = LB; IDIFF < LB + Ndiff; IDIFF++ ) {

                if (constants::PERIODIC_Y) {
                    if      (IDIFF < 0    ) { Idiff = IDIFF + Nlat; }
                    else if (IDIFF >= Nlat) { Idiff = IDIFF - Nlat; }
                    else                    { Idiff = IDIFF; }
                } else {
                    Idiff = IDIFF;
                }

                diff_index = Index(0, 0, Idiff, Ilon,
                                   1, 1, Nlat,  Nlon);

                if ( (diff_index >= Nlat*Nlon) or (diff_index < 0) ) {
                    fprintf(stdout, 
                            "  LAT(1) index out of bounds"
                            " : (lat,lon) = (%d, %d)"
                            " : (index, diff_index) = (%d, %d)"
                            " : IDIFF, Idiff = (%d, %d)"
                            " : LB = %d\n",
                            Ilat, Ilon, index, diff_index, IDIFF, Idiff, LB);
                    fflush(stdout);
                }

                old_val = sparseget(Lap, index, diff_index);

                tmp = - diff_vec.at(IDIFF-LB) * tan_lat * R2_inv;

                alglib::sparseset(Lap, index, diff_index, old_val + tmp);
            }

        }
    }

}
