#include "../constants.hpp"
#include "../functions.hpp"
#include "../differentiation_tools.hpp"
#include <algorithm>
#include <vector>
#include <omp.h>
#include <math.h>
#include <cassert>
#include "../ALGLIB/stdafx.h"
#include "../ALGLIB/linalg.h"
#include "../ALGLIB/solvers.h"

void build_Helm_matrix(
        alglib::sparsematrix & matr,
        const dataset & source_data,
        const int Itime,
        const int Idepth
        ) {

    //      Ordering is: 
    // [ 1   (sec(lat))^2 d^2dlon^2 - tan(lat) ddlat   - sec(lat)( tan(lat) ddlon + d^2dlatdlon )                      ]   [ v_r   ]   [ u_lon * u_lon ]
    // [ 0   sec(lat)( tan(lat) ddon + d^2dlatdlon )   0.5 * ( -tan(lat)ddlat - d^2ddlat^2 + (sec(lat))^2 d^2ddlon^2 ) ] * [ Phi_v ] = [ u_lon * u_lat ]
    // [ 1   d^2ddlat^2                                  sec(lat)( tan(lat) ddon + d^2dlatdlon )                       ]   [ Psi_v ]   [ u_lat * u_lat ]

    // Create some tidy names for variables
    const std::vector<double>   &latitude   = source_data.latitude,
                                &longitude  = source_data.longitude,
                                &dAreas     = source_data.areas;
    const std::vector<bool> &mask = source_data.mask;
    const std::vector<bool> unmask(mask.size(), true);

    const std::vector<int>  &myCounts = source_data.myCounts;

    const int   Ntime   = myCounts.at(0),
                Ndepth  = myCounts.at(1),
                Nlat    = myCounts.at(2),
                Nlon    = myCounts.at(3);

    const size_t Npts = Nlat * Nlon;

    double tmp_val, weight_val;
    size_t column_skip, row_skip, index_sub;
    int IDIFF, Idiff, Idiff_lat, Idiff_lon, 
               LB,    LB_lon,    LB_lat, 
               Ndiff, Ndiff_lon, Ndiff_lat;
    std::vector<double> diff_vec, diff_vec_lon, diff_vec_lat;
    for ( int Ilat = 0; Ilat < Nlat; Ilat++ ) {

        double tan_lat = tan(latitude.at(Ilat));
        double cos_lat = cos(latitude.at(Ilat));

        for ( int Ilon = 0; Ilon < Nlon; Ilon++ ) {
            index_sub = Index(0, 0, Ilat, Ilon, 1, 1, Nlat, Nlon);

            weight_val = 1.;  // don't weight this matrix

            // First longitude derivatives
            LB = - 2 * Nlon;
            get_diff_vector(diff_vec, LB, longitude, "lon", Itime, Idepth, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon, unmask, 1, constants::DiffOrd);
            Ndiff = diff_vec.size();

            if ( LB != - 2 * Nlon) {
                for ( IDIFF = LB; IDIFF < LB + Ndiff; IDIFF++ ) {
                    if (constants::PERIODIC_X) { Idiff = ( IDIFF % Nlon + Nlon ) % Nlon; }
                    else                       { Idiff = IDIFF;                          }
                    size_t diff_index = Index(0, 0, Ilat, Idiff, 1, 1, Nlat, Nlon);

                    // (0,2) entry
                    row_skip    = 0 * Npts;
                    column_skip = 2 * Npts;
                    tmp_val     = - tan_lat * diff_vec.at( IDIFF - LB ) / cos_lat;
                    tmp_val    *= weight_val;
                    alglib::sparseadd(  matr, row_skip + index_sub, column_skip + diff_index, tmp_val);

                    // (1,1) entry
                    row_skip    = 1 * Npts;
                    column_skip = 1 * Npts;
                    tmp_val     = tan_lat * diff_vec.at( IDIFF - LB ) / cos_lat;
                    tmp_val    *= weight_val;
                    alglib::sparseadd(  matr, row_skip + index_sub, column_skip + diff_index, tmp_val);

                    // (2,2) entry
                    row_skip    = 2 * Npts;
                    column_skip = 2 * Npts;
                    tmp_val     = tan_lat * diff_vec.at( IDIFF - LB ) / cos_lat;
                    tmp_val    *= weight_val;
                    alglib::sparseadd(  matr, row_skip + index_sub, column_skip + diff_index, tmp_val);
                }
            }

            // Second longitude derivatives
            LB = - 2 * Nlon;
            get_diff_vector(diff_vec, LB, longitude, "lon", Itime, Idepth, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon, unmask, 2, constants::DiffOrd);
            Ndiff = diff_vec.size();

            if ( LB != - 2 * Nlon) {
                for ( IDIFF = LB; IDIFF < LB + Ndiff; IDIFF++ ) {
                    if (constants::PERIODIC_X) { Idiff = ( IDIFF % Nlon + Nlon ) % Nlon; }
                    else                       { Idiff = IDIFF;                          }
                    size_t diff_index = Index(0, 0, Ilat, Idiff, 1, 1, Nlat, Nlon);

                    // (0,1) entry
                    row_skip    = 0 * Npts;
                    column_skip = 1 * Npts;
                    tmp_val     = diff_vec.at( IDIFF - LB ) / pow(cos_lat, 2.);
                    tmp_val    *= weight_val;
                    alglib::sparseadd(  matr, row_skip + index_sub, column_skip + diff_index, tmp_val);

                    // (1,2) entry
                    row_skip    = 1 * Npts;
                    column_skip = 2 * Npts;
                    tmp_val     = 0.5 * diff_vec.at( IDIFF - LB ) / pow(cos_lat, 2.);
                    tmp_val    *= weight_val;
                    alglib::sparseadd(  matr, row_skip + index_sub, column_skip + diff_index, tmp_val);
                }
            }

            // First latitude derivatives
            LB = - 2 * Nlat;
            get_diff_vector(diff_vec, LB, latitude, "lat", Itime, Idepth, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon, unmask, 1, constants::DiffOrd);
            Ndiff = diff_vec.size();

            if ( LB != - 2 * Nlat) {
                for ( IDIFF = LB; IDIFF < LB + Ndiff; IDIFF++ ) {
                    if (constants::PERIODIC_Y) { Idiff = ( IDIFF % Nlat + Nlat ) % Nlat; }
                    else                       { Idiff = IDIFF;                          }
                    size_t diff_index = Index(0, 0, Idiff, Ilon, 1, 1, Nlat, Nlon);

                    // (0,1) entry  
                    row_skip    = 0 * Npts;
                    column_skip = 1 * Npts;
                    tmp_val     = - tan_lat * diff_vec.at( IDIFF - LB );
                    tmp_val    *= weight_val;
                    alglib::sparseadd(  matr, row_skip + index_sub, column_skip + diff_index, tmp_val);

                    // (1,2) entry
                    row_skip    = 1 * Npts;
                    column_skip = 2 * Npts;
                    tmp_val     = - 0.5 * tan_lat * diff_vec.at( IDIFF - LB );
                    tmp_val    *= weight_val;
                    alglib::sparseadd(  matr, row_skip + index_sub, column_skip + diff_index, tmp_val);
                }
            }

            // Second latitude derivatives
            LB = - 2 * Nlat;
            get_diff_vector(diff_vec, LB, latitude, "lat", Itime, Idepth, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon, unmask, 2, constants::DiffOrd);
            Ndiff = diff_vec.size();

            if ( LB != - 2 * Nlat) {
                for ( IDIFF = LB; IDIFF < LB + Ndiff; IDIFF++ ) {
                    if (constants::PERIODIC_Y) { Idiff = ( IDIFF % Nlat + Nlat ) % Nlat; }
                    else                       { Idiff = IDIFF;                          }
                    size_t diff_index = Index(0, 0, Idiff, Ilon, 1, 1, Nlat, Nlon);

                    // (1,2) entry  
                    row_skip    = 1 * Npts;
                    column_skip = 2 * Npts;
                    tmp_val     = - 0.5 * diff_vec.at( IDIFF - LB );
                    tmp_val    *= weight_val;
                    alglib::sparseadd(  matr, row_skip + index_sub, column_skip + diff_index, tmp_val);

                    // (2,1) entry
                    row_skip    = 2 * Npts;
                    column_skip = 1 * Npts;
                    tmp_val     = diff_vec.at( IDIFF - LB );
                    tmp_val    *= weight_val;
                    alglib::sparseadd(  matr, row_skip + index_sub, column_skip + diff_index, tmp_val);
                }
            }

            // Mixed partial (first lon and first lat)
            LB_lon = - 2 * Nlon;
            get_diff_vector(diff_vec_lon, LB_lon, longitude, "lon", Itime, Idepth, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon, unmask, 1, constants::DiffOrd);
            Ndiff_lon = diff_vec_lon.size();

            LB_lat = - 2 * Nlat;
            get_diff_vector(diff_vec_lat, LB_lat, latitude, "lat", Itime, Idepth, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon, unmask, 1, constants::DiffOrd);
            Ndiff_lat = diff_vec_lat.size();

            if ( ( LB_lat != - 2 * Nlat) and ( LB_lon != - 2 * Nlon) ) {
                for ( int IDIFF_lat = LB_lat; IDIFF_lat < LB_lat + Ndiff_lat; IDIFF_lat++ ) {

                    if (constants::PERIODIC_Y) { Idiff_lat = ( IDIFF_lat % Nlat + Nlat ) % Nlat; }
                    else                       { Idiff_lat = IDIFF_lat;                          }

                    for ( int IDIFF_lon = LB_lon; IDIFF_lon < LB_lon + Ndiff_lon; IDIFF_lon++ ) {

                        if (constants::PERIODIC_X) { Idiff_lon = ( IDIFF_lon % Nlon + Nlon ) % Nlon; }
                        else                       { Idiff_lon = IDIFF_lon;                          }

                        size_t diff_index = Index(0, 0, Idiff_lat, Idiff_lon, 1, 1, Nlat, Nlon);

                        // (0,2) entry  
                        row_skip    = 0 * Npts;
                        column_skip = 2 * Npts;
                        tmp_val     = - diff_vec_lon.at( IDIFF_lon - LB_lon ) * diff_vec_lat.at( IDIFF_lat - LB_lat ) / cos_lat;
                        tmp_val    *= weight_val;
                        alglib::sparseadd(  matr, row_skip + index_sub, column_skip + diff_index, tmp_val);

                        // (1,1) entry  
                        row_skip    = 1 * Npts;
                        column_skip = 1 * Npts;
                        tmp_val     = diff_vec_lon.at( IDIFF_lon - LB_lon ) * diff_vec_lat.at( IDIFF_lat - LB_lat ) / cos_lat;
                        tmp_val    *= weight_val;
                        alglib::sparseadd(  matr, row_skip + index_sub, column_skip + diff_index, tmp_val);

                        // (2,2) entry  
                        row_skip    = 2 * Npts;
                        column_skip = 2 * Npts;
                        tmp_val     = diff_vec_lon.at( IDIFF_lon - LB_lon ) * diff_vec_lat.at( IDIFF_lat - LB_lat ) / cos_lat;
                        tmp_val    *= weight_val;
                        alglib::sparseadd(  matr, row_skip + index_sub, column_skip + diff_index, tmp_val);

                    }
                }
            }


            // (0,0)
            row_skip    = 0 * Npts;
            column_skip = 0 * Npts;
            tmp_val     = 1.;
            tmp_val    *= weight_val;
            alglib::sparseadd(  matr, row_skip + index_sub, column_skip + index_sub, tmp_val);
            
            // (2,0)
            row_skip    = 2 * Npts;
            column_skip = 0 * Npts;
            tmp_val     = 1.;
            tmp_val    *= weight_val;
            alglib::sparseadd(  matr, row_skip + index_sub, column_skip + index_sub, tmp_val);

        }
    }

    size_t lower_count = alglib::sparsegetlowercount( matr );
    size_t upper_count = alglib::sparsegetuppercount( matr );

    alglib::sparseconverttocrs( matr );

}

void uiuj_from_Helmholtz(  
        std::vector<double> & ulon_ulon,
        std::vector<double> & ulon_ulat,
        std::vector<double> & ulat_ulat,
        const std::vector<double> & v_r,
        const std::vector<double> & Phi_v,
        const std::vector<double> & Psi_v,
        const dataset & source_data
    ) {

    const std::vector<double>   &latitude   = source_data.latitude,
                                &longitude  = source_data.longitude;

    const int   Ntime   = source_data.Ntime,    // this is the MPI-local Ntime, not the full Ntime
                Ndepth  = source_data.Ndepth,   // this is the MPI-local Ndepth, not the full Ndepth
                Nlat    = source_data.Nlat,
                Nlon    = source_data.Nlon;

    const size_t Npts = Nlat * Nlon;

    const std::vector<bool> &mask = source_data.mask;

    assert( not(constants::CARTESIAN) );

    alglib::sparsematrix proj_matr;
    alglib::sparsecreate(3*Npts, 3*Npts, proj_matr);

    alglib::real_1d_array Helm_input, uiuj_output;

    std::vector<double> Helm_input_vector(3*Npts), uiuj_output_vector(3*Npts);
    Helm_input.attach_to_ptr(  3*Npts, &Helm_input_vector[0] );
    uiuj_output.attach_to_ptr( 3*Npts, &uiuj_output_vector[0] );

    for (int Itime = 0; Itime < Ntime; ++Itime) {
        for (int Idepth = 0; Idepth < Ndepth; ++Idepth) {

            // Build matrix to get uu from Helmholtz components of v
            build_Helm_matrix( proj_matr, source_data, Itime, Idepth );

            // Copy Helm components into input vector
            for (int Ilat = 0; Ilat < Nlat; ++Ilat) {
                for (int Ilon = 0; Ilon < Nlon; ++Ilon) {
                    size_t index = Index(Itime, Idepth, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon);
                    size_t index_sub = Index(0, 0, Ilat, Ilon, 1, 1, Nlat, Nlon);
                    Helm_input_vector.at( index_sub + 0 * Npts ) = v_r.at( index );
                    Helm_input_vector.at( index_sub + 1 * Npts ) = Phi_v.at( index );
                    Helm_input_vector.at( index_sub + 2 * Npts ) = Psi_v.at( index );
                }
            }

            // Apply proj_matr to Helm components
            alglib::sparsemv( proj_matr, Helm_input, uiuj_output );

            // Copy uiuj results into output arrays
            for (int Ilat = 0; Ilat < Nlat; ++Ilat) {
                for (int Ilon = 0; Ilon < Nlon; ++Ilon) {
                    size_t index = Index(Itime, Idepth, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon);
                    size_t index_sub = Index(0, 0, Ilat, Ilon, 1, 1, Nlat, Nlon);

                    ulon_ulon.at( index ) = uiuj_output_vector.at( index_sub + 0 * Npts );
                    ulon_ulat.at( index ) = uiuj_output_vector.at( index_sub + 1 * Npts );
                    ulat_ulat.at( index ) = uiuj_output_vector.at( index_sub + 2 * Npts );
                }
            }

        }
    }

}

/*
void uiuj_from_Helmholtz(  
        std::vector<double> & ulon_ulon,
        std::vector<double> & ulon_ulat,
        std::vector<double> & ulat_ulat,
        const std::vector<double> & v_r,
        const std::vector<double> & v_lon,
        const std::vector<double> & v_lat,
        const dataset & source_data
    ) {

    const std::vector<double>   &latitude   = source_data.latitude,
                                &longitude  = source_data.longitude;

    const int   Ntime   = source_data.Ntime,    // this is the MPI-local Ntime, not the full Ntime
                Ndepth  = source_data.Ndepth,   // this is the MPI-local Ndepth, not the full Ndepth
                Nlat    = source_data.Nlat,
                Nlon    = source_data.Nlon;

    const std::vector<bool> &mask = source_data.mask;

    assert( not(constants::CARTESIAN) );

    int Itime, Idepth, Ilat, Ilon;
    size_t index;
    double dvlon_dlon, dvlon_dlat, dvlat_dlon, dvlat_dlat, cos_lat, cos2_lat, tan_lat, 
           d2vlon_dlon2, d2vlon_dlat2, d2vlat_dlon2, d2vlat_dlat2,
           tmp_lon, tmp_lat, tmp_v_r, tmp_v_lon, tmp_v_lat, tmp_uu, tmp_uv, tmp_vv;
    std::vector<double*> lon_deriv_vals, lat_deriv_vals, lon2_deriv_vals, lat2_deriv_vals;
    std::vector<const std::vector<double>*> deriv_fields;
    bool is_pole;

    deriv_fields.push_back(&v_lon);
    deriv_fields.push_back(&v_lat);

    #pragma omp parallel \
    default(none) \
    shared( latitude, longitude, mask, v_r, v_lon, v_lat, ulon_ulon, ulon_ulat, ulat_ulat, deriv_fields)\
    private(Itime, Idepth, Ilat, Ilon, index, cos_lat, cos2_lat, tan_lat, tmp_lon, tmp_lat, \
            dvlon_dlon, dvlon_dlat, dvlat_dlon, dvlat_dlat, lon_deriv_vals, lat_deriv_vals, is_pole,\
            d2vlon_dlon2, d2vlon_dlat2, d2vlat_dlon2, d2vlat_dlat2, \
            tmp_uu, tmp_uv, tmp_vv, tmp_v_r, tmp_v_lon, tmp_v_lat )
    {

        lon_deriv_vals.push_back(&dvlon_dlon);
        lon_deriv_vals.push_back(&dvlat_dlon);

        lat_deriv_vals.push_back(&dvlon_dlat);
        lat_deriv_vals.push_back(&dvlat_dlat);

        lon2_deriv_vals.push_back(&d2vlon_dlon2);
        lon2_deriv_vals.push_back(&d2vlat_dlon2);

        lat2_deriv_vals.push_back(&d2vlon_dlat2);
        lat2_deriv_vals.push_back(&d2vlat_dlat2);

        #pragma omp for collapse(1) schedule(guided)
        for (index = 0; index < v_r.size(); ++index) {

            tmp_uu = 0.;
            tmp_uv = 0.;
            tmp_vv = 0.;

            if (mask.at(index)) { // Skip land areas

                Index1to4(index, Itime, Idepth, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon);

                spher_derivative_at_point( lon_deriv_vals, deriv_fields, longitude, "lon",
                        Itime, Idepth, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon, mask);

                spher_derivative_at_point( lat_deriv_vals, deriv_fields, latitude, "lat",
                        Itime, Idepth, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon, mask);

                spher_derivative_at_point( lon2_deriv_vals, deriv_fields, longitude, "lon",
                        Itime, Idepth, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon, mask, 2);

                spher_derivative_at_point( lat2_deriv_vals, deriv_fields, latitude, "lat",
                        Itime, Idepth, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon, mask, 2);

                tmp_v_r   = v_r.at(  index);
                tmp_v_lon = v_lon.at(index);
                tmp_v_lat = v_lat.at(index);

                // If we're too close to the pole (less than 0.01 degrees), bad things happen
                is_pole = std::fabs( std::fabs( latitude.at(Ilat) * 180.0 / M_PI ) - 90 ) < 0.01;
                cos_lat = cos(latitude.at(Ilat));
                cos2_lat = pow( cos(latitude.at(Ilat)), 2);
                tan_lat = tan(latitude.at(Ilat));

                // Ordering is: [ 1        (1/cos(lat)) ddlon          -tan(lat)               ]     [ v_r   ]        [  u_lon * u_lon  ]
                //              [ 0        0.5 ( tan(lat) + ddlat)     (0.5/cos(lat)) ddlon    ]  *  [ v_lon ]   =    [  u_lon * u_lat  ]
                //              [ 1        0                           ddlat                   ]     [ v_lat ]        [  u_lat * u_lat  ]
                tmp_uu = is_pole ? 0. : tmp_v_r + dvlon_dlon / cos_lat - tan_lat * tmp_v_lat;
                tmp_uv = is_pole ? 0. : 0.5 * ( tan_lat * tmp_v_lon + dvlon_dlat + dvlat_dlon / cos_lat  );
                tmp_vv = is_pole ? 0. : tmp_v_r + dvlat_dlat;

            }
            ulon_ulon.at(index) = tmp_uu;
            ulon_ulat.at(index) = tmp_uv;
            ulat_ulat.at(index) = tmp_vv;
        }
    }
}
*/
