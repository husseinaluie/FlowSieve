#include "../constants.hpp"
#include "../functions.hpp"
#include "../netcdf_io.hpp"
#include "../preprocess.hpp"
#include "../differentiation_tools.hpp"
#include <algorithm>
#include <vector>
#include <omp.h>
#include <math.h>
#include "../ALGLIB/stdafx.h"
#include "../ALGLIB/linalg.h"
#include "../ALGLIB/solvers.h"

void build_main_projection_matrix(
        alglib::sparsematrix & matr,
        const dataset & source_data,
        const int Itime,
        const int Idepth,
        const bool weight_err,
        const double v_r_noise_damp,
        const MPI_Comm comm
        ) {

    int wRank, wSize;
    MPI_Comm_rank( comm, &wRank );
    MPI_Comm_size( comm, &wSize );

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

        double  tan_lat = tan(latitude.at(Ilat)),
                cos_lat = cos(latitude.at(Ilat));

        //if (wRank == 0) {
        //    fprintf( stdout, "lat = %.4g, tan(lat) = %.4g,  cos(lat) = %.4g\n", latitude.at(Ilat), tan_lat, cos_lat );
        //}
        
        for ( int Ilon = 0; Ilon < Nlon; Ilon++ ) {
            index_sub = Index(0, 0, Ilat, Ilon, 1, 1, Nlat, Nlon);

            weight_val = weight_err ? dAreas.at(index_sub) : 1.;


            // Unfortunately, 1 is quite a bit smaller than the magnitude of the second derivative terms.
            //      so terms that are 1 (the v_r terms) get swamped out. So, let's just scale them by a comparable factor
            //      We'll take the mean absolute value of the second latitudinal derivative at the equator, for kicks
            // This same scale factor is then removed from the resulting solution afterwards
            LB = - 2 * Nlat;
            get_diff_vector(diff_vec, LB, latitude, "lat", Itime, Idepth, Nlat/2, 0, Ntime, Ndepth, Nlat, Nlon, unmask, 2, constants::DiffOrd);
            Ndiff = diff_vec.size();
            double Lap_comp_factor = 0;
            for ( IDIFF = 0; IDIFF < (int) diff_vec.size(); IDIFF++ ) { Lap_comp_factor += std::fabs( diff_vec.at(IDIFF) ) / Ndiff; }

            // These are for the v_r term.

            // (0,0)
            row_skip    = 0 * Npts;
            column_skip = 0 * Npts;
            tmp_val     = 1.;
            tmp_val    *= weight_val * Lap_comp_factor;
            alglib::sparseadd(  matr, row_skip + index_sub, column_skip + index_sub, tmp_val);
            
            // (2,0)
            row_skip    = 2 * Npts;
            column_skip = 0 * Npts;
            tmp_val     = 1.;
            tmp_val    *= weight_val * Lap_comp_factor;
            alglib::sparseadd(  matr, row_skip + index_sub, column_skip + index_sub, tmp_val);

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

                    //
                    //// Along southern pole-most latitude, force constant value. This is to eliminate spurious modes from the kernel
                    //
                    if ( (Ilat == 0) or (Ilat == Nlat - 1) ) {
                        row_skip    = 3 * Npts;
                        column_skip = 0 * Npts;
                        tmp_val     = tan_lat * diff_vec.at( IDIFF - LB ) / cos_lat;
                        tmp_val    *= weight_val;
                        alglib::sparseadd(  matr, row_skip + index_sub, column_skip + diff_index, tmp_val);

                        row_skip    = 4 * Npts;
                        column_skip = 1 * Npts;
                        tmp_val     = tan_lat * diff_vec.at( IDIFF - LB ) / cos_lat;
                        tmp_val    *= weight_val;
                        alglib::sparseadd(  matr, row_skip + index_sub, column_skip + diff_index, tmp_val);

                        row_skip    = 5 * Npts;
                        column_skip = 2 * Npts;
                        tmp_val     = tan_lat * diff_vec.at( IDIFF - LB ) / cos_lat;
                        tmp_val    *= weight_val;
                        //alglib::sparseadd(  matr, row_skip + index_sub, column_skip + diff_index, tmp_val);
                        alglib::sparseset(  matr, row_skip + index_sub, column_skip + index_sub, Lap_comp_factor);
                    }
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

                    // Try to remove jagged noise (Laplace weight)
                    if ( not( (Ilat == 0) or (Ilat == Nlat - 1) ) ) {
                        for ( int II = 0; II < 3; II++ ) {
                            row_skip    = (3+II) * Npts;
                            column_skip = II * Npts;
                            tmp_val     = v_r_noise_damp * diff_vec.at( IDIFF - LB ) / pow(cos_lat, 2.);
                            tmp_val    *= weight_val;
                            alglib::sparseadd(  matr, row_skip + index_sub, column_skip + diff_index, tmp_val);
                        }
                    }
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

                    // Try to remove jagged noise from (Laplace weight)
                    if ( not( (Ilat == 0) or (Ilat == Nlat - 1) ) ) {
                        for ( int II = 0; II < 3; II++ ) {
                            row_skip    = (3+II) * Npts;
                            column_skip = II * Npts;
                            tmp_val     = - v_r_noise_damp * diff_vec.at( IDIFF - LB ) * tan_lat;
                            tmp_val    *= weight_val;
                            alglib::sparseadd(  matr, row_skip + index_sub, column_skip + diff_index, tmp_val);
                        }
                    }
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

                    // Try to remove jagged noise (Laplace weight)
                    if ( not( (Ilat == 0) or (Ilat == Nlat - 1) ) ) {
                        for ( int II = 0; II < 3; II++ ) {
                            row_skip    = (3+II) * Npts;
                            column_skip = II * Npts;
                            tmp_val     = v_r_noise_damp * diff_vec.at( IDIFF - LB );
                            tmp_val    *= weight_val;
                            alglib::sparseadd(  matr, row_skip + index_sub, column_skip + diff_index, tmp_val);
                        }
                    }
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

        }
    }

    //size_t lower_count = alglib::sparsegetlowercount( matr );
    //size_t upper_count = alglib::sparsegetuppercount( matr );

    alglib::sparseconverttocrs( matr );

}

void Apply_Helmholtz_Projection_uiuj(
        const std::string output_fname,
        dataset & source_data,
        const std::vector<double> & seed_v_r,
        const std::vector<double> & seed_Phi_v,
        const std::vector<double> & seed_Psi_v,
        const bool single_seed,
        const double Tikhov_Lambda,
        const double Tikhov_Laplace,
        const double rel_tol,
        const int max_iters,
        const bool weight_err,
        const bool use_mask,
        const MPI_Comm comm
        ) {

    int wRank, wSize;
    MPI_Comm_rank( comm, &wRank );
    MPI_Comm_size( comm, &wSize );

    // Create some tidy names for variables
    const std::vector<double>   &latitude   = source_data.latitude,
                                &dAreas     = source_data.areas;

    const std::vector<bool> &mask = source_data.mask;

    const std::vector<int>  &myCounts = source_data.myCounts,
                            &myStarts = source_data.myStarts;

    std::vector<double>   &u_lat = source_data.variables.at("u_lat"),
                          &u_lon = source_data.variables.at("u_lon");

    // Create a 'no mask' mask variable
    //   we'll treat land values as zero velocity
    //   We do this because including land seems
    //   to introduce strong numerical issues
    const std::vector<bool> unmask(mask.size(), true);

    const int   Ntime   = myCounts.at(0),
                Ndepth  = myCounts.at(1),
                Nlat    = myCounts.at(2),
                Nlon    = myCounts.at(3);

    const size_t Npts = Nlat * Nlon;

    int Itime = 0, Idepth = 0, Ilat, Ilon;
    size_t index, index_sub, iters_used;

    // Fill in the land areas with zero velocity
    #pragma omp parallel default(none) shared( u_lon, u_lat, mask ) private( index )
    {
        #pragma omp for collapse(1) schedule(guided)
        for (index = 0; index < u_lon.size(); index++) {
            if (not(mask.at(index))) {
                u_lon.at(index) = 0.;
                u_lat.at(index) = 0.;
            }
        }
    }

    // Storage vectors
    std::vector<double> 
        full_v_r(   u_lon.size(), 0. ),
        full_Phi_v( u_lon.size(), 0. ),
        full_Psi_v( u_lon.size(), 0. ),
        full_uu(    u_lon.size(), 0. ),
        full_uv(    u_lon.size(), 0. ),
        full_vv(    u_lon.size(), 0. );

    // alglib variables
    alglib::real_1d_array rhs, rhs_seed, rhs_result, lhs_seed;
    std::vector<double> 
        RHS_vector( 6 * Npts, 0.),
        RHS_seed(   6 * Npts, 0.),
        LHS_seed(   3 * Npts, 0.),
        RHS_result( 3 * Npts, 0.);

    rhs.attach_to_ptr(              6 * Npts, &RHS_vector[0] );
    rhs_seed.attach_to_ptr(         6 * Npts, &RHS_seed[0] );
    lhs_seed.attach_to_ptr(         3 * Npts, &LHS_seed[0] );
    rhs_result.attach_to_ptr(       3 * Npts, &RHS_result[0] );
    
    // Laplace comparison scale factor, for v_r
    int LB = - 2 * Nlat;
    std::vector<double> diff_vec;
    get_diff_vector(diff_vec, LB, latitude, "lat", Itime, Idepth, Nlat/2, 0, Ntime, Ndepth, Nlat, Nlon, unmask, 2, constants::DiffOrd);
    int Ndiff = diff_vec.size();
    double Lap_comp_factor = 0;
    for ( size_t IDIFF = 0; IDIFF < diff_vec.size(); IDIFF++ ) { 
        Lap_comp_factor += std::fabs( diff_vec.at(IDIFF) ) / Ndiff; 
    }
    if (wRank == 0) { fprintf( stdout, "Lap_comp_factor = %g\n", Lap_comp_factor ); }

    // Copy the starting seed.
    if (single_seed) {
        #pragma omp parallel default(none) \
        shared( LHS_seed, seed_v_r, seed_Phi_v, seed_Psi_v, Lap_comp_factor ) \
        private( Ilat, Ilon, index ) \
        firstprivate( Nlon, Nlat, Npts )
        {
            #pragma omp for collapse(2) schedule(static)
            for (Ilat = 0; Ilat < Nlat; ++Ilat) {
                for (Ilon = 0; Ilon < Nlon; ++Ilon) {
                    index = Index(0, 0, Ilat, Ilon, 1, 1, Nlat, Nlon);
                    LHS_seed.at( index + 0 * Npts ) = seed_v_r.at( index ) / Lap_comp_factor;
                    LHS_seed.at( index + 1 * Npts ) = seed_Phi_v.at( index );
                    LHS_seed.at( index + 2 * Npts ) = seed_Psi_v.at( index );
                }
            }
        }
    }

    alglib::linlsqrstate state;
    alglib::linlsqrreport report;

    alglib::real_1d_array F_alglib;

    //
    //// Build the LHS part of the problem
    //      Ordering is: [ 1       (1/cos(lat))^2 d^2dlon^2 - tan(lat) ddlat        - (1/cos(lat))( tan(lat) ddlon + d^2dlatdlon )                      ]     [ v_r   ]        [  u_lon * u_lon  ]
    //                   [ 0       (1/cos(lat))( tan(lat) ddon + d^2dlatdlon )      0.5 * ( -tan(lat)ddlat - d^2ddlat^2 + (1/cos(lat))^2 d^2ddlon^2 )   ]  *  [ Phi_v ]   =    [  u_lon * u_lat  ]
    //                   [ 1       d^2ddlat^2                                       (1/cos(lat))( tan(lat) ddon + d^2dlatdlon )                         ]     [ Psi_v ]        [  u_lat * u_lat  ]
    //
    //      Ordering is: [ (1/cos(lat))^2 d^2dlon^2 - tan(lat) ddlat - d^2ddlat^2        - (2/cos(lat))( tan(lat) ddlon + d^2dlatdlon )                      ]     [ Phi_v ]        [  u_lon * u_lon  - ulat * ulat ]
    //                   [ (1/cos(lat))( tan(lat) ddon + d^2dlatdlon )                   0.5 * ( -tan(lat)ddlat - d^2ddlat^2 + (1/cos(lat))^2 d^2ddlon^2 )   ]  *  [ Psi_v ]   =    [  u_lon * u_lat  ]
    //
    #if DEBUG >= 1
    if (wRank == 0) {
        fprintf(stdout, "Building the LHS of the least squares problem.\n");
        fflush(stdout);
    }
    #endif

    alglib::sparsematrix proj_matr;
    alglib::sparsecreate(6*Npts, 3*Npts, proj_matr);

    const double v_r_noise_damp = Tikhov_Laplace; // 0.05
    build_main_projection_matrix(    proj_matr,   source_data, Itime, Idepth, weight_err, v_r_noise_damp, comm);

    #if DEBUG >= 1
    if (wRank == 0) {
        fprintf(stdout, "Declaring the least squares problem.\n");
        fflush(stdout);
    }
    #endif
    alglib::linlsqrcreate(6*Npts, 3*Npts, state);
    alglib::linlsqrsetcond(state, rel_tol, rel_tol, max_iters);
    alglib::linlsqrsetlambdai( state, Tikhov_Lambda );

    // Counters to track termination types
    int terminate_count_abs_tol = 0,
        terminate_count_rel_tol = 0,
        terminate_count_max_iter = 0,
        terminate_count_rounding = 0,
        terminate_count_other = 0;

    // Now do the solve!
    for (int Itime = 0; Itime < Ntime; ++Itime) {
        for (int Idepth = 0; Idepth < Ndepth; ++Idepth) {

            if (not(single_seed)) {
                // If single_seed == false, then we were provided seed values, pull out the appropriate values here
                #pragma omp parallel \
                default(none) \
                shared( LHS_seed, seed_v_r, seed_Phi_v, seed_Psi_v, Itime, Idepth, Lap_comp_factor ) \
                private( Ilat, Ilon, index, index_sub ) \
                firstprivate( Nlon, Nlat, Ndepth, Ntime, Npts )
                {
                    #pragma omp for collapse(2) schedule(static)
                    for (Ilat = 0; Ilat < Nlat; ++Ilat) {
                        for (Ilon = 0; Ilon < Nlon; ++Ilon) {
                            index = Index(Itime, Idepth, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon);
                            index_sub = Index(0, 0, Ilat, Ilon, 1, 1, Nlat, Nlon);
                            LHS_seed.at( index_sub + 0 * Npts ) = seed_v_r.at( index ) / Lap_comp_factor;
                            LHS_seed.at( index_sub + 1 * Npts ) = seed_Phi_v.at( index );
                            LHS_seed.at( index_sub + 2 * Npts ) = seed_Psi_v.at( index );
                        }
                    }
                }
            }

            #if DEBUG >= 2
            if (wRank == 0) {
                fprintf(stdout, "Computing RHS-equivalent of seed.\n");
                fflush(stdout);
            }
            #endif

            // Get velocity from seed
            alglib::sparsemv( proj_matr, lhs_seed, rhs_seed );

            #if DEBUG >= 2
            if ( (wRank == 0) and (Itime == 0) ) {
                fprintf(stdout, "Building the RHS of the least squares problem.\n");
                fflush(stdout);
            }
            #endif

            //
            //// Set up the RHS_vector
            //
            //  Subtract off the velocity from the seed
            //
            
            double weight_val, u_lon_loc, u_lat_loc;
            #pragma omp parallel default(none) \
            shared( dAreas, Itime, Idepth, RHS_vector, u_lon, u_lat, RHS_seed ) \
            private( Ilat, Ilon, index, index_sub, weight_val, u_lon_loc, u_lat_loc ) \
            firstprivate( Nlon, Nlat, Ndepth, Ntime, Npts, weight_err )
            {
                #pragma omp for collapse(2) schedule(static)
                for (Ilat = 0; Ilat < Nlat; ++Ilat) {
                    for (Ilon = 0; Ilon < Nlon; ++Ilon) {
                        index_sub = Index( 0,     0,      Ilat, Ilon, 1,     1,      Nlat, Nlon);
                        index     = Index( Itime, Idepth, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon);

                        u_lon_loc = u_lon.at(index);
                        u_lat_loc = u_lat.at(index);
                        RHS_vector.at(0 * Npts + index_sub) = pow( u_lon_loc, 2. );
                        RHS_vector.at(1 * Npts + index_sub) = u_lon_loc * u_lat_loc;
                        RHS_vector.at(2 * Npts + index_sub) = pow( u_lat_loc, 2. );

                        weight_val = weight_err ? dAreas.at(index_sub) : 1.;
                        RHS_vector.at(0 * Npts + index_sub) *= weight_val; 
                        RHS_vector.at(1 * Npts + index_sub) *= weight_val; 
                        RHS_vector.at(2 * Npts + index_sub) *= weight_val; 

                        // Seed already scaled by weight from multiplying with LHS_matr
                        RHS_vector.at(0 * Npts + index_sub) += -RHS_seed.at( 0 * Npts + index_sub );
                        RHS_vector.at(1 * Npts + index_sub) += -RHS_seed.at( 1 * Npts + index_sub );
                        RHS_vector.at(2 * Npts + index_sub) += -RHS_seed.at( 2 * Npts + index_sub );
                    }
                }
            }

            //
            //// Now apply the least-squares solver
            //
            #if DEBUG >= 2
            if ( (wRank == 0) and (Itime == 0) ) {
                fprintf(stdout, "Solving the least squares problem.\n");
                fflush(stdout);
            }
            #endif
            alglib::linlsqrsolvesparse(state, proj_matr, rhs);
            alglib::linlsqrresults(state, F_alglib, report);

            /*    Rep     -   optimization report:
                * Rep.TerminationType completetion code:
                    *  1    ||Rk||<=EpsB*||B||
                    *  4    ||A^T*Rk||/(||A||*||Rk||)<=EpsA
                    *  5    MaxIts steps was taken
                    *  7    rounding errors prevent further progress,
                            X contains best point found so far.
                            (sometimes returned on singular systems)
                    *  8    user requested termination via calling
                            linlsqrrequesttermination()
                * Rep.IterationsCount contains iterations count
                * NMV countains number of matrix-vector calculations
            */

            #if DEBUG >= 1
            if      (report.terminationtype == 1) { fprintf(stdout, "Termination type: absolulte tolerance reached.\n"); }
            else if (report.terminationtype == 4) { fprintf(stdout, "Termination type: relative tolerance reached.\n"); }
            else if (report.terminationtype == 5) { fprintf(stdout, "Termination type: maximum number of iterations reached.\n"); }
            else if (report.terminationtype == 7) { fprintf(stdout, "Termination type: round-off errors prevent further progress.\n"); }
            else if (report.terminationtype == 8) { fprintf(stdout, "Termination type: user requested (?)\n"); }
            else                                  { fprintf(stdout, "Termination type: unknown\n"); }
            #endif
            if      (report.terminationtype == 1) { terminate_count_abs_tol++; }
            else if (report.terminationtype == 4) { terminate_count_rel_tol++; }
            else if (report.terminationtype == 5) { terminate_count_max_iter++; }
            else if (report.terminationtype == 7) { terminate_count_rounding++; }
            else if (report.terminationtype == 8) { terminate_count_other++; }
            else                                  { terminate_count_other++; }

            iters_used = linlsqrpeekiterationscount( state );

            #if DEBUG >= 2
            if ( (wRank == 0) and (Itime == 0) ) {
                fprintf(stdout, " Done solving the least squares problem.\n");
                fflush(stdout);
            }
            #endif

            // Add the seed back in to the solution
            double *LHS_ptr = F_alglib.getcontent();
            for (size_t ii = 0; ii < 3 * Npts; ++ii) { LHS_ptr[ii] += LHS_seed.at(ii); }

            // Get velocity associated to computed F field
            #if DEBUG >= 2
            if ( (wRank == 0) and (Itime == 0) ) {
                fprintf(stdout, " Extracting uiuj from full Helmholtz vector.\n");
                fflush(stdout);
            }
            #endif
            alglib::real_1d_array lhs_result;
            lhs_result.attach_to_ptr( 3 * Npts, &LHS_ptr[0] );
            alglib::sparsemv( proj_matr, lhs_result, rhs_result );
            double *RHS_result_ptr = rhs_result.getcontent();

            //
            //// Store into the full arrays
            //
            #if DEBUG >= 2
            if ( (wRank == 0) and (Itime == 0) ) {
                fprintf(stdout, " Storing values into output arrays\n");
                fflush(stdout);
            }
            #endif
            #pragma omp parallel \
            default(none) \
            shared( full_v_r, full_Phi_v, full_Psi_v, full_uu, full_uv, full_vv, \
                    dAreas, LHS_ptr, RHS_result, RHS_result_ptr, Itime, Idepth, Lap_comp_factor ) \
            private( Ilat, Ilon, index, index_sub, weight_val ) \
            firstprivate( Nlon, Nlat, Ndepth, Ntime, Npts, weight_err )
            {
                #pragma omp for collapse(2) schedule(static)
                for (Ilat = 0; Ilat < Nlat; ++Ilat) {
                    for (Ilon = 0; Ilon < Nlon; ++Ilon) {
                        index = Index(Itime, Idepth, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon);

                        index_sub = Index(0, 0, Ilat, Ilon, 1, 1, Nlat, Nlon);

                        full_v_r.at(  index) = LHS_ptr[ index_sub + 0 * Npts ] * Lap_comp_factor;
                        full_Phi_v.at(index) = LHS_ptr[ index_sub + 1 * Npts ];
                        full_Psi_v.at(index) = LHS_ptr[ index_sub + 2 * Npts ];

                        weight_val = weight_err ? dAreas.at(index_sub) : 1.;
                        //full_uu.at( index ) = RHS_result.at( index_sub + 0 * Npts );
                        //full_uv.at( index ) = RHS_result.at( index_sub + 1 * Npts );
                        //full_vv.at( index ) = RHS_result.at( index_sub + 2 * Npts );
                        full_uu.at( index ) = RHS_result_ptr[ index_sub + 0 * Npts ] / weight_val;
                        full_uv.at( index ) = RHS_result_ptr[ index_sub + 1 * Npts ] / weight_val;
                        full_vv.at( index ) = RHS_result_ptr[ index_sub + 2 * Npts ] / weight_val;
                    }
                }
            }

            // If we don't have a seed for the next iteration, use this solution as the seed
            if (single_seed) { for (size_t ii = 0; ii < 3 * Npts; ++ii) { LHS_seed.at(ii) = LHS_ptr[ii]; } }

            #if DEBUG >= 0
            if ( source_data.full_Ndepth > 1 ) {
                fprintf(stdout, "  --  --  Rank %d done depth %d after %'zu iterations\n", wRank, Idepth + myStarts.at(1), iters_used );
                fflush(stdout);
            }
            #endif

        }

        #if DEBUG >= 0
        if ( source_data.full_Ntime > 1 ) {
            fprintf(stdout, " -- Rank %d done time %d after %'zu iterations\n", wRank, Itime + myStarts.at(0), iters_used );
            fflush(stdout);
        }
        #endif
    }

    //
    //// Print termination counts
    //

    int total_count_abs_tol, total_count_rel_tol, total_count_max_iter, total_count_rounding, total_count_other;

    MPI_Reduce( &terminate_count_abs_tol,  &total_count_abs_tol,  1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD );
    MPI_Reduce( &terminate_count_rel_tol,  &total_count_rel_tol,  1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD );
    MPI_Reduce( &terminate_count_max_iter, &total_count_max_iter, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD );
    MPI_Reduce( &terminate_count_rounding, &total_count_rounding, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD );
    MPI_Reduce( &terminate_count_other,    &total_count_other,    1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD );

    #if DEBUG >= 0
    if (wRank == 0) {
        fprintf( stdout, "\n" );
        fprintf( stdout, "Termination counts: %'d from absolute tolerance\n", total_count_abs_tol );
        fprintf( stdout, "                    %'d from relative tolerance\n", total_count_rel_tol );
        fprintf( stdout, "                    %'d from iteration maximum\n", total_count_max_iter );
        fprintf( stdout, "                    %'d from rounding errors \n", total_count_rounding );
        fprintf( stdout, "                    %'d from other causes \n", total_count_other );
        fprintf( stdout, "\n" );
    }
    #endif
    if      (report.terminationtype == 1) { terminate_count_abs_tol++; }
    else if (report.terminationtype == 4) { terminate_count_rel_tol++; }
    else if (report.terminationtype == 5) { terminate_count_max_iter++; }
    else if (report.terminationtype == 7) { terminate_count_rounding++; }
    else if (report.terminationtype == 8) { terminate_count_other++; }
    else                                  { terminate_count_other++; }


    //
    //// Write the output
    //

    const int ndims = 4;
    size_t starts[ndims] = {
        size_t(myStarts.at(0)), size_t(myStarts.at(1)), 
        size_t(myStarts.at(2)), size_t(myStarts.at(3))
    };
    size_t counts[ndims] = {
        size_t(Ntime), size_t(Ndepth), 
        size_t(Nlat),  size_t(Nlon)
    };

    std::vector<std::string> vars_to_write;
    vars_to_write.push_back("v_r");
    vars_to_write.push_back("Phi_v");
    vars_to_write.push_back("Psi_v");

    if (not(constants::MINIMAL_OUTPUT)) {
        vars_to_write.push_back("uu");
        vars_to_write.push_back("uv");
        vars_to_write.push_back("vv");
    }

    initialize_output_file( source_data, vars_to_write, output_fname.c_str(), -1);

    write_field_to_output(full_v_r,   "v_r",   starts, counts, output_fname.c_str(), &unmask);
    write_field_to_output(full_Phi_v, "Phi_v", starts, counts, output_fname.c_str(), &unmask);
    write_field_to_output(full_Psi_v, "Psi_v", starts, counts, output_fname.c_str(), &unmask);

    if (not(constants::MINIMAL_OUTPUT)) {
        write_field_to_output(full_uu, "uu", starts, counts, output_fname.c_str(), &unmask);
        write_field_to_output(full_uv, "uv", starts, counts, output_fname.c_str(), &unmask);
        write_field_to_output(full_vv, "vv", starts, counts, output_fname.c_str(), &unmask);
    }

    // Store some solver information
    add_attr_to_file("rel_tol",         rel_tol,                        output_fname.c_str());
    add_attr_to_file("max_iters",       (double) max_iters,             output_fname.c_str());
    add_attr_to_file("diff_order",      (double) constants::DiffOrd,    output_fname.c_str());
    add_attr_to_file("use_mask",        (double) use_mask,              output_fname.c_str());
    add_attr_to_file("weight_err",      (double) weight_err,            output_fname.c_str());
    add_attr_to_file("Tikhov_Lambda",   Tikhov_Lambda,                  output_fname.c_str());
    add_attr_to_file("Tikhov_Laplace",  Tikhov_Laplace,                 output_fname.c_str());


    //
    //// At the very end, compute the L2 and Linf error for each time/depth
    //

    std::vector<double> uu_2errors( Ntime * Ndepth, 0. ),
                        uv_2errors( Ntime * Ndepth, 0. ),
                        vv_2errors( Ntime * Ndepth, 0. ),
                        uu_2norms(  Ntime * Ndepth, 0. ),
                        uv_2norms(  Ntime * Ndepth, 0. ),
                        vv_2norms(  Ntime * Ndepth, 0. ),
                        uu_Inferrors( Ntime * Ndepth, 0. ),
                        uv_Inferrors( Ntime * Ndepth, 0. ),
                        vv_Inferrors( Ntime * Ndepth, 0. ),
                        uu_Infnorms(  Ntime * Ndepth, 0. ),
                        uv_Infnorms(  Ntime * Ndepth, 0. ),
                        vv_Infnorms(  Ntime * Ndepth, 0. );
    double total_area, uu_2error,   uv_2error,   vv_2error,   uu_2norm,   uv_2norm,   vv_2norm,
                       uu_Inferror, uv_Inferror, vv_Inferror, uu_Infnorm, uv_Infnorm, vv_Infnorm;
    for (int Itime = 0; Itime < Ntime; ++Itime) {
        for (int Idepth = 0; Idepth < Ndepth; ++Idepth) {

            total_area = 0.;
            uu_2error  = 0.;
            uv_2error  = 0.;
            vv_2error  = 0.;
            uu_2norm   = 0.;
            uv_2norm   = 0.;
            vv_2norm   = 0.;

            uu_Inferror = 0.;
            uv_Inferror = 0.;
            vv_Inferror = 0.;
            uu_Infnorm  = 0.;
            uv_Infnorm  = 0.;
            vv_Infnorm  = 0.;

            #pragma omp parallel \
            default(none) \
            shared( u_lon, u_lat, full_uu, full_uv, full_vv, Itime, Idepth, dAreas ) \
            reduction( + : total_area, uu_2error, uv_2error, vv_2error, uu_2norm, uv_2norm, vv_2norm ) \
            reduction( max : uu_Inferror, uv_Inferror, vv_Inferror, uu_Infnorm, uv_Infnorm, vv_Infnorm ) \
            private( Ilat, Ilon, index, index_sub ) \
            firstprivate( Nlon, Nlat, Ndepth, Ntime )
            {
                #pragma omp for collapse(2) schedule(static)
                for (Ilat = 0; Ilat < Nlat; ++Ilat) {
                    for (Ilon = 0; Ilon < Nlon; ++Ilon) {
                        index_sub = Index( 0,     0,      Ilat, Ilon, 1,     1,      Nlat, Nlon);
                        index     = Index( Itime, Idepth, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon);

                        total_area += dAreas.at(index_sub);

                        uu_2error += dAreas.at(index_sub) * ( pow( u_lon.at(index) * u_lon.at(index) - full_uu.at(index) , 2.) );
                        uv_2error += dAreas.at(index_sub) * ( pow( u_lon.at(index) * u_lat.at(index) - full_uv.at(index) , 2.) );
                        vv_2error += dAreas.at(index_sub) * ( pow( u_lat.at(index) * u_lat.at(index) - full_vv.at(index) , 2.) );

                        uu_2norm += dAreas.at(index_sub) * ( pow( u_lon.at(index) * u_lon.at(index) , 2.) );
                        uv_2norm += dAreas.at(index_sub) * ( pow( u_lon.at(index) * u_lat.at(index) , 2.) );
                        vv_2norm += dAreas.at(index_sub) * ( pow( u_lat.at(index) * u_lat.at(index) , 2.) );

                        uu_Inferror = std::fmax( std::fabs( u_lon.at(index) * u_lon.at(index) - full_uu.at(index) ), uu_Inferror );
                        uv_Inferror = std::fmax( std::fabs( u_lon.at(index) * u_lat.at(index) - full_uv.at(index) ), uv_Inferror );
                        vv_Inferror = std::fmax( std::fabs( u_lat.at(index) * u_lat.at(index) - full_vv.at(index) ), vv_Inferror );

                        uu_Infnorm = std::fmax( std::fabs( u_lon.at(index) * u_lon.at(index) ), uu_Infnorm );
                        uv_Infnorm = std::fmax( std::fabs( u_lon.at(index) * u_lat.at(index) ), uv_Infnorm );
                        vv_Infnorm = std::fmax( std::fabs( u_lat.at(index) * u_lat.at(index) ), vv_Infnorm );
                    }
                }
            }
            size_t int_index = Index( Itime, Idepth, 0, 0, Ntime, Ndepth, 1, 1);
            uu_2errors.at( int_index ) = sqrt( uu_2error / total_area );
            uv_2errors.at( int_index ) = sqrt( uv_2error / total_area );
            vv_2errors.at( int_index ) = sqrt( vv_2error / total_area );

            uu_2norms.at( int_index ) = sqrt( uu_2norm / total_area );
            uv_2norms.at( int_index ) = sqrt( uv_2norm / total_area );
            vv_2norms.at( int_index ) = sqrt( vv_2norm / total_area );

            uu_Inferrors.at( int_index ) = uu_Inferror;
            uv_Inferrors.at( int_index ) = uv_Inferror;
            vv_Inferrors.at( int_index ) = vv_Inferror;

            uu_Infnorms.at( int_index ) = uu_Infnorm;
            uv_Infnorms.at( int_index ) = uv_Infnorm;
            vv_Infnorms.at( int_index ) = vv_Infnorm;
        }
    }

    const char* dim_names[] = {"time", "depth"};
    const int ndims_error = 2;
    if (wRank == 0) {
        add_var_to_file( "uu_2error", dim_names, ndims_error, output_fname.c_str() );
        add_var_to_file( "uv_2error", dim_names, ndims_error, output_fname.c_str() );
        add_var_to_file( "vv_2error", dim_names, ndims_error, output_fname.c_str() );

        add_var_to_file( "uu_2norm", dim_names, ndims_error, output_fname.c_str() );
        add_var_to_file( "uv_2norm", dim_names, ndims_error, output_fname.c_str() );
        add_var_to_file( "vv_2norm", dim_names, ndims_error, output_fname.c_str() );

        add_var_to_file( "uu_Inferror", dim_names, ndims_error, output_fname.c_str() );
        add_var_to_file( "uv_Inferror", dim_names, ndims_error, output_fname.c_str() );
        add_var_to_file( "vv_Inferror", dim_names, ndims_error, output_fname.c_str() );

        add_var_to_file( "uu_Infnorm", dim_names, ndims_error, output_fname.c_str() );
        add_var_to_file( "uv_Infnorm", dim_names, ndims_error, output_fname.c_str() );
        add_var_to_file( "vv_Infnorm", dim_names, ndims_error, output_fname.c_str() );
    }
    MPI_Barrier(MPI_COMM_WORLD);

    size_t starts_error[ndims_error] = { size_t(myStarts.at(0)), size_t(myStarts.at(1)) };
    size_t counts_error[ndims_error] = { size_t(Ntime), size_t(Ndepth) };

    write_field_to_output( uu_2errors, "uu_2error", starts_error, counts_error, output_fname.c_str() );
    write_field_to_output( uv_2errors, "uv_2error", starts_error, counts_error, output_fname.c_str() );
    write_field_to_output( vv_2errors, "vv_2error", starts_error, counts_error, output_fname.c_str() );

    write_field_to_output( uu_2norms, "uu_2norm", starts_error, counts_error, output_fname.c_str() );
    write_field_to_output( uv_2norms, "uv_2norm", starts_error, counts_error, output_fname.c_str() );
    write_field_to_output( vv_2norms, "vv_2norm", starts_error, counts_error, output_fname.c_str() );

    write_field_to_output( uu_Inferrors, "uu_Inferror", starts_error, counts_error, output_fname.c_str() );
    write_field_to_output( uv_Inferrors, "uv_Inferror", starts_error, counts_error, output_fname.c_str() );
    write_field_to_output( vv_Inferrors, "vv_Inferror", starts_error, counts_error, output_fname.c_str() );

    write_field_to_output( uu_Infnorms, "uu_Infnorm", starts_error, counts_error, output_fname.c_str() );
    write_field_to_output( uv_Infnorms, "uv_Infnorm", starts_error, counts_error, output_fname.c_str() );
    write_field_to_output( vv_Infnorms, "vv_Infnorm", starts_error, counts_error, output_fname.c_str() );

}
