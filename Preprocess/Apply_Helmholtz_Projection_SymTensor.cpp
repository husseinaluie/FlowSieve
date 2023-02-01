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


void Apply_Helmholtz_Projection_SymTensor(
        const std::string output_fname,
        dataset & source_data,
        const std::vector<double> & seed_v_r,
        const std::vector<double> & seed_v_lon,
        const std::vector<double> & seed_v_lat,
        const bool single_seed,
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
                                &longitude  = source_data.longitude,
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
    size_t index, index_sub;

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
        full_v_lon( u_lon.size(), 0. ),
        full_v_lat( u_lon.size(), 0. ),
        full_uu(    u_lon.size(), 0. ),
        full_uv(    u_lon.size(), 0. ),
        full_vv(    u_lon.size(), 0. );

    // alglib variables
    alglib::real_1d_array rhs, rhs_seed, rhs_result, lhs_seed;
    std::vector<double> 
        RHS_vector( 3 * Npts, 0.),
        RHS_seed(   3 * Npts, 0.),
        RHS_result( 3 * Npts, 0.),
        LHS_seed(   3 * Npts, 0.);

    rhs.attach_to_ptr(          3 * Npts, &RHS_vector[0] );
    rhs_seed.attach_to_ptr(     3 * Npts, &RHS_seed[0] );
    rhs_result.attach_to_ptr(   3 * Npts, &RHS_result[0] );
    lhs_seed.attach_to_ptr(     3 * Npts, &LHS_seed[0] );
    

    // Copy the starting seed.
    if (single_seed) {
        #pragma omp parallel default(none) shared( LHS_seed, seed_v_r, seed_v_lon, seed_v_lat ) \
        private( Ilat, Ilon, index )
        {
            #pragma omp for collapse(2) schedule(static)
            for (Ilat = 0; Ilat < Nlat; ++Ilat) {
                for (Ilon = 0; Ilon < Nlon; ++Ilon) {
                    index = Index(0, 0, Ilat, Ilon, 1, 1, Nlat, Nlon);
                    LHS_seed.at( index + 0 * Npts ) = seed_v_r.at(   index );
                    LHS_seed.at( index + 1 * Npts ) = seed_v_lon.at( index );
                    LHS_seed.at( index + 2 * Npts ) = seed_v_lat.at( index );
                }
            }
        }
    }

    alglib::linlsqrstate state;
    alglib::linlsqrreport report;

    alglib::real_1d_array F_alglib;

    //
    //// Build the LHS part of the problem
    //      Ordering is: [ 1        (1/cos(lat)) ddlon          -tan(lat)               ]     [ v_r   ]        [  u_lon * u_lon  ]
    //                   [ 0        0.5 ( tan(lat) + ddlat)     (0.5/cos(lat)) ddlon    ]  *  [ v_lon ]   =    [  u_lon * u_lat  ]
    //                   [ 1        0                           ddlat                   ]     [ v_lat ]        [  u_lat * u_lat  ]
    //
    #if DEBUG >= 1
    if (wRank == 0) {
        fprintf(stdout, "Building the LHS of the least squares problem.\n");
        fflush(stdout);
    }
    #endif

    alglib::sparsematrix LHS_matr;
    alglib::sparsecreate(3*Npts, 3*Npts, LHS_matr);

    double tmp_val, weight_val;
    size_t column_skip, row_skip;
    int IDIFF, Idiff, LB, Ndiff;
    std::vector<double> diff_vec;
    for ( Ilat = 0; Ilat < Nlat; Ilat++ ) {
        for ( Ilon = 0; Ilon < Nlon; Ilon++ ) {
            index_sub = Index(0, 0, Ilat, Ilon, 1, 1, Nlat, Nlon);

            weight_val = weight_err ? dAreas.at(index_sub) : 1.;

            // Now handle longitude derivatives
            LB = - 2 * Nlon;
            get_diff_vector(diff_vec, LB, longitude, "lon", Itime, Idepth, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon, mask, 1, constants::DiffOrd);
            Ndiff = diff_vec.size();

            if ( LB != - 2 * Nlon) {
                for ( IDIFF = LB; IDIFF < LB + Ndiff; IDIFF++ ) {
                    if (constants::PERIODIC_X) { Idiff = ( IDIFF % Nlon + Nlon ) % Nlon; }
                    else                       { Idiff = IDIFF;                          }
                    size_t diff_index = Index(0, 0, Ilat, Idiff, 1, 1, Nlat, Nlon);

                    // (0,1) entry
                    row_skip    = 0 * Npts;
                    column_skip = 1 * Npts;
                    tmp_val     = diff_vec.at( IDIFF - LB ) / cos(latitude.at(Ilat));
                    tmp_val    *= weight_val;
                    alglib::sparseadd(  LHS_matr, row_skip + index_sub, column_skip + diff_index, tmp_val);

                    // (1,2) entry
                    row_skip    = 1 * Npts;
                    column_skip = 2 * Npts;
                    tmp_val     = 0.5 * diff_vec.at( IDIFF - LB ) / cos(latitude.at(Ilat));
                    tmp_val    *= weight_val;
                    alglib::sparseadd(  LHS_matr, row_skip + index_sub, column_skip + diff_index, tmp_val);
                }
            }

            // Now handle latitude derivatives
            LB = - 2 * Nlat;
            get_diff_vector(diff_vec, LB, latitude, "lat", Itime, Idepth, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon, mask, 1, constants::DiffOrd);
            Ndiff = diff_vec.size();

            if ( LB != - 2 * Nlat) {
                for ( IDIFF = LB; IDIFF < LB + Ndiff; IDIFF++ ) {
                    if (constants::PERIODIC_Y) { Idiff = ( IDIFF % Nlat + Nlat ) % Nlat; }
                    else                       { Idiff = IDIFF;                          }
                    size_t diff_index = Index(0, 0, Idiff, Ilon, 1, 1, Nlat, Nlon);

                    // (1,1) entry   -  tan(lat) term included later
                    row_skip    = 1 * Npts;
                    column_skip = 1 * Npts;
                    tmp_val     = 0.5 * diff_vec.at( IDIFF - LB );
                    tmp_val    *= weight_val;
                    alglib::sparseadd(  LHS_matr, row_skip + index_sub, column_skip + diff_index, tmp_val);

                    // (2,2) entry
                    row_skip    = 2 * Npts;
                    column_skip = 2 * Npts;
                    tmp_val     = diff_vec.at( IDIFF - LB );
                    tmp_val    *= weight_val;
                    alglib::sparseadd(  LHS_matr, row_skip + index_sub, column_skip + diff_index, tmp_val);
                }
            }

            // (0,0)
            row_skip    = 0 * Npts;
            column_skip = 0 * Npts;
            tmp_val     = 1.;
            tmp_val    *= weight_val;
            alglib::sparseadd(  LHS_matr, row_skip + index_sub, column_skip + index_sub, tmp_val);
            
            // (1,0)
            // This entry is zero
            
            // (2,0)
            row_skip    = 2 * Npts;
            column_skip = 0 * Npts;
            tmp_val     = 1.;
            tmp_val    *= weight_val;
            alglib::sparseadd(  LHS_matr, row_skip + index_sub, column_skip + index_sub, tmp_val);

            // (0,1)
            // Handled by derivatives above

            // (1,1)        -   ddlat term include earlier
            row_skip    = 1 * Npts;
            column_skip = 1 * Npts;
            tmp_val     = 0.5 * tan(latitude.at(Ilat));
            tmp_val    *= weight_val;
            alglib::sparseadd(  LHS_matr, row_skip + index_sub, column_skip + index_sub, tmp_val);
            
            // (2,1)
            // This entry is zero

            // (0,2)
            row_skip    = 0 * Npts;
            column_skip = 2 * Npts;
            tmp_val     = - tan(latitude.at(Ilat));
            tmp_val    *= weight_val;
            alglib::sparseadd(  LHS_matr, row_skip + index_sub, column_skip + index_sub, tmp_val);

            // (1,2)
            // Handled by derivatives above

            // (2,2)
            // Handled by derivatives above

        }
    }

    size_t lower_count = alglib::sparsegetlowercount( LHS_matr );
    size_t upper_count = alglib::sparsegetuppercount( LHS_matr );

    fprintf( stdout, "counts: %'zu %'zu\n", lower_count, upper_count );

    alglib::sparseconverttocrs(LHS_matr);

    #if DEBUG >= 1
    if (wRank == 0) {
        fprintf(stdout, "Declaring the least squares problem.\n");
        fflush(stdout);
    }
    #endif
    alglib::linlsqrcreate(3*Npts, 3*Npts, state);
    alglib::linlsqrsetcond(state, rel_tol, rel_tol, max_iters);

    // Now do the solve!
    for (int Itime = 0; Itime < Ntime; ++Itime) {
        for (int Idepth = 0; Idepth < Ndepth; ++Idepth) {

            if (not(single_seed)) {
                // If single_seed == false, then we were provided seed values, pull out the appropriate values here
                #pragma omp parallel \
                default(none) \
                shared( LHS_seed, seed_v_r, seed_v_lon, seed_v_lat, Itime, Idepth ) \
                private( Ilat, Ilon, index, index_sub )
                {
                    #pragma omp for collapse(2) schedule(static)
                    for (Ilat = 0; Ilat < Nlat; ++Ilat) {
                        for (Ilon = 0; Ilon < Nlon; ++Ilon) {
                            index = Index(Itime, Idepth, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon);
                            index_sub = Index(0, 0, Ilat, Ilon, 1, 1, Nlat, Nlon);
                            LHS_seed.at( index_sub + 0 * Npts ) = seed_v_r.at(   index );
                            LHS_seed.at( index_sub + 1 * Npts ) = seed_v_lon.at( index );
                            LHS_seed.at( index_sub + 2 * Npts ) = seed_v_lat.at( index );
                        }
                    }
                }
            }

            // Get velocity from seed
            alglib::sparsemv( LHS_matr, lhs_seed, rhs_seed );

            #if DEBUG >= 1
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
            
            double weight_val;
            #pragma omp parallel default(none) \
            shared( dAreas, Itime, Idepth, RHS_vector, u_lon, u_lat, RHS_seed ) \
            private( Ilat, Ilon, index, index_sub, weight_val )
            {
                #pragma omp for collapse(2) schedule(static)
                for (Ilat = 0; Ilat < Nlat; ++Ilat) {
                    for (Ilon = 0; Ilon < Nlon; ++Ilon) {
                        index_sub = Index( 0,     0,      Ilat, Ilon, 1,     1,      Nlat, Nlon);
                        index     = Index( Itime, Idepth, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon);

                        weight_val = weight_err ? dAreas.at(index_sub) : 1.;

                        // Seed already scaled by weight from multiplying with LHS_matr
                        RHS_vector.at(0 * Npts + index_sub) = weight_val * u_lon.at(index) * u_lon.at(index) - RHS_seed.at( 0 * Npts + index_sub );
                        RHS_vector.at(1 * Npts + index_sub) = weight_val * u_lon.at(index) * u_lat.at(index) - RHS_seed.at( 1 * Npts + index_sub );
                        RHS_vector.at(2 * Npts + index_sub) = weight_val * u_lat.at(index) * u_lat.at(index) - RHS_seed.at( 2 * Npts + index_sub );
                    }
                }
            }

            //
            //// Now apply the least-squares solver
            //
            #if DEBUG >= 1
            if ( (wRank == 0) and (Itime == 0) ) {
                fprintf(stdout, "Solving the least squares problem.\n");
                fflush(stdout);
            }
            #endif
            alglib::linlsqrsolvesparse(state, LHS_matr, rhs);
            alglib::linlsqrresults(state, F_alglib, report);

            #if DEBUG >= 1
            if      (report.terminationtype == 1) { fprintf(stdout, "Termination type: absolulte tolerance reached.\n"); }
            else if (report.terminationtype == 4) { fprintf(stdout, "Termination type: relative tolerance reached.\n"); }
            else if (report.terminationtype == 5) { fprintf(stdout, "Termination type: maximum number of iterations reached.\n"); }
            else if (report.terminationtype == 7) { fprintf(stdout, "Termination type: round-off errors prevent further progress.\n"); }
            else if (report.terminationtype == 8) { fprintf(stdout, "Termination type: user requested (?)\n"); }
            else                                  { fprintf(stdout, "Termination type: unknown\n"); }
            #endif

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
            if ( (wRank == 0) and (Itime == 0) ) {
                fprintf(stdout, " Done solving the least squares problem.\n");
                fflush(stdout);
            }
            #endif

            // Extract the solution and add the seed back in
            double *LHS_ptr = F_alglib.getcontent();
            std::vector<double> F_array( LHS_ptr, LHS_ptr + 3 * Npts);
            //std::copy ( LHS_ptr, LHS_ptr + 3 * Npts, F_array.begin() );
            for (size_t ii = 0; ii < Npts; ++ii) { F_array.at(ii) += LHS_seed.at(ii); }

            // Get velocity associated to computed F field
            #if DEBUG >= 1
            if ( (wRank == 0) and (Itime == 0) ) {
                fprintf(stdout, " Extracting velocities and divergence from toroidal field.\n");
                fflush(stdout);
            }
            #endif
            alglib::sparsemv( LHS_matr, F_alglib, rhs_result );

            //
            //// Store into the full arrays
            //
            #if DEBUG >= 1
            if ( (wRank == 0) and (Itime == 0) ) {
                fprintf(stdout, " Storing values into output arrays\n");
                fflush(stdout);
            }
            #endif
            #pragma omp parallel \
            default(none) \
            shared( full_v_r, full_v_lon, full_v_lat, full_uu, full_uv, full_vv, \
                    dAreas, F_array, RHS_result, RHS_seed, Itime, Idepth ) \
            private( Ilat, Ilon, index, index_sub, weight_val )
            {
                #pragma omp for collapse(2) schedule(static)
                for (Ilat = 0; Ilat < Nlat; ++Ilat) {
                    for (Ilon = 0; Ilon < Nlon; ++Ilon) {
                        index = Index(Itime, Idepth, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon);

                        index_sub = Index(0, 0, Ilat, Ilon, 1, 1, Nlat, Nlon);

                        full_v_r.at(  index) = F_array.at( index_sub + 0 * Npts );
                        full_v_lon.at(index) = F_array.at( index_sub + 1 * Npts );
                        full_v_lat.at(index) = F_array.at( index_sub + 2 * Npts );

                        full_uu.at( index ) = RHS_result.at( index_sub + 0 * Npts ) + RHS_seed.at( index_sub + 0 * Npts );
                        full_uv.at( index ) = RHS_result.at( index_sub + 1 * Npts ) + RHS_seed.at( index_sub + 1 * Npts );
                        full_vv.at( index ) = RHS_result.at( index_sub + 2 * Npts ) + RHS_seed.at( index_sub + 2 * Npts );

                        // Need to scale weight back out afterwards
                        weight_val = weight_err ? dAreas.at(index_sub) : 1.;

                        full_uu.at( index ) *= 1. / weight_val;
                        full_uv.at( index ) *= 1. / weight_val;
                        full_vv.at( index ) *= 1. / weight_val;
                    }
                }
            }

            // If we don't have a seed for the next iteration, use this solution as the seed
            if (single_seed) { LHS_seed = F_array; }

            #if DEBUG >= 1
            fprintf(stdout, "  --  --  Rank %d done depth %d\n", wRank, Idepth + myStarts.at(1) );
            fflush(stdout);
            #endif

        }

        #if DEBUG >= 1
        fprintf(stdout, " -- Rank %d done time %d\n", wRank, Itime + myStarts.at(0) );
        fflush(stdout);
        #endif
    }

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
    vars_to_write.push_back("v_lon");
    vars_to_write.push_back("v_lat");

    vars_to_write.push_back("uu");
    vars_to_write.push_back("uv");
    vars_to_write.push_back("vv");

    initialize_output_file( source_data, vars_to_write, output_fname.c_str(), -1);

    write_field_to_output(full_v_r,   "v_r",   starts, counts, output_fname.c_str(), &mask);
    write_field_to_output(full_v_lon, "v_lon", starts, counts, output_fname.c_str(), &mask);
    write_field_to_output(full_v_lat, "v_lat", starts, counts, output_fname.c_str(), &mask);

    write_field_to_output(full_uu, "uu", starts, counts, output_fname.c_str(), &mask);
    write_field_to_output(full_uv, "uv", starts, counts, output_fname.c_str(), &mask);
    write_field_to_output(full_vv, "vv", starts, counts, output_fname.c_str(), &mask);

    // Store some solver information
    add_attr_to_file("rel_tol",    rel_tol,                     output_fname.c_str());
    add_attr_to_file("max_iters",  (double) max_iters,          output_fname.c_str());
    add_attr_to_file("diff_order", (double) constants::DiffOrd, output_fname.c_str());
    add_attr_to_file("use_mask",   (double) use_mask,           output_fname.c_str());
    add_attr_to_file("weight_err", (double) weight_err,         output_fname.c_str());


    //
    //// At the very end, compute the L2 error for each time/depth
    //

    std::vector<double> uu_errors( Ntime * Ndepth, 0. ),
                        uv_errors( Ntime * Ndepth, 0. ),
                        vv_errors( Ntime * Ndepth, 0. ),
                        uu_norms(  Ntime * Ndepth, 0. ),
                        uv_norms(  Ntime * Ndepth, 0. ),
                        vv_norms(  Ntime * Ndepth, 0. );
    double total_area = 0, uu_error = 0, uv_error = 0, 
           vv_error = 0, uu_norm = 0, uv_norm = 0, vv_norm = 0;
    for (int Itime = 0; Itime < Ntime; ++Itime) {
        for (int Idepth = 0; Idepth < Ndepth; ++Idepth) {
            #pragma omp parallel \
            default(none) \
            shared( u_lon, u_lat, full_uu, full_uv, full_vv, Itime, Idepth, dAreas ) \
            reduction( + : total_area, uu_error, uv_error, vv_error, uu_norm, uv_norm, vv_norm ) \
            private( Ilat, Ilon, index, index_sub )
            {
                #pragma omp for collapse(2) schedule(static)
                for (Ilat = 0; Ilat < Nlat; ++Ilat) {
                    for (Ilon = 0; Ilon < Nlon; ++Ilon) {
                        index_sub = Index( 0,     0,      Ilat, Ilon, 1,     1,      Nlat, Nlon);
                        index     = Index( Itime, Idepth, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon);

                        total_area += dAreas.at(index_sub);

                        uu_error += dAreas.at(index_sub) * ( pow( u_lon.at(index) * u_lon.at(index) - full_uu.at(index) , 2.) );
                        uv_error += dAreas.at(index_sub) * ( pow( u_lon.at(index) * u_lat.at(index) - full_uv.at(index) , 2.) );
                        vv_error += dAreas.at(index_sub) * ( pow( u_lat.at(index) * u_lat.at(index) - full_vv.at(index) , 2.) );

                        uu_norm += dAreas.at(index_sub) * ( pow( u_lon.at(index) * u_lon.at(index) , 2.) );
                        uv_norm += dAreas.at(index_sub) * ( pow( u_lon.at(index) * u_lat.at(index) , 2.) );
                        vv_norm += dAreas.at(index_sub) * ( pow( u_lat.at(index) * u_lat.at(index) , 2.) );
                    }
                }
            }
            size_t int_index = Index( Itime, Idepth, 0, 0, Ntime, Ndepth, 1, 1);
            uu_errors.at( int_index ) = uu_error / total_area;
            uv_errors.at( int_index ) = uv_error / total_area;
            vv_errors.at( int_index ) = vv_error / total_area;

            uu_norms.at( int_index ) = uu_norm / total_area;
            uv_norms.at( int_index ) = uv_norm / total_area;
            vv_norms.at( int_index ) = vv_norm / total_area;
        }
    }

    const char* dim_names[] = {"time", "depth"};
    const int ndims_error = 2;
    add_var_to_file( "uu_error", dim_names, ndims_error, output_fname.c_str() );
    add_var_to_file( "uv_error", dim_names, ndims_error, output_fname.c_str() );
    add_var_to_file( "vv_error", dim_names, ndims_error, output_fname.c_str() );

    add_var_to_file( "uu_norm", dim_names, ndims_error, output_fname.c_str() );
    add_var_to_file( "uv_norm", dim_names, ndims_error, output_fname.c_str() );
    add_var_to_file( "vv_norm", dim_names, ndims_error, output_fname.c_str() );

    size_t starts_error[ndims_error] = { size_t(myStarts.at(0)), size_t(myStarts.at(1)) };
    size_t counts_error[ndims_error] = { size_t(Ntime), size_t(Ndepth) };

    write_field_to_output( uu_errors, "uu_error", starts_error, counts_error, output_fname.c_str() );
    write_field_to_output( uv_errors, "uv_error", starts_error, counts_error, output_fname.c_str() );
    write_field_to_output( vv_errors, "vv_error", starts_error, counts_error, output_fname.c_str() );

    write_field_to_output( uu_norms, "uu_norm", starts_error, counts_error, output_fname.c_str() );
    write_field_to_output( uv_norms, "uv_norm", starts_error, counts_error, output_fname.c_str() );
    write_field_to_output( vv_norms, "vv_norm", starts_error, counts_error, output_fname.c_str() );

}
