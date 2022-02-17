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

void build_LHS_matrix(
        alglib::sparsematrix & matr,
        const dataset & source_data,
        const std::vector<bool> &mask,
        const int Itime,
        const int Idepth
        ) {

    // Create some tidy names for variables
    const std::vector<double>   &latitude   = source_data.latitude,
                                &longitude  = source_data.longitude,
                                &dAreas     = source_data.areas;

    const std::vector<int>  &myCounts = source_data.myCounts;

    const int   Ntime   = myCounts.at(0),
                Ndepth  = myCounts.at(1),
                Nlat    = myCounts.at(2),
                Nlon    = myCounts.at(3);

    const size_t Npts = Nlat * Nlon;

    double tmp_val, weight_val;
    size_t column_skip, row_skip, index_sub;
    int IDIFF, Idiff, LB, Ndiff;
    std::vector<double> diff_vec;
    for ( int Ilat = 0; Ilat < Nlat; Ilat++ ) {

        for ( int Ilon = 0; Ilon < Nlon; Ilon++ ) {
            index_sub = Index(0, 0, Ilat, Ilon, 1, 1, Nlat, Nlon);

            // (0,0)
            row_skip    = 0 * Npts;
            column_skip = 0 * Npts;
            tmp_val     = 1.;
            alglib::sparseadd(  matr, row_skip + index_sub, column_skip + index_sub, tmp_val);
            
            // (1,1)
            row_skip    = 1 * Npts;
            column_skip = 1 * Npts;
            tmp_val     = 1.;
            alglib::sparseadd(  matr, row_skip + index_sub, column_skip + index_sub, tmp_val);

            // Only using this within 8 degrees of equator, so don't waste time solving the rest of the domain.
            //if ( fabs( latitude.at(Ilat) ) < 8 ) {
            if ( true ) {
                // First latitude derivative
                LB = - 2 * Nlat;
                get_diff_vector(diff_vec, LB, latitude, "lat", Itime, Idepth, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon, mask, 1, constants::DiffOrd);
                Ndiff = diff_vec.size();

                if ( LB != - 2 * Nlat) {
                    for ( IDIFF = LB; IDIFF < LB + Ndiff; IDIFF++ ) {
                        if (constants::PERIODIC_Y) { Idiff = ( IDIFF % Nlat + Nlat ) % Nlat; }
                        else                       { Idiff = IDIFF;                          }
                        size_t diff_index = Index(0, 0, Idiff, Ilon, 1, 1, Nlat, Nlon);

                        // (0,0) entry  
                        row_skip    = 0 * Npts;
                        column_skip = 0 * Npts;
                        tmp_val     = latitude.at(Ilat) * diff_vec.at( IDIFF - LB );
                        alglib::sparseadd(  matr, row_skip + index_sub, column_skip + diff_index, tmp_val);

                        // (1,1) entry
                        row_skip    = 1 * Npts;
                        column_skip = 1 * Npts;
                        tmp_val     = latitude.at(Ilat) * diff_vec.at( IDIFF - LB );
                        alglib::sparseadd(  matr, row_skip + index_sub, column_skip + diff_index, tmp_val);
                    }
                }
            }
        }
    }

    alglib::sparseconverttocrs( matr );
}

void build_RHS_matrix(
        alglib::sparsematrix & matr,
        const dataset & source_data,
        const std::vector<bool> &mask,
        const int Itime,
        const int Idepth
        ) {

    // Create some tidy names for variables
    const std::vector<double>   &latitude   = source_data.latitude,
                                &longitude  = source_data.longitude,
                                &dAreas     = source_data.areas;

    const std::vector<int>  &myCounts = source_data.myCounts;

    const int   Ntime   = myCounts.at(0),
                Ndepth  = myCounts.at(1),
                Nlat    = myCounts.at(2),
                Nlon    = myCounts.at(3);

    const size_t Npts = Nlat * Nlon;
    const double Omega = 2. * M_PI / (24. * 60. * 60.);

    double tmp_val, weight_val;
    size_t column_skip, row_skip, index_sub;
    int IDIFF, Idiff, Idiff_lat, Idiff_lon, 
               LB,    LB_lon,    LB_lat, 
               Ndiff, Ndiff_lon, Ndiff_lat;
    std::vector<double> diff_vec, diff_vec_lon, diff_vec_lat;
    for ( int Ilat = 0; Ilat < Nlat; Ilat++ ) {

        double cos_lat = cos(latitude.at(Ilat));
        double beta_Coriolis = 2 * Omega * cos( latitude.at(Ilat) ) / constants::R_earth;

        for ( int Ilon = 0; Ilon < Nlon; Ilon++ ) {
            index_sub = Index(0, 0, Ilat, Ilon, 1, 1, Nlat, Nlon);

            // Only using this within 5 degrees of equator, so don't waste time solving the rest of the domain.
            //if ( fabs( latitude.at(Ilat) ) < 8 ) {
            if ( true ) {
                // Second latitude derivatives
                LB = - 2 * Nlat;
                get_diff_vector(diff_vec, LB, latitude, "lat", Itime, Idepth, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon, mask, 2, constants::DiffOrd);
                Ndiff = diff_vec.size();

                if ( LB != - 2 * Nlat) {
                    for ( IDIFF = LB; IDIFF < LB + Ndiff; IDIFF++ ) {
                        if (constants::PERIODIC_Y) { Idiff = ( IDIFF % Nlat + Nlat ) % Nlat; }
                        else                       { Idiff = IDIFF;                          }
                        size_t diff_index = Index(0, 0, Idiff, Ilon, 1, 1, Nlat, Nlon);

                        // (0,0) entry  
                        row_skip    = 0 * Npts;
                        column_skip = 0 * Npts;
                        tmp_val     = - diff_vec.at( IDIFF - LB );
                        tmp_val    *= constants::g / beta_Coriolis / pow( constants::R_earth, 2 );
                        alglib::sparseadd(  matr, row_skip + index_sub, column_skip + diff_index, tmp_val);
                    }
                }

                // Mixed partial (first lon and first lat)
                LB_lon = - 2 * Nlon;
                get_diff_vector(diff_vec_lon, LB_lon, longitude, "lon", Itime, Idepth, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon, mask, 1, constants::DiffOrd);
                Ndiff_lon = diff_vec_lon.size();

                LB_lat = - 2 * Nlat;
                get_diff_vector(diff_vec_lat, LB_lat, latitude, "lat", Itime, Idepth, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon, mask, 1, constants::DiffOrd);
                Ndiff_lat = diff_vec_lat.size();

                if ( ( LB_lat != - 2 * Nlat) and ( LB_lon != - 2 * Nlon) ) {
                    for ( int IDIFF_lat = LB_lat; IDIFF_lat < LB_lat + Ndiff_lat; IDIFF_lat++ ) {

                        if (constants::PERIODIC_Y) { Idiff_lat = ( IDIFF_lat % Nlat + Nlat ) % Nlat; }
                        else                       { Idiff_lat = IDIFF_lat;                          }

                        for ( int IDIFF_lon = LB_lon; IDIFF_lon < LB_lon + Ndiff_lon; IDIFF_lon++ ) {

                            if (constants::PERIODIC_X) { Idiff_lon = ( IDIFF_lon % Nlon + Nlon ) % Nlon; }
                            else                       { Idiff_lon = IDIFF_lon;                          }

                            size_t diff_index = Index(0, 0, Idiff_lat, Idiff_lon, 1, 1, Nlat, Nlon);

                            // (1,0) entry  
                            row_skip    = 1 * Npts;
                            column_skip = 0 * Npts;
                            tmp_val     = diff_vec_lon.at( IDIFF_lon - LB_lon ) * diff_vec_lat.at( IDIFF_lat - LB_lat ) / cos_lat;
                            tmp_val    *= constants::g / beta_Coriolis / pow( constants::R_earth, 2 );
                            alglib::sparseadd(  matr, row_skip + index_sub, column_skip + diff_index, tmp_val);
                        }
                    }
                }
            }
        }
    }

    alglib::sparseconverttocrs( matr );

}

void Extract_Beta_Geos_Vel(
        std::vector<double> & u_beta,
        std::vector<double> & v_beta,
        const std::vector<double> & ssh,
        const std::vector<bool> &mask,
        dataset & source_data,
        const double rel_tol,
        const int max_iters,
        const MPI_Comm comm
        ) {

    int wRank, wSize;
    MPI_Comm_rank( comm, &wRank );
    MPI_Comm_size( comm, &wSize );

    // Create some tidy names for variables
    const std::vector<double>   &time       = source_data.time,
                                &depth      = source_data.depth,
                                &latitude   = source_data.latitude,
                                &longitude  = source_data.longitude,
                                &dAreas     = source_data.areas;

    const std::vector<int>  &myCounts = source_data.myCounts,
                            &myStarts = source_data.myStarts;

    const int   Ntime   = myCounts.at(0),
                Ndepth  = myCounts.at(1),
                Nlat    = myCounts.at(2),
                Nlon    = myCounts.at(3);

    const size_t Npts = Nlat * Nlon;
    #if DEBUG >= 1
    if (wRank == 0) { fprintf(stdout, "Working on %'zu points in lat-lon space.\n", Npts); fflush(stdout); }
    #endif

    int Itime, Idepth, Ilat, Ilon;
    size_t index, index_sub;

    // alglib variables
    alglib::real_1d_array ssh_subset, lhs_result, rhs_input;
    std::vector<double> SSH_subset( 1 * Npts, 0.);
    double *LHS_result;

    ssh_subset.attach_to_ptr( 1 * Npts, &SSH_subset[0] );
    
    alglib::linlsqrstate state;
    alglib::linlsqrreport report;

    #if DEBUG >= 1
    if (wRank == 0) { fprintf(stdout, "Building the LHS of the least squares problem.\n"); fflush(stdout); }
    #endif

    alglib::sparsematrix LHS_matr, RHS_matr;
    alglib::sparsecreate(2*Npts, 2*Npts, LHS_matr);
    alglib::sparsecreate(2*Npts, 1*Npts, RHS_matr);

    //build_LHS_matrix( LHS_matr, source_data, mask, Itime, Idepth );
    build_RHS_matrix( RHS_matr, source_data, mask, Itime, Idepth );

    alglib::linlsqrcreate(2*Npts, 2*Npts, state);
    alglib::linlsqrsetcond(state, rel_tol, rel_tol, max_iters);

    // Now do the solve!
    for (int Itime = 0; Itime < Ntime; ++Itime) {
        for (int Idepth = 0; Idepth < Ndepth; ++Idepth) {

            // Get RHS from SSH
            #if DEBUG >= 2
            if ( (wRank == 0) and (Itime == 0) ) { fprintf(stdout, "Preparing the RHS for the problem.\n"); fflush(stdout); }
            #endif
            for (Ilat = 0; Ilat < Nlat; ++Ilat) {
                for (Ilon = 0; Ilon < Nlon; ++Ilon) {
                    index_sub = Index( 0,     0,      Ilat, Ilon, 1,     1,      Nlat, Nlon);
                    index     = Index( Itime, Idepth, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon);
                    SSH_subset.at(index_sub) = ssh.at(index);
                }
            }

            alglib::sparsemv( RHS_matr, ssh_subset, rhs_input );

            // Just take the rhs_input as the result
            //  [ i.e. drop the lat * ddlat term ]
            LHS_result = rhs_input.getcontent();

            /*
            // Now apply the least-squares solver
            #if DEBUG >= 2
            if ( (wRank == 0) and (Itime == 0) ) { fprintf(stdout, "Solving the least squares problem.\n"); fflush(stdout); }
            #endif
            alglib::linlsqrsolvesparse(state, LHS_matr, rhs_input);
            alglib::linlsqrresults(state, lhs_result, report);
            LHS_result = lhs_result.getcontent();
            */

            // Store into the full arrays
            #if DEBUG >= 2
            if ( (wRank == 0) and (Itime == 0) ) { fprintf(stdout, " Storing values into output arrays\n"); fflush(stdout); }
            #endif
            for (Ilat = 0; Ilat < Nlat; ++Ilat) {
                for (Ilon = 0; Ilon < Nlon; ++Ilon) {
                    index = Index(Itime, Idepth, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon);

                    index_sub = Index(0, 0, Ilat, Ilon, 1, 1, Nlat, Nlon);

                    u_beta.at(index) = LHS_result[ index_sub + 0 * Npts ];
                    v_beta.at(index) = LHS_result[ index_sub + 1 * Npts ];
                }
            }
        }
    }
}
