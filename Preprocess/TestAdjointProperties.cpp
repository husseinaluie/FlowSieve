#include "../constants.hpp"
#include "../functions.hpp"
#include "../netcdf_io.hpp"
#include "../preprocess.hpp"
#include "../differentiation_tools.hpp"
#include <algorithm>
#include <vector>
#include <omp.h>
#include <math.h>
#include <numeric>   // std::accumulate
#include <random>


void TestAdjointProperties( const dataset & grid ) {

    
    const int Nlat = grid.Nlat;
    const int Nlon = grid.Nlon;
    const size_t Npts = ((size_t) Nlat) * ( (size_t) Nlon);

    std::vector<short int> mask( Npts, true );

    size_t Ipt;

    std::random_device rd;  // Non-deterministic seed engine
    std::mt19937_64 gen(rd());  // replace rd() with a number to fix the seed
    std::uniform_real_distribution<> get_rand(-0.5, 0.5);

    std::vector<double> F( Npts ), G( Npts ), H( Npts ), I( Npts ),
        U1_lon( Npts ), U1_lat( Npts ),
        U2_lon( Npts ), U2_lat( Npts );

    std::generate(F.begin(), F.end(), [&]() { return get_rand(gen); });
    std::generate(G.begin(), G.end(), [&]() { return get_rand(gen); });

    std::generate(U1_lon.begin(), U1_lon.end(), [&]() { return get_rand(gen); });

    long double LHS, RHS, error;
    long double eps = 1e-14;

    //
    //// Test ddlon / ddlon_T
    //
    double dFdlon, dGdlonT;
    std::vector<double*> lon_deriv_vals  { &dFdlon },
                         lonT_deriv_vals { &dGdlonT };
    std::vector<const std::vector<double>*> deriv_fields = { &F },
                                            derivT_fields = { &G };

    // build the derivative fields
    for ( Ipt = 0; Ipt < Npts; Ipt++ ) {
        int Itime, Idepth, Ilat, Ilon;
        grid.index1to4_local( Ipt, Itime, Idepth, Ilat, Ilon );

        spher_derivative_at_point(
                lon_deriv_vals, deriv_fields, grid.longitude, "lon", grid,
                Itime, Idepth, Ilat, Ilon, mask);

        spher_derivative_transpose_at_point(
                lonT_deriv_vals, derivT_fields, grid.longitude, "lon", grid,
                Itime, Idepth, Ilat, Ilon, mask);

        H[Ipt] = dFdlon;
        I[Ipt] = dGdlonT;
    }
    LHS = 0, RHS = 0;
    for ( Ipt = 0; Ipt < Npts; Ipt++ ) {
        // LHS = < F, ddlon_T(G) >
        LHS += F[Ipt] * I[Ipt];

        // RHS = < ddlon(F), G >
        RHS += H[Ipt] * G[Ipt];
    }
    error = std::fabs( LHS - RHS ) / std::fmax( std::fabs( LHS ), std::fabs(RHS) );
    if ( error < eps ) {
        fprintf( stdout, "\033[32mddlon/ddlon_T adjoint test passed\033[0m\n" );
    } else {
        fprintf( stdout, "\033[31mddlon/ddlon_T adjoint test failed with error %.5Le (%.5Le vs %.5Le)\033[0m\n", error, LHS, RHS );
    }



    //
    //// Test ddlat / ddlat_T
    //

    double dFdlat, dGdlatT;
    std::vector<double*> lat_deriv_vals  { &dFdlat },
                         latT_deriv_vals { &dGdlatT };

    // build the derivative fields
    for ( Ipt = 0; Ipt < Npts; Ipt++ ) {
        int Itime, Idepth, Ilat, Ilon;
        grid.index1to4_local( Ipt, Itime, Idepth, Ilat, Ilon );

        spher_derivative_at_point(
                lat_deriv_vals, deriv_fields, grid.latitude, "lat", grid,
                Itime, Idepth, Ilat, Ilon, mask);

        spher_derivative_transpose_at_point(
                latT_deriv_vals, derivT_fields, grid.latitude, "lat", grid,
                Itime, Idepth, Ilat, Ilon, mask);

        H[Ipt] = dFdlat;
        I[Ipt] = dGdlatT;
    }
    LHS = 0, RHS = 0;
    for ( Ipt = 0; Ipt < Npts; Ipt++ ) {
        // LHS = < F, ddlat_T(G) >
        LHS += F[Ipt] * I[Ipt];

        // RHS = < ddlat(F), G >
        RHS += H[Ipt] * G[Ipt];
    }
    error = std::fabs( LHS - RHS ) / std::fmax( std::fabs( LHS ), std::fabs(RHS) );
    if ( error < eps ) {
        fprintf( stdout, "\033[32mddlat/ddlat_T adjoint test passed\033[0m\n" );
    } else {
        fprintf( stdout, "\033[31mddlat/ddlat_T adjoint test failed with error %.5Le (%.5Le vs %.5Le)\033[0m\n", error, LHS, RHS );
    }






    //
    //// Test grad / grad_T
    //

    grad_scalar( U2_lon, U2_lat, F, grid, mask );
    grad_scalar_T( H, U1_lon, U1_lat, grid, mask, false );

    LHS = 0, RHS = 0;
    for ( Ipt = 0; Ipt < Npts; Ipt++ ) {
        // LHS = < u, grad(F) >
        LHS += U1_lon[Ipt] * U2_lon[Ipt] + U1_lat[Ipt] * U2_lat[Ipt];

        // RHS = < grad_T(u), F >
        RHS += H[Ipt] * F[Ipt];
    }
    error = std::fabs( LHS - RHS ) / std::fmax( std::fabs( LHS ), std::fabs(RHS) );
    if ( error < eps ) {
        fprintf( stdout, "\033[32mgrad/grad_T adjoint test passed\033[0m\n" );
    } else {
        fprintf( stdout, "\033[31mgrad/grad_T adjoint test failed with error %.5Le (%.5Le vs %.5Le)\033[0m\n", error, LHS, RHS );
    }

    //
    //// Test perp_grad / perp_grad_T
    //

    perp_grad_scalar( U2_lon, U2_lat, F, grid, mask );
    perp_grad_scalar_T( H, U1_lon, U1_lat, grid, mask, false );

    LHS = 0, RHS = 0;
    for ( Ipt = 0; Ipt < Npts; Ipt++ ) {
        // LHS = < u, perp_grad(F) >
        LHS += U1_lon[Ipt] * U2_lon[Ipt] + U1_lat[Ipt] * U2_lat[Ipt];

        // RHS = < perp_grad_T(u), F >
        RHS += H[Ipt] * F[Ipt];
    }
    error = std::fabs( LHS - RHS ) / std::fmax( std::fabs( LHS ), std::fabs(RHS) );
    if ( error < eps ) {
        fprintf( stdout, "\033[32mperp_grad/perp_grad_T adjoint test passed\033[0m\n" );
    } else {
        fprintf( stdout, "\033[31mperp_grad/perp_grad_T adjoint test failed with error %.5Le (%.5Le vs %.5Le)\033[0m\n", error, LHS, RHS );
    }


    //
    //// Test grad / div
    //

    grad_scalar( U2_lon, U2_lat, F, grid, mask );
    div_vec( H, U1_lon, U1_lat, grid, mask );

    LHS = 0, RHS = 0;
    for ( Ipt = 0; Ipt < Npts; Ipt++ ) {
        long double dA = grid.areas[Ipt];

        // LHS = < u, grad(F) >
        LHS += (U1_lon[Ipt] * U2_lon[Ipt] + U1_lat[Ipt] * U2_lat[Ipt]) * dA;

        // RHS = < -div(u), F >
        RHS += (-H[Ipt] * F[Ipt]) * dA;
    }
    error = std::fabs( LHS - RHS ) / std::fmax( std::fabs( LHS ), std::fabs(RHS) );
    if ( error < eps ) {
        fprintf( stdout, "\033[32mgrad/div adjoint test passed\033[0m\n" );
    } else {
        fprintf( stdout, "\033[31mgrad/div adjoint test failed with error %.5Le (%.5Le vs %.5Le)\033[0m\n", error, LHS, RHS );
    }


    //
    //// Test perp_grad / vort
    //

    perp_grad_scalar( U2_lon, U2_lat, F, grid, mask );
    curl_vec_dot_er( H, U1_lon, U1_lat, grid, mask );

    LHS = 0, RHS = 0;
    for ( Ipt = 0; Ipt < Npts; Ipt++ ) {
        long double dA = grid.areas[Ipt];

        // LHS = < u, perp_grad(F) >
        LHS += (U1_lon[Ipt] * U2_lon[Ipt] + U1_lat[Ipt] * U2_lat[Ipt]) * dA;

        // RHS = < -vort(u), F >
        RHS += (-H[Ipt] * F[Ipt]) * dA;
    }
    error = std::fabs( LHS - RHS ) / std::fmax( std::fabs( LHS ), std::fabs(RHS) );
    if ( error < eps ) {
        fprintf( stdout, "\033[32mperp_grad/vort adjoint test passed\033[0m\n" );
    } else {
        fprintf( stdout, "\033[31mperp_grad/vort adjoint test failed with error %.5Le (%.5Le vs %.5Le)\033[0m\n", error, LHS, RHS );
    }


    //
    //// Test Laplacian
    //

    scalar_laplacian( H, F, grid, mask );
    scalar_laplacian( I, G, grid, mask );

    LHS = 0, RHS = 0;
    for ( Ipt = 0; Ipt < Npts; Ipt++ ) {
        long double dA = grid.areas[Ipt];

        // LHS = < Lap(F), G >
        LHS += (H[Ipt] * G[Ipt]) * dA;

        // RHS = < F, Lap(G) >
        RHS += (F[Ipt] * I[Ipt]) * dA;
    }
    error = std::fabs( LHS - RHS ) / std::fmax( std::fabs( LHS ), std::fabs(RHS) );
    if ( error < eps ) {
        fprintf( stdout, "\033[32mLaplacian adjoint test passed\033[0m\n" );
    } else {
        fprintf( stdout, "\033[31mLaplacian adjoint test failed with error %.5Le (%.5Le vs %.5Le)\033[0m\n", error, LHS, RHS );
    }


    //
    //// Test nullspace
    //

    std::fill( F.begin(), F.end(), 0. );
    scalar_laplacian( H, F, grid, mask );

    LHS = 0;
    RHS = 0;
    for ( Ipt = 0; Ipt < Npts; Ipt++ ) {
        long double dA = grid.areas[Ipt];

        // LHS = < Lap(const), Lap(const) > = || Lap(const) ||
        LHS += dA * H[Ipt] * H[Ipt];

        // RHS = 0
    }
    error = LHS;
    if ( error < eps ) {
        fprintf( stdout, "\033[32mLaplacian null-space test passed\033[0m\n" );
    } else {
        fprintf( stdout, "\033[31mLaplacian null-space test failed with error %.5Le (%.5Le vs %.5Le)\033[0m\n", error, LHS, RHS );
    }


}
