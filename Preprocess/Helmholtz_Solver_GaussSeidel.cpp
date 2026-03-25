#include "../constants.hpp"
#include "../functions.hpp"
#include "../netcdf_io.hpp"
#include "../preprocess.hpp"
#include "../differentiation_tools.hpp"
#include <algorithm>
#include <utility> // std::swap
#include <vector>
#include <deque>
#include <omp.h>
#include <math.h>
#include <numeric> // std::accumulate
//#include <Eigen/Sparse>
//#include <Eigen/IterativeLinearSolvers>
//#include <Eigen/QR>

// Pragma-magic to allow reduction over vector operator
//      thanks to: https://stackoverflow.com/questions/43168661/openmp-and-reduction-on-stdvector
#pragma omp declare reduction(vec_long_plus : std::vector<long double> : \
                              std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<long double>())) \
                    initializer(omp_priv = decltype(omp_orig)(omp_orig.size()))

void vector_linear_combine( 
        const double a,
        const std::vector<double> & x,
        const double b,
        const std::vector<double> & y,
        std::vector<double> & output
        ){

    // computes a*x + b*y

    const size_t Npts = x.size();
    #ifndef NDEBUG
    if ( Npts != y.size() ) {
        throw std::runtime_error("Size mismatch in linear_combine (x vs y).");
    }
    if ( Npts != output.size() ) {
        throw std::runtime_error("Size mismatch in linear_combine (x vs output).");
    }
    #endif

    long double accumulant = 0.;
    size_t Ipt;
    #pragma omp parallel default(none) private( Ipt ) \
    firstprivate( Npts, a, b ) \
    shared( x, y, output )
    {
        #pragma omp for collapse(1) schedule(static)
        for ( Ipt = 0; Ipt < Npts; Ipt++ ) {
            output[Ipt] = a * x[Ipt] + b * y[Ipt];
        }
    }

}

double vector_dot( 
        const std::vector<double> & a,
        const std::vector<double> & b,
        dataset * coarsened_grid
        ){

    const size_t Npts = a.size();
    #ifndef NDEBUG
    if ( Npts != b.size() ) {
        throw std::runtime_error("Size mismatch in dot.");
    }
    #endif

    long double accumulant = 0., denominator = 0.;
    size_t Ipt;
    #pragma omp parallel default(none) private( Ipt ) \
    firstprivate( Npts ) \
    shared( a, b, coarsened_grid ) \
    reduction( +:accumulant,denominator )
    {
        #pragma omp for collapse(1) schedule(static)
        for ( Ipt = 0; Ipt < Npts; Ipt++ ) {
            long double dA = coarsened_grid->areas[Ipt];
            accumulant += a[Ipt] * b[Ipt] * dA;
            denominator += dA;
        }
    }

    return (double) (accumulant / denominator);

}

void get_true_vort_div(
        dataset * coarsened_grid,
        dataset & PsiPhi_grid,
        std::vector<double> & cos_ip,
        std::vector<double> & cos_im
        ) {

    /*
    // Vorticity
    curl_vec_dot_er( 
            coarsened_grid->variables.at("vort"),
            coarsened_grid->variables.at("u_lon"),
            coarsened_grid->variables.at("u_lat"),
            PsiPhi_grid,
            coarsened_grid->mask
            );

    // Divergence
    div_vec( 
            coarsened_grid->variables.at("div"),
            coarsened_grid->variables.at("u_lon"),
            coarsened_grid->variables.at("u_lat"),
            PsiPhi_grid,
            coarsened_grid->mask
            );
    */

    //
    // get vorticity and divergence
    //

    const std::vector<double>   &latitude   = coarsened_grid->latitude,
                                &longitude  = coarsened_grid->longitude;
    const std::vector<int>  &myCounts = coarsened_grid->myCounts;

    const int   Ntime  = myCounts.at(0),
                Ndepth = myCounts.at(1),
                Nlat   = myCounts.at(2),
                Nlon   = myCounts.at(3);

    std::vector<double> &uo = coarsened_grid->variables.at("u_lon"),
                        &vo = coarsened_grid->variables.at("u_lat"),
                        &vort = coarsened_grid->variables.at("vort"),
                        &div  = coarsened_grid->variables.at("div" );

    const size_t Npts = div.size();

    // for now, assume uniform grid
    const double dlon_inv = 1. / ( coarsened_grid->longitude[1] - coarsened_grid->longitude[0]);
    const double dlat_inv = 1. / ( coarsened_grid->latitude[1]  - coarsened_grid->latitude[0] );
    const double R_E_inv = 1. / constants::R_earth;

    size_t Ipt;

    /*
    // this time using south and west neighbours
    #pragma omp parallel default(none) private( Ipt ) \
    firstprivate( Ntime, Ndepth, Nlat, Nlon, Npts, dlat_inv, dlon_inv ) \
    shared( coarsened_grid, PsiPhi_grid, vort, div, uo, vo )
    {

        #pragma omp for collapse(1) schedule(static)
        for ( Ipt = 0; Ipt < Npts; Ipt++ ) {
            int Itime, Idepth, Ilat, Ilon;
            coarsened_grid->index1to4_local(Ipt, Itime, Idepth, Ilat, Ilon);

            double cos_lat = 
                coarsened_grid->sin_cos_lat_pre_comp ?
                coarsened_grid->cos_lat[Ilat] :
                cos( coarsened_grid->latitude[Ilat] );
            const double R_E_cos_lat_inv = 1. / ( constants::R_earth * cos_lat );

            // point to right index
            int Ilon_W = ( Ilon > 0 ) ? Ilon - 1 : Nlon - 1;
            size_t E_ind = coarsened_grid->local_index( Itime, Idepth, Ilat, Ilon   );
            size_t W_ind = coarsened_grid->local_index( Itime, Idepth, Ilat, Ilon_W );

            // get point to north index (unless at north-most)
            int Ilat_S;
            if ( Ilat > 0 ) {
                Ilat_S = Ilat - 1;
            } else {
                Ilat_S = 0;
                Ilat = Ilat_S + 1;
            }
            size_t N_ind = coarsened_grid->local_index( Itime, Idepth, Ilat,   Ilon );
            size_t S_ind = coarsened_grid->local_index( Itime, Idepth, Ilat_S, Ilon );

            double cos_lat_N = 
                coarsened_grid->sin_cos_lat_pre_comp ?
                coarsened_grid->cos_lat[Ilat] :
                cos( coarsened_grid->latitude[Ilat] );

            double cos_lat_S = 
                coarsened_grid->sin_cos_lat_pre_comp ?
                coarsened_grid->cos_lat[Ilat_S] :
                cos( coarsened_grid->latitude[Ilat_S] );

            #ifdef NDEBUG
            vort[Ipt] = (
                  ( vo[E_ind]             - vo[W_ind]             ) * dlon_inv
                - ( uo[N_ind] * cos_lat_N - uo[S_ind] * cos_lat_S ) * dlat_inv
                ) * R_E_cos_lat_inv;
            div[ Ipt] = (
                  ( uo[E_ind]             - uo[W_ind]             ) * dlon_inv
                + ( vo[N_ind] * cos_lat_N - vo[S_ind] * cos_lat_S ) * dlat_inv
                ) * R_E_cos_lat_inv;
            #else
            vort.at(Ipt) = (
                  ( vo.at(E_ind)             - vo.at(W_ind)             ) * dlon_inv
                - ( uo.at(N_ind) * cos_lat_N - uo.at(S_ind) * cos_lat_S ) * dlat_inv
                ) * R_E_cos_lat_inv;
            div.at( Ipt) = (
                  ( uo.at(E_ind)             - uo.at(W_ind)             ) * dlon_inv
                + ( vo.at(N_ind) * cos_lat_N - vo.at(S_ind) * cos_lat_S ) * dlat_inv
                ) * R_E_cos_lat_inv;
            #endif
        }
    }
    */

    // get mean polar values of vo and uo
    double S_pole_vo = 0, N_pole_vo = 0., S_pole_uo = 0, N_pole_uo = 0;
    const double per_Nlon = 1. / Nlon;
    for ( int Ilon = 0; Ilon < Nlon; Ilon++ ) {
        size_t Ipt_N = coarsened_grid->local_index( 0, 0, Nlat-1, Ilon ); 
        size_t Ipt_S = coarsened_grid->local_index( 0, 0, 0,      Ilon ); 

        S_pole_uo += uo[Ipt_S] * per_Nlon;
        N_pole_uo += uo[Ipt_N] * per_Nlon;

        S_pole_vo += vo[Ipt_S] * per_Nlon;
        N_pole_vo += vo[Ipt_N] * per_Nlon;
    }

    //
    // Second loop sets gets vorticity and divergence
    //
    // this time using south and west neighbours

    #pragma omp parallel default(none) private( Ipt ) \
    firstprivate( Ntime, Ndepth, Nlat, Nlon, Npts, dlat_inv, dlon_inv,\
           S_pole_uo, S_pole_vo, N_pole_uo, N_pole_vo ) \
    shared( coarsened_grid, PsiPhi_grid, vort, div, uo, vo, cos_im, cos_ip )
    {

        #pragma omp for collapse(1) schedule(static)
        for ( Ipt = 0; Ipt < Npts; Ipt++ ) {
            int Itime, Idepth, Ilat, Ilon;
            coarsened_grid->index1to4_local(Ipt, Itime, Idepth, Ilat, Ilon);

            double cos_lat = 
                coarsened_grid->sin_cos_lat_pre_comp ?
                coarsened_grid->cos_lat[Ilat] :
                cos( coarsened_grid->latitude[Ilat] );
            const double R_E_cos_lat_inv = 1. / ( constants::R_earth * cos_lat );
            double  cos_plus  = cos_ip[Ilat],
                    cos_minus = cos_im[Ilat];

            // point to right index
            int Ilon_E = ( Ilon < Nlon - 1 ) ? Ilon + 1 : 0;
            int Ilon_W = ( Ilon > 0        ) ? Ilon - 1 : Nlon - 1;
            size_t E_ind = coarsened_grid->local_index( Itime, Idepth, Ilat, Ilon_E );
            size_t W_ind = coarsened_grid->local_index( Itime, Idepth, Ilat, Ilon_W );

            // get point to north index (unless at north-most)
            int Ilat_N = ( Ilat < Nlat - 1 ) ? Ilat + 1 : Ilat;
            int Ilat_S = ( Ilat > 0        ) ? Ilat - 1 : Ilat;
            size_t N_ind = coarsened_grid->local_index( Itime, Idepth, Ilat_N, Ilon );
            size_t S_ind = coarsened_grid->local_index( Itime, Idepth, Ilat_S, Ilon );

            #ifdef NDEBUG
            double uo_N = ( Ilat < Nlat - 1 ) ? uo[N_ind] : N_pole_uo;
            double uo_S = ( Ilat > 0        ) ? uo[S_ind] : S_pole_uo;
            double uo_E = uo[E_ind];
            double uo_W = uo[W_ind];

            double vo_N = ( Ilat < Nlat - 1 ) ? vo[N_ind] : N_pole_vo;
            double vo_S = ( Ilat > 0        ) ? vo[S_ind] : S_pole_vo;
            double vo_E = vo[E_ind];
            double vo_W = vo[W_ind];

            double uo_centre = uo[Ipt];
            double vo_centre = vo[Ipt];
            #else
            double uo_N = ( Ilat < Nlat - 1 ) ? uo.at(N_ind) : N_pole_uo;
            double uo_S = ( Ilat > 0        ) ? uo.at(S_ind) : S_pole_uo;
            double uo_E = uo.at(E_ind);
            double uo_W = uo.at(W_ind);

            double vo_N = ( Ilat < Nlat - 1 ) ? vo.at(N_ind) : N_pole_vo;
            double vo_S = ( Ilat > 0        ) ? vo.at(S_ind) : S_pole_vo;
            double vo_E = vo.at(E_ind);
            double vo_W = vo.at(W_ind);

            double uo_centre = uo.at(Ipt);
            double vo_centre = vo.at(Ipt);
            #endif

            // Divergence
            // zonal
            double tmp_div = ( uo_E - uo_W ) * dlon_inv * 0.5;

            // meridional
            tmp_div += (   cos_plus  * 0.5 * ( vo_N + vo_centre )
                         - cos_minus * 0.5 * ( vo_S + vo_centre )
                       ) * dlat_inv * 0.5;
            //tmp_div += (   cos_plus  * ( vo_N + vo_centre )
            //             - cos_minus * ( vo_S + vo_centre )
            //           ) * dlat_inv * R_E_cos_lat_inv * 0.5;
            tmp_div *= R_E_cos_lat_inv;

            // Vorticity
            // zonal
            double tmp_vort = ( vo_E - vo_W ) * dlon_inv * 0.5;

            // meridional
            tmp_vort -= (   cos_plus  * 0.5 * ( uo_N + uo_centre )
                          - cos_minus * 0.5 * ( uo_S + uo_centre )
                        ) * dlat_inv * 0.5;
            tmp_vort *= R_E_cos_lat_inv;
            //tmp_vort -= (   cos_plus  * ( uo_N + uo_centre )
            //              - cos_minus * ( uo_S + uo_centre )
            //            ) * dlat_inv * R_E_cos_lat_inv * 0.5;
            

            #ifdef NDEBUG
            vort[Ipt] = tmp_vort;
            div[ Ipt] = tmp_div;
            #else
            vort.at(Ipt) = tmp_vort;
            div.at( Ipt) = tmp_div;
            #endif
        }
    }
}


void get_vort_div_of_pot_tor(
        dataset * coarsened_grid,
        std::vector<double> & cos_ip,
        std::vector<double> & cos_im
        ) {

    const std::vector<double>   &latitude   = coarsened_grid->latitude,
                                &longitude  = coarsened_grid->longitude;
    const std::vector<int>  &myCounts = coarsened_grid->myCounts;

    const int   Ntime  = myCounts.at(0),
                Ndepth = myCounts.at(1),
                Nlat   = myCounts.at(2),
                Nlon   = myCounts.at(3);

    std::vector<double> &uo_tor = coarsened_grid->variables.at("proj_uo_tor"),
                        &uo_pot = coarsened_grid->variables.at("proj_uo_pot"),
                        &vo_tor = coarsened_grid->variables.at("proj_vo_tor"),
                        &vo_pot = coarsened_grid->variables.at("proj_vo_pot"),
                        &vort = coarsened_grid->variables.at("proj_vort_of_pot"),
                        &div  = coarsened_grid->variables.at("proj_div_of_tor" );

    const size_t Npts = div.size();

    // for now, assume uniform grid
    const double dlon_inv = 1. / ( coarsened_grid->longitude[1] - coarsened_grid->longitude[0]);
    const double dlat_inv = 1. / ( coarsened_grid->latitude[1]  - coarsened_grid->latitude[0] );
    const double R_E_inv = 1. / constants::R_earth;

    size_t Ipt;

    /*
    #pragma omp parallel default(none) private( Ipt ) \
    firstprivate( Ntime, Ndepth, Nlat, Nlon, Npts, dlat_inv, dlon_inv ) \
    shared( coarsened_grid, vort, div, uo_tor, uo_pot, vo_tor, vo_pot )
    {

        #pragma omp for collapse(1) schedule(static)
        for ( Ipt = 0; Ipt < Npts; Ipt++ ) {
            int Itime, Idepth, Ilat, Ilon;
            coarsened_grid->index1to4_local(Ipt, Itime, Idepth, Ilat, Ilon);

            double cos_lat = 
                coarsened_grid->sin_cos_lat_pre_comp ?
                coarsened_grid->cos_lat[Ilat] :
                cos( coarsened_grid->latitude[Ilat] );
            const double R_E_cos_lat_inv = 1. / ( constants::R_earth * cos_lat );

            // point to right index
            int Ilon_W = ( Ilon > 0 ) ? Ilon - 1 : Nlon - 1;
            size_t E_ind = coarsened_grid->local_index( Itime, Idepth, Ilat, Ilon   );
            size_t W_ind = coarsened_grid->local_index( Itime, Idepth, Ilat, Ilon_W );

            // get point to north index (unless at north-most)
            int Ilat_S;
            if ( Ilat > 0 ) {
                Ilat_S = Ilat - 1;
            } else {
                Ilat_S = 0;
                Ilat = Ilat_S + 1;
            }
            size_t N_ind = coarsened_grid->local_index( Itime, Idepth, Ilat,   Ilon );
            size_t S_ind = coarsened_grid->local_index( Itime, Idepth, Ilat_S, Ilon );

            double cos_lat_N = 
                coarsened_grid->sin_cos_lat_pre_comp ?
                coarsened_grid->cos_lat[Ilat] :
                cos( coarsened_grid->latitude[Ilat] );

            double cos_lat_S = 
                coarsened_grid->sin_cos_lat_pre_comp ?
                coarsened_grid->cos_lat[Ilat_S] :
                cos( coarsened_grid->latitude[Ilat_S] );

            #ifdef NDEBUG
            vort[Ipt] = (
                  ( vo_pot[E_ind]             - vo_pot[W_ind]             ) * dlon_inv
                - ( uo_pot[N_ind] * cos_lat_N - uo_pot[S_ind] * cos_lat_S ) * dlat_inv
                ) * R_E_cos_lat_inv;
            div[ Ipt] = (
                  ( uo_tor[E_ind]             - uo_tor[W_ind]             ) * dlon_inv
                + ( vo_tor[N_ind] * cos_lat_N - vo_tor[S_ind] * cos_lat_S ) * dlat_inv
                ) * R_E_cos_lat_inv;
            #else
            vort.at(Ipt) = (
                  ( vo_pot.at(E_ind)             - vo_pot.at(W_ind)             ) * dlon_inv
                - ( uo_pot.at(N_ind) * cos_lat_N - uo_pot.at(S_ind) * cos_lat_S ) * dlat_inv
                ) * R_E_cos_lat_inv;
            div.at( Ipt) = (
                  ( uo_tor.at(E_ind)             - uo_tor.at(W_ind)             ) * dlon_inv
                + ( vo_tor.at(N_ind) * cos_lat_N - vo_tor.at(S_ind) * cos_lat_S ) * dlat_inv
                ) * R_E_cos_lat_inv;
            #endif
        }
    }
    */


    #pragma omp parallel default(none) private( Ipt ) \
    firstprivate( Ntime, Ndepth, Nlat, Nlon, Npts, dlat_inv, dlon_inv ) \
    shared( coarsened_grid, vort, div, uo_tor, uo_pot, vo_tor, vo_pot, cos_ip, cos_im )
    {

        #pragma omp for collapse(1) schedule(static)
        for ( Ipt = 0; Ipt < Npts; Ipt++ ) {
            int Itime, Idepth, Ilat, Ilon;
            coarsened_grid->index1to4_local(Ipt, Itime, Idepth, Ilat, Ilon);

            double cos_lat = 
                coarsened_grid->sin_cos_lat_pre_comp ?
                coarsened_grid->cos_lat[Ilat] :
                cos( coarsened_grid->latitude[Ilat] );
            const double R_E_cos_lat_inv = 1. / ( constants::R_earth * cos_lat );
            double  cos_plus  = cos_ip[Ilat],
                    cos_minus = cos_im[Ilat];

            // point to right index
            int Ilon_E = ( Ilon < Nlon - 1 ) ? Ilon + 1 : 0;
            int Ilon_W = ( Ilon > 0        ) ? Ilon - 1 : Nlon - 1;
            size_t E_ind = coarsened_grid->local_index( Itime, Idepth, Ilat, Ilon_E );
            size_t W_ind = coarsened_grid->local_index( Itime, Idepth, Ilat, Ilon_W );

            // get point to north index (unless at north-most)
            int Ilat_N = ( Ilat < Nlat - 1 ) ? Ilat + 1 : Ilat;
            int Ilat_S = ( Ilat > 0        ) ? Ilat - 1 : Ilat;
            size_t N_ind = coarsened_grid->local_index( Itime, Idepth, Ilat_N, Ilon );
            size_t S_ind = coarsened_grid->local_index( Itime, Idepth, Ilat_S, Ilon );

            #ifdef NDEBUG
            double uo_N = uo_pot[N_ind];
            double uo_S = uo_pot[S_ind];
            double uo_E = uo_tor[E_ind];
            double uo_W = uo_tor[W_ind];

            double vo_N = vo_tor[N_ind];
            double vo_S = vo_tor[S_ind];
            double vo_E = vo_pot[E_ind];
            double vo_W = vo_pot[W_ind];

            double uo_centre = uo_pot[Ipt];
            double vo_centre = vo_tor[Ipt];
            #else
            double uo_N = uo_pot.at(N_ind);
            double uo_S = uo_pot.at(S_ind);
            double uo_E = uo_tor.at(E_ind);
            double uo_W = uo_tor.at(W_ind);

            double vo_N = vo_tor.at(N_ind);
            double vo_S = vo_tor.at(S_ind);
            double vo_E = vo_pot.at(E_ind);
            double vo_W = vo_pot.at(W_ind);

            double uo_centre = uo_pot.at(Ipt);
            double vo_centre = vo_tor.at(Ipt);
            #endif

            // Divergence
            // zonal
            long double tmp_div = ( uo_E - uo_W ) * dlon_inv * 0.5;

            // meridional
            tmp_div += (   cos_plus  * 0.5 * ( vo_N + vo_centre )
                         - cos_minus * 0.5 * ( vo_S + vo_centre )
                       ) * dlat_inv * 0.5;
            //tmp_div += (   cos_plus  * ( vo_N + vo_centre )
            //             - cos_minus * ( vo_S + vo_centre )
            //           ) * dlat_inv * R_E_cos_lat_inv * 0.5;
            tmp_div *= R_E_cos_lat_inv;

            // Vorticity
            // zonal
            long double tmp_vort = ( vo_E - vo_W ) * dlon_inv * 0.5;

            // meridional
            tmp_vort -= (   cos_plus  * 0.5 * ( uo_N + uo_centre )
                          - cos_minus * 0.5 * ( uo_S + uo_centre )
                        ) * dlat_inv * 0.5;
            //tmp_vort -= (   cos_plus  * ( uo_N + uo_centre )
            //              - cos_minus * ( uo_S + uo_centre )
            //            ) * dlat_inv * R_E_cos_lat_inv * 0.5;
            tmp_vort *= R_E_cos_lat_inv;
            

            #ifdef NDEBUG
            vort[Ipt] = tmp_vort;
            div[ Ipt] = tmp_div;
            #else
            vort.at(Ipt) = tmp_vort;
            div.at( Ipt) = tmp_div;
            #endif
        }
    }

}

void get_vels_vort_and_div(
        dataset * coarsened_grid,
        dataset & PsiPhi_grid,
        std::vector<double> & cos_ip,
        std::vector<double> & cos_im
            ) {

    /*
    // Toroidal velocity
    perp_grad_scalar( 
            coarsened_grid->variables.at("proj_uo_tor"),
            coarsened_grid->variables.at("proj_vo_tor"),
            PsiPhi_grid.variables.at("Psi"),
            *coarsened_grid,
            coarsened_grid->mask );

    // Vorticity
    curl_vec_dot_er( 
            coarsened_grid->variables.at("proj_vort"),
            coarsened_grid->variables.at("proj_uo_tor"),
            coarsened_grid->variables.at("proj_vo_tor"),
            *coarsened_grid,
            coarsened_grid->mask
            );

    // Potential velocity
    grad_scalar( 
            coarsened_grid->variables.at("proj_uo_pot"),
            coarsened_grid->variables.at("proj_vo_pot"),
            PsiPhi_grid.variables.at("Phi"),
            *coarsened_grid,
            coarsened_grid->mask );
    
    // Divergence
    div_vec( 
            coarsened_grid->variables.at("proj_div"),
            coarsened_grid->variables.at("proj_uo_pot"),
            coarsened_grid->variables.at("proj_vo_pot"),
            *coarsened_grid,
            coarsened_grid->mask
           );
    */

    /*
    // Similar to get_vels_vort_and_div, except this manually implements the derivatives
    // By lumping all of the derivative calculations together we can be more efficient, especially
    // with the threading.


    const std::vector<double>   &latitude   = coarsened_grid->latitude,
                                &longitude  = coarsened_grid->longitude;
    const std::vector<int>  &myCounts = coarsened_grid->myCounts;

    const int   Ntime  = myCounts.at(0),
                Ndepth = myCounts.at(1),
                Nlat   = myCounts.at(2),
                Nlon   = myCounts.at(3);

    size_t Ipt, In;

    // some short-hand vars
    std::vector<double> &proj_uo_tor = coarsened_grid->variables.at("proj_uo_tor"),
                        &proj_uo_pot = coarsened_grid->variables.at("proj_uo_pot"),
                        &proj_vo_tor = coarsened_grid->variables.at("proj_vo_tor"),
                        &proj_vo_pot = coarsened_grid->variables.at("proj_vo_pot"),
                        &proj_vort   = coarsened_grid->variables.at("proj_vort"),
                        &proj_div    = coarsened_grid->variables.at("proj_div" ),
                        &Psi         = PsiPhi_grid.variables["Psi"],
                        &Phi         = PsiPhi_grid.variables["Phi"];

    const size_t Npts = Psi.size();

    // zero out the fields
    #pragma omp parallel default(none) private( Ipt ) \
    shared( proj_uo_tor, proj_vo_tor, proj_uo_pot, proj_vo_pot, proj_vort, proj_div ) \
    firstprivate( Npts )
    {
        #pragma omp for collapse(1) schedule(static)
        for ( Ipt = 0; Ipt < Npts; Ipt++ ) {

            #ifdef NDEBUG
            proj_uo_tor[Ipt] = 0; 
            proj_vo_tor[Ipt] = 0;  
            proj_uo_pot[Ipt] = 0;
            proj_vo_pot[Ipt] = 0;

            proj_vort[Ipt] = 0;
            proj_div[ Ipt] = 0;
            #else
            proj_uo_tor.at(Ipt) = 0; 
            proj_vo_tor.at(Ipt) = 0;  
            proj_uo_pot.at(Ipt) = 0;
            proj_vo_pot.at(Ipt) = 0;

            proj_vort.at(Ipt) = 0;
            proj_div.at( Ipt) = 0;
            #endif
        }
    }

    const double R_E_inv = 1. / constants::R_earth;


    //
    // First loop set gets the velocity components
    //

    // for now, assume uniform grid
    const double dlon_inv = 1. / ( coarsened_grid->longitude[1] - coarsened_grid->longitude[0]);
    const double dlat_inv = 1. / ( coarsened_grid->latitude[1]  - coarsened_grid->latitude[0] );


    #pragma omp parallel default(none) private( Ipt, In ) \
    shared( coarsened_grid, PsiPhi_grid, Psi, Phi, \
         proj_uo_tor, proj_vo_tor, proj_uo_pot, proj_vo_pot ) \
    firstprivate( R_E_inv, Npts, Ntime, Ndepth, Nlat, Nlon, dlat_inv, dlon_inv )
    {

        #pragma omp for collapse(1) schedule(static)
        for ( Ipt = 0; Ipt < Npts; Ipt++ ) {

            int Itime, Idepth, Ilat, Ilon;
            coarsened_grid->index1to4_local(Ipt, Itime, Idepth, Ilat, Ilon);

            double cos_lat = 
                coarsened_grid->sin_cos_lat_pre_comp ?
                coarsened_grid->cos_lat[Ilat] :
                cos( coarsened_grid->latitude[Ilat] );
            const double R_E_cos_lat_inv = 1. / ( constants::R_earth * cos_lat );

            // point to right index
            int Ilon_E = ( Ilon < Nlon - 1 ) ? Ilon + 1 : 0;
            size_t E_ind = coarsened_grid->local_index( Itime, Idepth, Ilat, Ilon_E );
            size_t W_ind = coarsened_grid->local_index( Itime, Idepth, Ilat,   Ilon );

            // get point to north index (unless at north-most)
            int Ilat_N;
            if ( Ilat < Nlat - 1) {
                Ilat_N = Ilat + 1;
            } else {
                Ilat_N = Ilat;
                Ilat = Ilat_N - 1;
            }
            size_t N_ind = coarsened_grid->local_index( Itime, Idepth, Ilat_N, Ilon );
            size_t S_ind = coarsened_grid->local_index( Itime, Idepth, Ilat,   Ilon );


            #ifdef NDEBUG
            const double Psi_N = Psi[N_ind],
                         Psi_S = Psi[S_ind],
                         Psi_E = Psi[E_ind],
                         Psi_W = Psi[W_ind];
            const double Phi_N = Phi[N_ind],
                         Phi_S = Phi[S_ind],
                         Phi_E = Phi[E_ind],
                         Phi_W = Phi[W_ind];
            #else
            const double Psi_N = Psi.at(N_ind),
                         Psi_S = Psi.at(S_ind),
                         Psi_E = Psi.at(E_ind),
                         Psi_W = Psi.at(W_ind);
            const double Phi_N = Phi.at(N_ind),
                         Phi_S = Phi.at(S_ind),
                         Phi_E = Phi.at(E_ind),
                         Phi_W = Phi.at(W_ind);
            #endif

            //#ifdef NDEBUG
            proj_uo_tor[Ipt] = - ( Psi_N - Psi_S ) * R_E_inv         * dlat_inv;
            proj_vo_tor[Ipt] =   ( Psi_E - Psi_W ) * R_E_cos_lat_inv * dlon_inv;

            proj_uo_pot[Ipt] =   ( Phi_E - Phi_W ) * R_E_cos_lat_inv * dlon_inv;
            proj_vo_pot[Ipt] =   ( Phi_N - Phi_S ) * R_E_inv         * dlat_inv;
            //#endif
        }
    }


    //
    // Second loop sets gets vorticity and divergence
    //
    // this time using south and west neighbours

    #pragma omp parallel default(none) private( Ipt, In ) \
    firstprivate( Ntime, Ndepth, Nlat, Nlon, Npts, dlat_inv, dlon_inv ) \
    shared( coarsened_grid, PsiPhi_grid, proj_vort, proj_div, proj_uo_tor, \
           proj_vo_tor, proj_uo_pot, proj_vo_pot )
    {

        #pragma omp for collapse(1) schedule(static)
        for ( Ipt = 0; Ipt < Npts; Ipt++ ) {
            int Itime, Idepth, Ilat, Ilon;
            coarsened_grid->index1to4_local(Ipt, Itime, Idepth, Ilat, Ilon);

            double cos_lat = 
                coarsened_grid->sin_cos_lat_pre_comp ?
                coarsened_grid->cos_lat[Ilat] :
                cos( coarsened_grid->latitude[Ilat] );
            const double R_E_cos_lat_inv = 1. / ( constants::R_earth * cos_lat );

            // point to right index
            int Ilon_W = ( Ilon > 0 ) ? Ilon - 1 : Nlon - 1;
            size_t E_ind = coarsened_grid->local_index( Itime, Idepth, Ilat, Ilon   );
            size_t W_ind = coarsened_grid->local_index( Itime, Idepth, Ilat, Ilon_W );

            // get point to north index (unless at north-most)
            int Ilat_S;
            if ( Ilat > 0 ) {
                Ilat_S = Ilat - 1;
            } else {
                Ilat_S = 0;
                Ilat = Ilat_S + 1;
            }
            size_t N_ind = coarsened_grid->local_index( Itime, Idepth, Ilat,   Ilon );
            size_t S_ind = coarsened_grid->local_index( Itime, Idepth, Ilat_S, Ilon );

            double cos_lat_N = 
                coarsened_grid->sin_cos_lat_pre_comp ?
                coarsened_grid->cos_lat[Ilat] :
                cos( coarsened_grid->latitude[Ilat] );

            double cos_lat_S = 
                coarsened_grid->sin_cos_lat_pre_comp ?
                coarsened_grid->cos_lat[Ilat_S] :
                cos( coarsened_grid->latitude[Ilat_S] );

            #ifdef NDEBUG
            proj_vort[Ipt] = (
                  ( proj_vo_tor[E_ind]             - proj_vo_tor[W_ind]             ) * dlon_inv
                - ( proj_uo_tor[N_ind] * cos_lat_N - proj_uo_tor[S_ind] * cos_lat_S ) * dlat_inv
                ) * R_E_cos_lat_inv;
            proj_div[ Ipt] = (
                  ( proj_uo_pot[E_ind]             - proj_uo_pot[W_ind]             ) * dlon_inv
                + ( proj_vo_pot[N_ind] * cos_lat_N - proj_vo_pot[S_ind] * cos_lat_S ) * dlat_inv
                ) * R_E_cos_lat_inv;
            #else
            proj_vort.at(Ipt) = (
                  ( proj_vo_tor.at(E_ind)             - proj_vo_tor.at(W_ind)             ) * dlon_inv
                - ( proj_uo_tor.at(N_ind) * cos_lat_N - proj_uo_tor.at(S_ind) * cos_lat_S ) * dlat_inv
                ) * R_E_cos_lat_inv;
            proj_div.at( Ipt) = (
                  ( proj_uo_pot.at(E_ind)             - proj_uo_pot.at(W_ind)             ) * dlon_inv
                + ( proj_vo_pot.at(N_ind) * cos_lat_N - proj_vo_pot.at(S_ind) * cos_lat_S ) * dlat_inv
                ) * R_E_cos_lat_inv;
            #endif
        }
    }
    */

    // Similar to get_vels_vort_and_div, except this manually implements the derivatives
    // By lumping all of the derivative calculations together we can be more efficient, especially
    // with the threading.


    const std::vector<double>   &latitude   = coarsened_grid->latitude,
                                &longitude  = coarsened_grid->longitude;
    const std::vector<int>  &myCounts = coarsened_grid->myCounts;

    const int   Ntime  = myCounts.at(0),
                Ndepth = myCounts.at(1),
                Nlat   = myCounts.at(2),
                Nlon   = myCounts.at(3);

    size_t Ipt, In;

    // some short-hand vars
    std::vector<double> &proj_uo_tor = coarsened_grid->variables.at("proj_uo_tor"),
                        &proj_uo_pot = coarsened_grid->variables.at("proj_uo_pot"),
                        &proj_vo_tor = coarsened_grid->variables.at("proj_vo_tor"),
                        &proj_vo_pot = coarsened_grid->variables.at("proj_vo_pot"),
                        &proj_vort   = coarsened_grid->variables.at("proj_vort"),
                        &proj_div    = coarsened_grid->variables.at("proj_div" ),
                        &Psi         = PsiPhi_grid.variables["Psi"],
                        &Phi         = PsiPhi_grid.variables["Phi"];

    const size_t Npts = Psi.size();

    // get mean polar values of Psi and Phi
    double S_pole_Phi = 0, N_pole_Phi = 0., S_pole_Psi = 0, N_pole_Psi = 0;
    const double per_Nlon = 1. / Nlon;
    for ( int Ilon = 0; Ilon < Nlon; Ilon++ ) {
        size_t Ipt_N = coarsened_grid->local_index( 0, 0, Nlat-1, Ilon ); 
        size_t Ipt_S = coarsened_grid->local_index( 0, 0, 0, Ilon ); 

        S_pole_Phi += Phi[Ipt_S] * per_Nlon;
        N_pole_Phi += Phi[Ipt_N] * per_Nlon;

        S_pole_Psi += Psi[Ipt_S] * per_Nlon;
        N_pole_Psi += Psi[Ipt_N] * per_Nlon;
    }

    // zero out the fields
    #pragma omp parallel default(none) private( Ipt ) \
    shared( proj_uo_tor, proj_vo_tor, proj_uo_pot, proj_vo_pot, proj_vort, proj_div ) \
    firstprivate( Npts )
    {
        #pragma omp for collapse(1) schedule(static)
        for ( Ipt = 0; Ipt < Npts; Ipt++ ) {

            #ifdef NDEBUG
            proj_uo_tor[Ipt] = 0; 
            proj_vo_tor[Ipt] = 0;  
            proj_uo_pot[Ipt] = 0;
            proj_vo_pot[Ipt] = 0;

            proj_vort[Ipt] = 0;
            proj_div[ Ipt] = 0;
            #else
            proj_uo_tor.at(Ipt) = 0; 
            proj_vo_tor.at(Ipt) = 0;  
            proj_uo_pot.at(Ipt) = 0;
            proj_vo_pot.at(Ipt) = 0;

            proj_vort.at(Ipt) = 0;
            proj_div.at( Ipt) = 0;
            #endif
        }
    }

    const double R_E_inv = 1. / constants::R_earth;


    //
    // First loop set gets the velocity components
    //

    // for now, assume uniform grid
    const double dlon_inv = 1. / ( coarsened_grid->longitude[1] - coarsened_grid->longitude[0]);
    const double dlat_inv = 1. / ( coarsened_grid->latitude[1]  - coarsened_grid->latitude[0] );


    #pragma omp parallel default(none) private( Ipt, In ) \
    shared( coarsened_grid, PsiPhi_grid, Psi, Phi, cos_ip, cos_im, \
         proj_uo_tor, proj_vo_tor, proj_uo_pot, proj_vo_pot ) \
    firstprivate( R_E_inv, Npts, Ntime, Ndepth, Nlat, Nlon, dlat_inv, dlon_inv,\
           N_pole_Psi, N_pole_Phi, S_pole_Psi, S_pole_Phi )
    {

        #pragma omp for collapse(1) schedule(static)
        for ( Ipt = 0; Ipt < Npts; Ipt++ ) {

            int Itime, Idepth, Ilat, Ilon;
            coarsened_grid->index1to4_local(Ipt, Itime, Idepth, Ilat, Ilon);

            double cos_lat = 
                coarsened_grid->sin_cos_lat_pre_comp ?
                coarsened_grid->cos_lat[Ilat] :
                cos( coarsened_grid->latitude[Ilat] );
            const double R_E_cos_lat_inv = 1. / ( constants::R_earth * cos_lat );
            double  cos_plus  = cos_ip[Ilat],
                    cos_minus = cos_im[Ilat];

            // E-W indices
            int Ilon_E = ( Ilon < Nlon - 1 ) ? Ilon + 1 : 0;
            int Ilon_W = ( Ilon > 0        ) ? Ilon - 1 : Nlon - 1;
            size_t E_ind = coarsened_grid->local_index( Itime, Idepth, Ilat, Ilon_E );
            size_t W_ind = coarsened_grid->local_index( Itime, Idepth, Ilat, Ilon_W );

            // N-S indices
            int Ilat_N = ( Ilat < Nlat - 1 ) ? Ilat + 1 : Ilat;
            int Ilat_S = ( Ilat > 0        ) ? Ilat - 1 : Ilat;
            size_t N_ind = coarsened_grid->local_index( Itime, Idepth, Ilat_N, Ilon );
            size_t S_ind = coarsened_grid->local_index( Itime, Idepth, Ilat_S, Ilon );


            #ifdef NDEBUG
            const double Psi_N = ( Ilat < Nlat - 1 ) ? Psi[N_ind] : N_pole_Psi,
                         Psi_S = ( Ilat > 0        ) ? Psi[S_ind] : S_pole_Psi,
                         Psi_E = Psi[E_ind],
                         Psi_W = Psi[W_ind],
                         Psi_centre = Psi[Ipt];
            const double Phi_N = ( Ilat < Nlat - 1 ) ? Phi[N_ind] : N_pole_Phi,
                         Phi_S = ( Ilat > 0        ) ? Phi[S_ind] : S_pole_Phi,
                         Phi_E = Phi[E_ind],
                         Phi_W = Phi[W_ind],
                         Phi_centre = Phi[Ipt];
            #else
            const double Psi_N = ( Ilat < Nlat - 1 ) ? Psi.at(N_ind) : N_pole_Psi,
                         Psi_S = ( Ilat > 0        ) ? Psi.at(S_ind) : S_pole_Psi,
                         Psi_E = Psi.at(E_ind),
                         Psi_W = Psi.at(W_ind),
                         Psi_centre = Psi.at(Ipt);
            const double Phi_N = ( Ilat < Nlat - 1 ) ? Phi.at(N_ind) : N_pole_Phi,
                         Phi_S = ( Ilat > 0        ) ? Phi.at(S_ind) : S_pole_Phi,
                         Phi_E = Phi.at(E_ind),
                         Phi_W = Phi.at(W_ind),
                         Phi_centre = Phi.at(Ipt);
            #endif

            //#ifdef NDEBUG
            //proj_uo_tor[Ipt] = - ( cos_plus  * ( Psi_N - Psi_centre ) + 
            //                       cos_minus * ( Psi_centre - Psi_S ) 
            //                       ) * R_E_cos_lat_inv * dlat_inv * 0.5;
            proj_uo_tor[Ipt] = - ( cos_plus  * ( Psi_N + Psi_centre ) * 0.5 -
                                   cos_minus * ( Psi_centre + Psi_S ) * 0.5 
                                   ) * R_E_cos_lat_inv * dlat_inv * 0.5;
            proj_vo_tor[Ipt] =   ( Psi_E - Psi_W ) * R_E_cos_lat_inv * dlon_inv * 0.5;

            proj_uo_pot[Ipt] =   ( Phi_E - Phi_W ) * R_E_cos_lat_inv * dlon_inv * 0.5;
            //proj_vo_pot[Ipt] =   ( cos_plus  * ( Phi_N - Phi_centre ) + 
            //                       cos_minus * ( Phi_centre - Phi_S ) 
            //                       ) * R_E_cos_lat_inv * dlat_inv * 0.5;
            proj_vo_pot[Ipt] =   ( cos_plus  * ( Phi_N + Phi_centre ) * 0.5 -
                                   cos_minus * ( Phi_centre + Phi_S ) * 0.5 
                                   ) * R_E_cos_lat_inv * dlat_inv * 0.5;
            //#endif
        }
    }

    // get mean polar values of vo and uo
    double S_pole_vo = 0, N_pole_vo = 0., S_pole_uo = 0, N_pole_uo = 0;
    for ( int Ilon = 0; Ilon < Nlon; Ilon++ ) {
        size_t Ipt_N = coarsened_grid->local_index( 0, 0, Nlat-1, Ilon ); 
        size_t Ipt_S = coarsened_grid->local_index( 0, 0, 0, Ilon ); 

        S_pole_uo += (proj_uo_tor[Ipt_S] + proj_uo_pot[Ipt_S]) * per_Nlon;
        N_pole_uo += (proj_uo_tor[Ipt_N] + proj_uo_pot[Ipt_N]) * per_Nlon;

        S_pole_vo += (proj_vo_tor[Ipt_S] + proj_vo_pot[Ipt_S]) * per_Nlon;
        N_pole_vo += (proj_vo_tor[Ipt_N] + proj_vo_pot[Ipt_N]) * per_Nlon;
    }


    //
    // Second loop sets gets vorticity and divergence
    //
    // this time using south and west neighbours

    #pragma omp parallel default(none) private( Ipt, In ) \
    firstprivate( Ntime, Ndepth, Nlat, Nlon, Npts, dlat_inv, dlon_inv, \
           S_pole_uo, S_pole_vo, N_pole_uo, N_pole_vo ) \
    shared( coarsened_grid, PsiPhi_grid, proj_vort, proj_div, proj_uo_tor, \
           proj_vo_tor, proj_uo_pot, proj_vo_pot, cos_ip, cos_im )
    {

        #pragma omp for collapse(1) schedule(static)
        for ( Ipt = 0; Ipt < Npts; Ipt++ ) {
            int Itime, Idepth, Ilat, Ilon;
            coarsened_grid->index1to4_local(Ipt, Itime, Idepth, Ilat, Ilon);

            double cos_lat = 
                coarsened_grid->sin_cos_lat_pre_comp ?
                coarsened_grid->cos_lat[Ilat] :
                cos( coarsened_grid->latitude[Ilat] );
            const double R_E_cos_lat_inv = 1. / ( constants::R_earth * cos_lat );
            double  cos_plus  = cos_ip[Ilat],
                    cos_minus = cos_im[Ilat];

            // point to right index
            int Ilon_E = ( Ilon < Nlon - 1 ) ? Ilon + 1 : 0;
            int Ilon_W = ( Ilon > 0        ) ? Ilon - 1 : Nlon - 1;
            size_t E_ind = coarsened_grid->local_index( Itime, Idepth, Ilat, Ilon_E );
            size_t W_ind = coarsened_grid->local_index( Itime, Idepth, Ilat, Ilon_W );

            // get point to north index (unless at north-most)
            int Ilat_N = ( Ilat < Nlat - 1 ) ? Ilat + 1 : Ilat;
            int Ilat_S = ( Ilat > 0        ) ? Ilat - 1 : Ilat;
            size_t N_ind = coarsened_grid->local_index( Itime, Idepth, Ilat_N, Ilon );
            size_t S_ind = coarsened_grid->local_index( Itime, Idepth, Ilat_S, Ilon );

            #ifdef NDEBUG
            double uo_N = ( Ilat < Nlat - 1 ) ? proj_uo_pot[N_ind] + proj_uo_tor[N_ind] : N_pole_uo;
            double uo_S = ( Ilat < 0        ) ? proj_uo_pot[S_ind] + proj_uo_tor[S_ind] : S_pole_uo;
            double uo_E = proj_uo_pot[E_ind] + proj_uo_tor[E_ind];
            double uo_W = proj_uo_pot[W_ind] + proj_uo_tor[W_ind];

            double vo_N = ( Ilat < Nlat - 1 ) ? proj_vo_pot[N_ind] + proj_vo_tor[N_ind] : N_pole_vo;
            double vo_S = ( Ilat < 0        ) ? proj_vo_pot[S_ind] + proj_vo_tor[S_ind] : S_pole_vo;
            double vo_E = proj_vo_pot[E_ind] + proj_vo_tor[E_ind];
            double vo_W = proj_vo_pot[W_ind] + proj_vo_tor[W_ind];

            double uo_centre = proj_uo_pot[Ipt] + proj_uo_tor[Ipt];
            double vo_centre = proj_vo_pot[Ipt] + proj_vo_tor[Ipt];
            #else
            double uo_N = ( Ilat < Nlat - 1 ) ? proj_uo_pot.at(N_ind) + proj_uo_tor.at(N_ind) : N_pole_uo;
            double uo_S = ( Ilat < 0        ) ? proj_uo_pot.at(S_ind) + proj_uo_tor.at(S_ind) : S_pole_uo;
            double uo_E = proj_uo_pot.at(E_ind) + proj_uo_tor.at(E_ind);
            double uo_W = proj_uo_pot.at(W_ind) + proj_uo_tor.at(W_ind);

            double vo_N = ( Ilat < Nlat - 1 ) ? proj_vo_pot.at(N_ind) + proj_vo_tor.at(N_ind) : N_pole_vo;
            double vo_S = ( Ilat < 0        ) ? proj_vo_pot.at(S_ind) + proj_vo_tor.at(S_ind) : S_pole_vo;
            double vo_E = proj_vo_pot.at(E_ind) + proj_vo_tor.at(E_ind);
            double vo_W = proj_vo_pot.at(W_ind) + proj_vo_tor.at(W_ind);

            double uo_centre = proj_uo_pot.at(Ipt) + proj_uo_tor.at(Ipt);
            double vo_centre = proj_vo_pot.at(Ipt) + proj_vo_tor.at(Ipt);
            #endif

            // Divergence
            // zonal
            double tmp_div = ( uo_E - uo_W ) * dlon_inv * 0.5;

            // meridional
            tmp_div += (   cos_plus  * 0.5 * ( vo_N + vo_centre )
                         - cos_minus * 0.5 * ( vo_S + vo_centre )
                       ) * dlat_inv * 0.5;
            //tmp_div += (   cos_plus  * ( vo_N + vo_centre )
            //             - cos_minus * ( vo_S + vo_centre )
            //           ) * dlat_inv * R_E_cos_lat_inv * 0.5;
            tmp_div *= R_E_cos_lat_inv;

            // Vorticity
            // zonal
            double tmp_vort = ( vo_E - vo_W ) * dlon_inv * 0.5;

            // meridional
            tmp_vort -= (   cos_plus  * 0.5 * ( uo_N + uo_centre )
                          - cos_minus * 0.5 * ( uo_S + uo_centre )
                        ) * dlat_inv * 0.5;
            //tmp_vort -= (   cos_plus  * ( uo_N + uo_centre )
            //              - cos_minus * ( uo_S + uo_centre )
            //            ) * dlat_inv * R_E_cos_lat_inv * 0.5;
            tmp_vort *= R_E_cos_lat_inv;
            

            #ifdef NDEBUG
            proj_vort[Ipt] = tmp_vort;
            proj_div[ Ipt] = tmp_div;
            #else
            proj_vort.at(Ipt) = tmp_vort;
            proj_div.at( Ipt) = tmp_div;
            #endif
        }
    }


}

void get_norms_GS(
        long double & vel_2_err,
        long double & vort_2_err,
        long double & div_2_err,
        long double & vel_2_norm,
        long double & vort_2_norm,
        long double & div_2_norm,
        long double & vel_inf_err,
        long double & vort_inf_err,
        long double & div_inf_err,
        long double & vel_inf_norm,
        long double & vort_inf_norm,
        long double & div_inf_norm,
        long double & vort_2_ener,
        long double & div_2_ener,
        long double & Psi_2_norm,
        long double & Phi_2_norm,
        dataset * coarsened_grid
        ) {

    size_t Ipt;
    const size_t & Npts_coarse = coarsened_grid->mask.size();

    // We lose precision when summing a lot of numbers, so
    // step into higher precision for the accumulants
    long double total_area = 0;
    long double vel_2_err_L = 0, vort_2_err_L = 0, div_2_err_L = 0, 
                vel_2_norm_L = 0, vort_2_norm_L = 0, div_2_norm_L = 0, vort_2_ener_L = 0, 
                div_2_ener_L = 0, Phi_2_norm_L = 0, Psi_2_norm_L = 0;

    // Make some short-hand references
    std::vector<double> &proj_uo_tor = coarsened_grid->variables.at("proj_uo_tor"),
                        &proj_vo_tor = coarsened_grid->variables.at("proj_vo_tor"),
                        &proj_uo_pot = coarsened_grid->variables.at("proj_uo_pot"),
                        &proj_vo_pot = coarsened_grid->variables.at("proj_vo_pot");

    std::vector<double> &u_lon = coarsened_grid->variables.at("u_lon"),
                        &u_lat = coarsened_grid->variables.at("u_lat");

    std::vector<double> &vort = coarsened_grid->variables.at("vort"),
                        &div  = coarsened_grid->variables.at("div"),
                        &proj_vort = coarsened_grid->variables.at("proj_vort"),
                        &proj_div  = coarsened_grid->variables.at("proj_div");

    std::vector<double> &Phi = coarsened_grid->variables.at("Phi"),
                        &Psi = coarsened_grid->variables.at("Psi");

    const long double vort_mean = get_mean( vort, coarsened_grid ),
                      proj_vort_mean = get_mean( proj_vort, coarsened_grid );

    #pragma omp parallel default(none) \
    shared( coarsened_grid, u_lon, u_lat, vort, div, Phi, Psi, \
         proj_uo_tor, proj_vo_tor, proj_uo_pot, proj_vo_pot, proj_vort, proj_div ) \
    private( Ipt ) \
    firstprivate( Npts_coarse, vort_mean, proj_vort_mean ) \
    reduction( +:total_area,vel_2_err_L,vort_2_err_L,div_2_err_L,vel_2_norm_L,vort_2_norm_L,div_2_norm_L ) \
    reduction( +:vort_2_ener_L,div_2_ener_L,Phi_2_norm_L,Psi_2_norm_L )\
    reduction( max:vel_inf_err,vort_inf_err,div_inf_err,vel_inf_norm,vort_inf_norm,div_inf_norm )
    {
        #pragma omp for collapse(1) schedule(static)
        for ( Ipt = 0; Ipt < Npts_coarse; Ipt++ ) {

            long double area = coarsened_grid->areas[Ipt];
            long double incr;

            total_area += area;

            // Velocity Errors
            #ifdef NDEBUG
            incr = sqrt(
                       pow( proj_uo_tor[Ipt] + proj_uo_pot[Ipt] - u_lon[Ipt], 2)
                    +  pow( proj_vo_tor[Ipt] + proj_vo_pot[Ipt] - u_lat[Ipt], 2)
                    );
            #else
            incr = sqrt(
                    pow(
                           coarsened_grid->variables.at("proj_uo_tor").at(Ipt)
                         + coarsened_grid->variables.at("proj_uo_pot").at(Ipt)
                         - coarsened_grid->variables.at("u_lon").at(Ipt)
                        , 2)
                    +  pow(
                           coarsened_grid->variables.at("proj_vo_tor").at(Ipt)
                         + coarsened_grid->variables.at("proj_vo_pot").at(Ipt)
                         - coarsened_grid->variables.at("u_lat").at(Ipt)
                        , 2)
                    );
            #endif
            vel_inf_err = std::fmax( incr, vel_inf_err );

            vel_2_err_L += area * pow( incr, 2 ); 

            // Vorticity errors
            #ifdef NDEBUG
            incr = std::fabs( proj_vort[Ipt] - proj_vort_mean - (vort[Ipt] - vort_mean) );
            #else
            incr = std::fabs(
                  coarsened_grid->variables.at("proj_vort").at(Ipt)
                    - coarsened_grid->variables.at("vort").at(Ipt)
                    - ( proj_vort_mean - vort_mean )
                    );
            #endif
            vort_inf_err = std::fmax( incr, vort_inf_err );
            vort_2_err_L += area * pow( incr, 2 );

            #ifdef NDEBUG
            incr = std::fabs( proj_vort[Ipt] );
            #else
            incr = std::fabs( coarsened_grid->variables.at("proj_vort").at(Ipt) );
            #endif
            vort_2_ener_L += area * pow( incr, 2 );

            // Divergence errors
            #ifdef NDEBUG
            incr = std::fabs( proj_div[Ipt] - div[Ipt] );
            #else
            incr = std::fabs(
                      coarsened_grid->variables.at("proj_div").at(Ipt)
                    - coarsened_grid->variables.at("div").at(Ipt)
                    );
            #endif
            div_inf_err = std::fmax( incr, div_inf_err );
            div_2_err_L += area * pow( incr, 2 );

            #ifdef NDEBUG
            incr = std::fabs( proj_div[Ipt] );
            #else
            incr = std::fabs( coarsened_grid->variables.at("proj_div").at(Ipt) );
            #endif
            div_2_ener_L += area * pow( incr, 2 );

            // Reference norms
            #ifdef NDEBUG
            incr = sqrt( pow( u_lon[Ipt], 2) + pow( u_lat[Ipt], 2) );
            #else
            incr = sqrt(   pow(  coarsened_grid->variables.at("u_lon").at(Ipt), 2)
                         + pow(  coarsened_grid->variables.at("u_lat").at(Ipt), 2)
                    );
            #endif
            vel_inf_norm = std::fmax( incr, vel_inf_norm );
            vel_2_norm_L += area * pow( incr, 2 ); 

            #ifdef NDEBUG
            incr = std::fabs( vort[Ipt] );
            #else
            incr = std::fabs( coarsened_grid->variables.at("vort").at(Ipt) );
            #endif
            vort_inf_norm = std::fmax( incr, vort_inf_norm );
            vort_2_norm_L += area * pow( incr, 2 );

            #ifdef NDEBUG
            incr = std::fabs( div.at(Ipt) );
            #else
            incr = std::fabs( coarsened_grid->variables.at("div")[Ipt] );
            #endif
            div_inf_norm = std::fmax( incr, div_inf_norm );
            div_2_norm_L += area * pow( incr, 2 );

            // Psi and Phi 2 norms
            #ifdef NDEBUG
            incr = std::fabs( Psi[Ipt] );
            #else
            incr = std::fabs( coarsened_grid->variables.at("Psi").at(Ipt) );
            #endif
            Psi_2_norm_L += area * pow( incr, 2 );

            #ifdef NDEBUG
            incr = std::fabs( Phi.at(Ipt) );
            #else
            incr = std::fabs( coarsened_grid->variables.at("Phi")[Ipt] );
            #endif
            Phi_2_norm_L += area * pow( incr, 2 );
        }
    }

    const long double TA_inv = 1. / total_area;

    vel_2_err  = sqrt( vel_2_err_L  * TA_inv );
    vort_2_err = sqrt( vort_2_err_L * TA_inv );
    div_2_err  = sqrt( div_2_err_L  * TA_inv );

    vort_2_ener = sqrt( vort_2_ener_L * TA_inv );
    div_2_ener  = sqrt( div_2_ener_L  * TA_inv );

    vel_2_norm  = sqrt( vel_2_norm_L  * TA_inv );
    vort_2_norm = sqrt( vort_2_norm_L * TA_inv );
    div_2_norm  = sqrt( div_2_norm_L  * TA_inv );

    Psi_2_norm = sqrt( Psi_2_norm_L * TA_inv );
    Phi_2_norm = sqrt( Phi_2_norm_L * TA_inv );
}


void Helmholtz_Solver_GaussSeidel(
        const std::string output_fname,
        dataset & source_data,
        const double rel_tol,
        const unsigned int max_iters,
        const unsigned int iters_per_batch,
        const int num_refinements,
        const bool seed_provided,
        const std::vector<double> & Psi_seed,
        const std::vector<double> & Phi_seed,
        const MPI_Comm comm
        ) {

    int wRank, wSize;
    MPI_Comm_rank( comm, &wRank );
    MPI_Comm_size( comm, &wSize );

    // If we've passed the DO_TIMING flag, then create some timing vars
    Timing_Records timing_records;
    double clock_on;

    const size_t Npts = source_data.mask.size();

    const int Ntime = source_data.Ntime;
    const int Ndepth = source_data.Ndepth;

    // Fill in the land areas with zero velocity
    size_t index, Ipt;
    #pragma omp parallel default(none) \
    shared( source_data ) \
    private( index ) \
    firstprivate( Npts )
    {
        #pragma omp for collapse(1) schedule(static)
        for (index = 0; index < Npts; index++) {
            if (not(source_data.mask[index])) {
                source_data.variables["u_lon"][index] = 0.;
                source_data.variables["u_lat"][index] = 0.;
            }
        }
    }

    double mean_uo = get_mean( source_data.variables.at("u_lon"), &source_data );
    double mean_vo = get_mean( source_data.variables.at("u_lat"), &source_data );
    fprintf( stdout, "Mean values of zonal and meridional velocity are %.2e, %.2e\n", mean_uo, mean_vo );

    #if DEBUG >= 2
    fprintf(stdout, "Initializing class objects for data storage.\n");
    #endif
    HelmholtzDataClass Helmholtz_data;
    dataset solution_grid;
    dataset *coarsened_grid = new dataset();

    // Beginning the v-cycle [or, more acurrately, just a / cycle?]
    // Start at the coarsest resolution, solve, refine, repeat
    int Nlat_coarse = 0, Nlon_coarse = 0;
    size_t Npts_coarse = 0;
    for ( int refine_level = num_refinements; refine_level >= 0; refine_level-- ) {

        #if DEBUG >= 0
        fprintf(stdout, "\n\nBeginning refinement level %d.\n", refine_level);
        #endif
        std::vector<short int> saved_mask;

        // If refinement is turned on, do it
        if ( refine_level > 0 ) {
            Nlat_coarse = floor(sqrt( 0.5 * (Npts / pow(4., 2*refine_level )) ));
            Nlon_coarse = 2 * Nlat_coarse;
            Npts_coarse = (size_t) Nlat_coarse * Nlon_coarse;

            #if DEBUG >= 2
            fprintf(stdout, "Refinement level has Nlat x Nlon of %d x %d.\n", Nlat_coarse, Nlon_coarse);
            #endif


            // Downsample the velocities onto the coarse grid
            std::vector<std::string> vars_to_map = { "u_lon", "u_lat", "vort", "div" };

            if (constants::DO_TIMING) { clock_on = MPI_Wtime(); }
            initialize_coarsened_grid( *coarsened_grid, source_data, Nlat_coarse, Nlon_coarse );
            coarsened_grid->initialize_index_functions();
            coarsened_grid->compute_sin_cos_lat();
            if (constants::DO_TIMING) { timing_records.add_to_record(MPI_Wtime() - clock_on, "InitializeCoarseGrid"); }

            Npts_coarse = (size_t) Nlat_coarse * Nlon_coarse;
            for ( size_t Ivar = 0; Ivar < vars_to_map.size(); Ivar++ ) {
                coarsened_grid->variables.insert( std::pair< std::string, std::vector<double> >(
                            vars_to_map[Ivar], std::vector<double>(Npts_coarse, 0.) ) );
            }

            #if DEBUG >= 2
            fprintf(stdout, "Mapping velocities onto the coarse grid.\n");
            #endif
            vars_to_map = { "u_lon", "u_lat" };
            if (constants::DO_TIMING) { clock_on = MPI_Wtime(); }
            map_grid_to_grid( source_data, *coarsened_grid, vars_to_map );
            if (constants::DO_TIMING) { timing_records.add_to_record(MPI_Wtime() - clock_on, "Downsampling"); }

            coarsened_grid->adjacency_indices_adjoint.resize(0);


        } else {
            // Otherwise we're on the original grid, so just use the source
            coarsened_grid = &source_data;
            Npts_coarse = Npts;
            Nlat_coarse = coarsened_grid->Nlat;
            Nlon_coarse = coarsened_grid->Nlon;

            coarsened_grid->variables.insert( std::pair< std::string, std::vector<double> >(
                        "vort", std::vector<double>(Npts_coarse, 0.) ) );
            coarsened_grid->variables.insert( std::pair< std::string, std::vector<double> >(
                        "div", std::vector<double>(Npts_coarse, 0.) ) );

            coarsened_grid->initialize_index_functions();
            coarsened_grid->compute_sin_cos_lat();
            coarsened_grid->compute_cell_areas();

        }

        fprintf( stdout, "Npts_coarse = %'zu, len(areas) = %'zu\n", Npts_coarse, coarsened_grid->areas.size() );

        // We're done with the land mask, so remove it. This causes vort and
        // div to be computed on land cells to ensure agreement
        saved_mask = coarsened_grid->mask; // store a copy
        std::fill( coarsened_grid->mask.begin(), coarsened_grid->mask.end(), true );
        
        // Solver coefficients
        const bool use_cos_square_scaling = true;
        std::vector<double> Aii(Nlat_coarse, 0.),
                            coef_phi(Nlat_coarse, 0.),
                            coef_lambda(Nlat_coarse, 0.),
                            cos_ip(Nlat_coarse, 0.),
                            cos_ii(Nlat_coarse, 0.),
                            cos_im(Nlat_coarse, 0.);
        const double dphi    = std::fabs(coarsened_grid->latitude[ 1] - coarsened_grid->latitude[ 0]);
        const double dlambda = std::fabs(coarsened_grid->longitude[1] - coarsened_grid->longitude[0]);

        const double R_E_sqr = pow( constants::R_earth, 2. );

        for ( int Ilat = 0; Ilat < Nlat_coarse; Ilat++ ) {
            const double cos_i = cos(coarsened_grid->latitude[Ilat]);
            const double cos_ip_val = ( Ilat < Nlat_coarse - 1 ) ?
                                      cos( 0.5 * (coarsened_grid->latitude[Ilat+1]
                                                  + coarsened_grid->latitude[Ilat]) ) :
                                      cos( 0.5 * ((M_PI/2)
                                                  + coarsened_grid->latitude[Ilat]) ) ;
                                      //0.5 * ( cos_i + cos(coarsened_grid->latitude[Ilat+1]) ) :
                                      //0.5 *   cos_i;
            const double cos_im_val = ( Ilat > 0 ) ?
                                      cos( 0.5 * (coarsened_grid->latitude[Ilat-1]
                                                  + coarsened_grid->latitude[Ilat]) ) :
                                      cos( 0.5 * (-(M_PI/2)
                                                  + coarsened_grid->latitude[Ilat]) ) ;
                                      //0.5 * ( cos_i + cos(coarsened_grid->latitude[Ilat-1]) ) :
                                      //0.5 *   cos_i;

            const double inv_cos = 1. / cos_i;
            const double inv_cos_2 = inv_cos * inv_cos;

            double alpha_lambda, alpha_phi;
            if (use_cos_square_scaling) {
                alpha_lambda = 1. / (4. * dlambda * dlambda);
                alpha_phi    = 1. / (dphi * dphi);
            } else {
                alpha_lambda = inv_cos_2 / (4. * dlambda * dlambda);
                alpha_phi    = inv_cos_2 / (dphi * dphi);
            }

            //const double Aii_val = 1. / ( coef_theta_val * ( cos_ip_val + cos_im_val ) + 2.0 * coef_phi_val );
            const double Aii_val = 1. / ( 2 * (alpha_lambda + alpha_phi * cos_ip_val * cos_im_val) );

            Aii[Ilat] = Aii_val;
            coef_phi[Ilat] = alpha_phi;
            coef_lambda[Ilat] = alpha_lambda;
            cos_ii[Ilat] = cos_i;
            cos_ip[Ilat] = cos_ip_val;
            cos_im[Ilat] = cos_im_val;
        }

        // Get vorticity and divergence
        #if DEBUG >= 2
        fprintf(stdout, "Computing vorticity and divergence on coarse grid.\n");
        #endif
        if (constants::DO_TIMING) { clock_on = MPI_Wtime(); }

        get_true_vort_div( coarsened_grid, *coarsened_grid, cos_ip, cos_im );

        double mean_vort, mean_div;
        demean( coarsened_grid->variables.at("vort"), coarsened_grid, mean_vort );
        demean( coarsened_grid->variables.at("div" ), coarsened_grid, mean_div  );
        fprintf( stdout, "Mean (coarse-grid) values of vorticity and divergence are %.2e, %.2e\n", mean_vort, mean_div );

        double mean_uo_coarse = get_mean( coarsened_grid->variables.at("u_lon"), coarsened_grid );
        double mean_vo_coarse = get_mean( coarsened_grid->variables.at("u_lat"), coarsened_grid );
        fprintf( stdout, "Mean (coarse-grid) values of zonal and meridional velocity are %.2e, %.2e\n", mean_uo_coarse, mean_vo_coarse );

        if (constants::DO_TIMING) { timing_records.add_to_record(MPI_Wtime() - clock_on, "Downsampling"); }

        // Add Psi and Phi variables to our coarsened grid
        coarsened_grid->variables.insert( std::pair< std::string, std::vector<double> >(
                    "Psi", std::vector<double>(Npts_coarse, 0.) ) );
        coarsened_grid->variables.insert( std::pair< std::string, std::vector<double> >(
                    "Phi", std::vector<double>(Npts_coarse, 0.) ) );

        // Now, up-sample the previous solutions onto the new grid
        if ( refine_level < num_refinements ) {

            std::vector<std::string> vars_to_map = { "Psi", "Phi" };

            if (constants::DO_TIMING) { clock_on = MPI_Wtime(); }
            map_grid_to_grid( solution_grid, *coarsened_grid, vars_to_map );
            if (constants::DO_TIMING) { timing_records.add_to_record(MPI_Wtime() - clock_on, "Upsampling"); }
        }

        // Re-initialize the solution grid on the new computational grid
        // To avoid having to re-build, just copy from the coarsened grid
        solution_grid.clear();
        solution_grid.copy_from_ptr( coarsened_grid );
        solution_grid.initialize_index_functions();
        solution_grid.compute_sin_cos_lat();
        std::fill( solution_grid.mask.begin(), solution_grid.mask.end(), true ); 

        solution_grid.variables.insert( std::pair< std::string, std::vector<double> >(
                    "Psi", std::vector<double>(Npts_coarse, 0.) ) );
        solution_grid.variables.insert( std::pair< std::string, std::vector<double> >(
                    "Phi", std::vector<double>(Npts_coarse, 0.) ) );


        // Reset the Helmholtz data for this resolution
        Helmholtz_data.clear();
        Helmholtz_data.tolerance = rel_tol;

        //Helmholtz_data.iterations_per_cycle = Nlon_coarse *   10;
        //Helmholtz_data.iterations_per_cycle = Nlon_coarse *   1;
        //Helmholtz_data.iteration_max        = Nlon_coarse * 1000;
        Helmholtz_data.iterations_per_cycle = 100;
        Helmholtz_data.iteration_max        = 10000;
        //Helmholtz_data.iterations_per_cycle = refine_level == 0 ?   1000 : Npts_coarse / 10;
        //Helmholtz_data.iteration_max        = refine_level == 0 ? 100000 : Npts_coarse * 10;
        //Helmholtz_data.iterations_per_cycle =   1000;
        //Helmholtz_data.iteration_max        = 100000;

        // Use spherical harmonics to set the seed
        coarsened_grid->variables.insert( std::pair< std::string, std::vector<double> >(
                    "Psi", std::vector<double>(Npts_coarse, 0.) ) );
        coarsened_grid->variables.insert( std::pair< std::string, std::vector<double> >(
                    "Phi", std::vector<double>(Npts_coarse, 0.) ) );
        if ( (num_refinements == 0) and (seed_provided) ) {
            coarsened_grid->variables["Psi"] = Psi_seed;
            coarsened_grid->variables["Phi"] = Phi_seed;
        } else if ( refine_level == num_refinements ) {
            if (constants::DO_TIMING) { clock_on = MPI_Wtime(); }
            /*
            SphericalHarmonicSolver(  
                    coarsened_grid->variables.at("Psi"),
                    coarsened_grid->variables.at("Phi"),
                   *coarsened_grid,
                    coarsened_grid->variables.at("vort"),
                    coarsened_grid->variables.at("div") 
                    );
            */
            if (constants::DO_TIMING) { timing_records.add_to_record(MPI_Wtime() - clock_on, "SphericalHarmonicSolver"); }
        }


        //
        //// Output the seed / starting guess
        //
        std::vector<std::string> vars_to_add = { 
            "projected_u_lon", "projected_u_lat", "projected_vort", "projected_div",
            "proj_uo_tor", "proj_uo_pot", "proj_vo_tor", "proj_vo_pot" ,
            "proj_vort", "proj_div"
        };
        for ( size_t Ivar = 0; Ivar < vars_to_add.size(); Ivar++ ) {
            coarsened_grid->variables.insert( std::pair< std::string, std::vector<double> >(
                        vars_to_add[Ivar], std::vector<double>(Npts_coarse, 0.) ) );
        }

        // Compute velocity components
        get_vels_vort_and_div( coarsened_grid, *coarsened_grid, cos_ip, cos_im );
        demean( coarsened_grid->variables.at("proj_vort"), coarsened_grid, mean_vort );
        demean( coarsened_grid->variables.at("proj_div" ), coarsened_grid, mean_div  );

        for ( size_t Ipt = 0; Ipt < Npts_coarse; Ipt++ ) {
            coarsened_grid->variables.at("projected_u_lon")[Ipt] = 
                  coarsened_grid->variables.at("proj_uo_tor")[Ipt] 
                + coarsened_grid->variables.at("proj_uo_pot")[Ipt];
            coarsened_grid->variables.at("projected_u_lat")[Ipt] = 
                  coarsened_grid->variables.at("proj_vo_tor")[Ipt] 
                + coarsened_grid->variables.at("proj_vo_pot")[Ipt];
        }

        coarsened_grid->variables.at("projected_vort") = coarsened_grid->variables.at("proj_vort");
        coarsened_grid->variables.at("projected_div") = coarsened_grid->variables.at("proj_div");

        // As a sanity check, get the divergence of the toroidal field and vorticity of the potential field
        coarsened_grid->variables.insert( std::pair< std::string, std::vector<double> >(
                    "proj_vort_of_pot", std::vector<double>(Npts_coarse, 0.) ) );
        coarsened_grid->variables.insert( std::pair< std::string, std::vector<double> >(
                    "proj_div_of_tor", std::vector<double>(Npts_coarse, 0.) ) );
        get_vort_div_of_pot_tor( coarsened_grid, cos_ip, cos_im );


        // Write the seed to output
        std::ostringstream seed_stream;
        seed_stream << "seed_R" << (refine_level) << ".nc";
        const std::string seed_filename = seed_stream.str();
        solution_grid.variables.at("Psi") = coarsened_grid->variables.at("Psi");
        solution_grid.variables.at("Phi") = coarsened_grid->variables.at("Phi");
        write_Helmholtz_output( seed_filename, Helmholtz_data, *coarsened_grid, solution_grid );

        // Make some short-hand references
        std::vector<double> &soln_Psi = solution_grid.variables.at("Psi"),
                            &soln_Phi = solution_grid.variables.at("Phi"); 

        std::vector<double> &cgrid_proj_vort = coarsened_grid->variables.at("proj_vort"),
                            &cgrid_proj_div  = coarsened_grid->variables.at("proj_div"),
                            &cgrid_vort = coarsened_grid->variables.at("vort"),
                            &cgrid_div  = coarsened_grid->variables.at("div");


        // compute errors and check convergence
        long double vel_2_err  = 0, vel_inf_err  = 0,
                    vort_2_err = 0, vort_inf_err = 0,
                    div_2_err  = 0, div_inf_err  = 0;
        long double vel_2_norm  = 0, vel_inf_norm  = 0,
                    vort_2_norm = 0, vort_inf_norm = 0,
                    div_2_norm  = 0, div_inf_norm  = 0;
        long double vort_2_ener = 0, div_2_ener = 0;
        long double Psi_2_norm = 0, Phi_2_norm = 0;

        get_norms_GS(  vel_2_err,  vort_2_err,  div_2_err,
                vel_2_norm, vort_2_norm, div_2_norm,
                vel_inf_err,  vort_inf_err,  div_inf_err,
                vel_inf_norm, vort_inf_norm, div_inf_norm,
                vort_2_ener, div_2_ener,
                Psi_2_norm, Phi_2_norm,
                coarsened_grid
                );
        fprintf( stdout,"Solver starting with relative 2-norm errors of vels, vort, and div at  %.2Le, %.2Le, and %.2Le\n", 
                vel_2_err / vel_2_norm,
                vort_2_err / vort_2_norm,
                div_2_err / div_2_norm
               );

        double min_div_err = 1e10, min_vort_err = 1e10;
        int vort_stagnation_count = 0, div_stagnation_count = 0;
        std::vector<double> dPsidt(Npts_coarse, 0), dPhidt(Npts_coarse, 0);

        bool keep_solving = true;
        unsigned long int iter_cycle = 0;

        // Storage for polar values
        std::vector<double> N_poles_Psi_vals( Ntime*Ndepth, 0. ),
                            S_poles_Psi_vals( Ntime*Ndepth, 0. ),
                            N_poles_Phi_vals( Ntime*Ndepth, 0. ),
                            S_poles_Phi_vals( Ntime*Ndepth, 0. );

        // We're going to pay the memory cost of storing
        // the 4-indices for each 1-index to get the speed-up
        std::vector<int> Itime_array(Npts_coarse, 0),
                         Idepth_array(Npts_coarse, 0),
                         Ilat_array(Npts_coarse, 0),
                         Ilon_array(Npts_coarse, 0);
        #pragma omp parallel default(none) private(Ipt) \
        shared( coarsened_grid, Itime_array, Idepth_array, Ilat_array, Ilon_array ) \
        firstprivate( Npts_coarse )
        {
            #pragma omp for collapse(1) schedule(static)
            for ( Ipt = 0; Ipt < Npts_coarse; Ipt++ ) {
                int Itime, Idepth, Ilat, Ilon;
                coarsened_grid->index1to4_local(Ipt, Itime, Idepth, Ilat, Ilon);

                Itime_array[ Ipt] = Itime;
                Idepth_array[Ipt] = Idepth;
                Ilat_array[  Ipt] = Ilat;
                Ilon_array[  Ipt] = Ilon;
            }
        }

        // BiCGSTAB vectors
        std::vector<double> 
            nu_Psi(Npts_coarse, 0.), nu_Phi( Npts_coarse, 0. ),
            h_Psi(Npts_coarse, 0.), h_Phi( Npts_coarse, 0. ),
            s_Psi(Npts_coarse, 0.), s_Phi(Npts_coarse, 0.),
            t_Psi(Npts_coarse, 0.), t_Phi(Npts_coarse, 0.),
            x_Psi(Npts_coarse, 0.), x_Phi(Npts_coarse, 0.),
            p_i_Psi(Npts_coarse, 0.), p_i_Phi(Npts_coarse, 0.),
            p_im1_Psi(Npts_coarse, 0.), p_im1_Phi(Npts_coarse, 0.),
            r_i_Psi(Npts_coarse, 0.), r_i_Phi(Npts_coarse, 0.),
            r_im1_Psi(Npts_coarse, 0.), r_im1_Phi(Npts_coarse, 0.);

        //std::vector<double> &nu_Psi = coarsened_grid->variables.at("proj_vort"),
        //                    &nu_Phi = coarsened_grid->variables.at("proj_div");

        // BiCGSTAB scalars
        double alpha_Psi, alpha_Phi,
               omega_Psi, omega_Phi,
               rho_i_Psi, rho_i_Phi,
               rho_im1_Psi, rho_im1_Phi,
               beta_Psi, beta_Phi;

        // BiCGSTAB initial values

        // x0, inherited from soln_Psi and soln_Phi
        x_Psi = soln_Psi;
        x_Phi = soln_Phi;
        get_vels_vort_and_div( coarsened_grid, solution_grid, cos_ip, cos_im );

        /*
        fprintf( stdout, "dot(cgrid_vort, cgrid_vort) = %.5e\n", vector_dot( cgrid_vort, cgrid_vort ) );
        fprintf( stdout, "dot(cgrid_div,  cgrid_div ) = %.5e\n", vector_dot( cgrid_div,  cgrid_div  ) );

        fprintf( stdout, "dot(cgrid_proj_vort, cgrid_proj_vort) = %.5e\n", 
                vector_dot( cgrid_proj_vort, cgrid_proj_vort ) );
        fprintf( stdout, "dot(cgrid_proj_div,  cgrid_proj_div ) = %.5e\n", 
                vector_dot( cgrid_proj_div,  cgrid_proj_div  ) );
        */

        // r_0 = b - Ax0
        vector_linear_combine( 1, cgrid_vort, -1, cgrid_proj_vort, r_im1_Psi );
        vector_linear_combine( 1, cgrid_div,  -1, cgrid_proj_div,  r_im1_Phi );

        // choose r_hat s.t. dot(r_hat, r_0) != 0
        std::vector<double> r_hat_Psi = r_im1_Psi;
        std::vector<double> r_hat_Phi = r_im1_Phi;

        // rho_0 = dot( r_hat, r_0 )
        rho_im1_Psi = vector_dot( r_hat_Psi, r_im1_Psi, coarsened_grid );
        rho_im1_Phi = vector_dot( r_hat_Phi, r_im1_Phi, coarsened_grid );

        // p_0 = r_0
        p_im1_Psi = r_im1_Psi;
        p_im1_Phi = r_im1_Phi;

        while ( keep_solving ) {

            if (constants::DO_TIMING) { clock_on = MPI_Wtime(); }

            // nu = A(p_{i-1})
            //  i.e. nu_Psi and nu_Phi are the projected vorticity and divergence
            soln_Psi = p_im1_Psi;
            soln_Phi = p_im1_Phi;
            get_vels_vort_and_div( coarsened_grid, solution_grid, cos_ip, cos_im );
            nu_Psi = cgrid_proj_vort;
            nu_Phi = cgrid_proj_div;
            
            // alpha = rho_{i-1} / dot( r_hat, nu )
            alpha_Psi = rho_im1_Psi / vector_dot( r_hat_Psi, nu_Psi, coarsened_grid );
            alpha_Phi = rho_im1_Phi / vector_dot( r_hat_Phi, nu_Phi, coarsened_grid );

            // h = x_{i-1} + alpha * p_{i-1}
            vector_linear_combine( 1, x_Psi, alpha_Psi, p_im1_Psi, h_Psi );
            vector_linear_combine( 1, x_Phi, alpha_Phi, p_im1_Phi, h_Phi );

            // s = r_{i-1} - alpha * nu
            vector_linear_combine( 1, r_im1_Psi, -alpha_Psi, nu_Psi, s_Psi );
            vector_linear_combine( 1, r_im1_Phi, -alpha_Phi, nu_Phi, s_Phi );

            //fprintf( stdout, "dot(s_Psi, s_Psi) = %.5e\n", vector_dot( s_Psi, s_Psi ) );
            //fprintf( stdout, "dot(s_Phi, s_Phi) = %.5e\n", vector_dot( s_Phi, s_Phi ) );

            // t = A(s)
            //  i.e. t_Psi and t_Phi are the Laplacians of s_Psi and s_Phi (their vort and div)
            soln_Psi = s_Psi;
            soln_Phi = s_Phi;
            get_vels_vort_and_div( coarsened_grid, solution_grid, cos_ip, cos_im );
            t_Psi = cgrid_proj_vort;
            t_Phi = cgrid_proj_div;
    
            // omega = dot( t, s ) / dot( t, t )
            omega_Psi = vector_dot( t_Psi, s_Psi, coarsened_grid ) / vector_dot( t_Psi, t_Psi, coarsened_grid );
            omega_Phi = vector_dot( t_Phi, s_Phi, coarsened_grid ) / vector_dot( t_Phi, t_Phi, coarsened_grid );

            // x_i = h + omega * s
            vector_linear_combine( 1., h_Psi, omega_Psi, s_Psi, x_Psi );
            vector_linear_combine( 1., h_Phi, omega_Phi, s_Phi, x_Phi );

            // r_i = s - omega * t
            vector_linear_combine( 1., s_Psi, -omega_Psi, t_Psi, r_i_Psi );
            vector_linear_combine( 1., s_Phi, -omega_Phi, t_Phi, r_i_Phi );

            if ( iter_cycle % Helmholtz_data.iterations_per_cycle == 0 ) {
                fprintf( stdout, "dot(r_i_Psi, r_i_Psi) = %.5e\n", vector_dot( r_i_Psi, r_i_Psi, coarsened_grid ) );
                fprintf( stdout, "dot(r_i_Phi, r_i_Phi) = %.5e\n", vector_dot( r_i_Phi, r_i_Phi, coarsened_grid ) );
            }

            // if ||r_i|| < tolerance, halt. x_i is solution
            //  this is done at the end

            // rho_i = dot( r_hat, r_i )
            rho_i_Psi = vector_dot( r_hat_Psi, r_i_Psi, coarsened_grid );
            rho_i_Phi = vector_dot( r_hat_Phi, r_i_Phi, coarsened_grid );

            // beta = ( rho_i / rho_{i-1} ) * ( alpha / omega )
            beta_Psi = ( rho_i_Psi / rho_im1_Psi ) * ( alpha_Psi / omega_Psi );
            beta_Phi = ( rho_i_Phi / rho_im1_Phi ) * ( alpha_Phi / omega_Phi );

            // p_i = r_i + beta * ( p_{i-1} - omega * nu )
            //  we will store [p_{i-1} - omega * nu] as t to re-use functions
            vector_linear_combine( 1., p_im1_Psi, -omega_Psi, nu_Psi, t_Psi );
            vector_linear_combine( 1., p_im1_Phi, -omega_Phi, nu_Phi, t_Phi );

            vector_linear_combine( 1., r_i_Psi, beta_Psi, t_Psi, p_i_Psi );
            vector_linear_combine( 1., r_i_Phi, beta_Phi, t_Phi, p_i_Phi );

        
            // re-index for next round i.e. current becomes previous
            std::swap( p_i_Psi, p_im1_Psi );
            std::swap( p_i_Phi, p_im1_Phi );

            std::swap( r_i_Psi, r_im1_Psi );
            std::swap( r_i_Phi, r_im1_Phi );

            std::swap( rho_i_Psi, rho_im1_Psi );
            std::swap( rho_i_Phi, rho_im1_Phi );
            


            // De-drift Psi, Phi by setting Psi[0] = Phi[0] = 0
            const double Psi_shift = soln_Psi[0],
                         Phi_shift = soln_Phi[0];
            #pragma omp parallel default(none) private(Ipt) \
            shared( soln_Phi, soln_Psi ) \
            firstprivate( Npts_coarse, Psi_shift, Phi_shift )
            {
                #pragma omp for collapse(1) schedule(static)
                for ( Ipt = 0; Ipt < Npts_coarse; Ipt++ ) {
                    soln_Phi[Ipt] -= Phi_shift;
                    soln_Psi[Ipt] -= Psi_shift;
                }
            }

            /*
            // Update the North and South polar values
            for ( Itime = 0; Itime < Ntime; Itime++ ) {
                for ( Idepth = 0; Idepth < Ndepth; Idepth++ ) {
                    long double S_sol_Psi = 0.,
                                N_sol_Psi = 0.,
                                S_sol_Phi = 0.,
                                N_sol_Phi = 0.;
                    for ( Ilon = 0; Ilon < Nlon_coarse; Ilon++ ) {
                        size_t Ipt_N = coarsened_grid->local_index( Itime, Idepth, Nlat_coarse-1, Ilon );
                        size_t Ipt_S = coarsened_grid->local_index( Itime, Idepth, 0,             Ilon );

                        N_sol_Psi += soln_Psi[Ipt_N] / Nlon_coarse;
                        N_sol_Phi += soln_Phi[Ipt_N] / Nlon_coarse;

                        S_sol_Psi += soln_Psi[Ipt_S] / Nlon_coarse;
                        S_sol_Phi += soln_Phi[Ipt_S] / Nlon_coarse;
                    }
                    N_poles_Psi_vals[Itime*Ndepth + Idepth] = N_sol_Psi;
                    N_poles_Phi_vals[Itime*Ndepth + Idepth] = N_sol_Phi;
                    S_poles_Psi_vals[Itime*Ndepth + Idepth] = S_sol_Psi;
                    S_poles_Phi_vals[Itime*Ndepth + Idepth] = S_sol_Phi;
                }
            }
            */

            iter_cycle++;
            if ( iter_cycle > Helmholtz_data.iteration_max ) { keep_solving = false; }
            if (constants::DO_TIMING) { timing_records.add_to_record(MPI_Wtime() - clock_on, "Solving"); }

            if ( iter_cycle % Helmholtz_data.iterations_per_cycle == 0 ) {
                // Get convergence metrics

                // Compute velocity components and vort,div
                get_vels_vort_and_div( coarsened_grid, solution_grid, cos_ip, cos_im );

                double mean_vort, mean_div;
                demean( coarsened_grid->variables.at("proj_vort"), coarsened_grid, mean_vort );
                demean( coarsened_grid->variables.at("proj_div" ), coarsened_grid, mean_div  );

                //double mean_Psi, mean_Phi;
                //demean( solution_grid.variables.at("Psi"),  coarsened_grid, mean_Psi );
                //demean( solution_grid.variables.at("Phi" ), coarsened_grid, mean_Phi );

                // Copy the solution over to the other storage set
                coarsened_grid->variables.at("Psi") = solution_grid.variables.at("Psi");
                coarsened_grid->variables.at("Phi") = solution_grid.variables.at("Phi");

                // compute errors and check convergence
                long double vel_2_err  = 0, vel_inf_err  = 0,
                            vort_2_err = 0, vort_inf_err = 0,
                            div_2_err  = 0, div_inf_err  = 0;
                long double vel_2_norm  = 0, vel_inf_norm  = 0,
                            vort_2_norm = 0, vort_inf_norm = 0,
                            div_2_norm  = 0, div_inf_norm  = 0;
                long double vort_2_ener = 0, div_2_ener = 0;
                long double Psi_2_norm = 0, Phi_2_norm = 0;

                get_norms_GS(  vel_2_err,  vort_2_err,  div_2_err,
                            vel_2_norm, vort_2_norm, div_2_norm,
                            vel_inf_err,  vort_inf_err,  div_inf_err,
                            vel_inf_norm, vort_inf_norm, div_inf_norm,
                            vort_2_ener, div_2_ener,
                            Psi_2_norm, Phi_2_norm,
                            coarsened_grid
                          );

                // If we've hit tolerance, halt
                if (      (vort_2_err / vort_2_norm < rel_tol) 
                      and (div_2_err  / div_2_norm  < rel_tol) ) {
                    keep_solving = false;
                }

                // If diverged or stagnated, halt
                if ( vort_2_err / vort_2_norm > 500 * min_vort_err ) {
                    keep_solving = false;
                } else if ( std::fabs( (vort_2_err / vort_2_norm) - min_vort_err ) < std::fmax(1e-4 * min_vort_err, 1e-14) ) {
                    vort_stagnation_count++;
                } else {
                    vort_stagnation_count = 0;
                }
                min_vort_err = std::fmin( min_vort_err, vort_2_err / vort_2_norm );

                // If diverged or stagnated, halt
                if ( div_2_err / div_2_norm > 500 * min_div_err ) {
                    keep_solving = false;
                } else if ( std::fabs( (div_2_err / div_2_norm) - min_div_err ) < std::fmax(1e-4 * min_div_err, 1e-14) ) {
                    div_stagnation_count++;
                } else {
                    div_stagnation_count = 0;
                }
                min_div_err = std::fmin( min_div_err, div_2_err / div_2_norm );

                // If the solver has stagnated, move on to the next resolution tier
                if ( ( div_stagnation_count > 500 ) and ( vort_stagnation_count > 500 ) ) { keep_solving = false; }

                // Save the convergence records
                Helmholtz_data.vel_2_errors.push_back(  (double) vel_2_err  );
                Helmholtz_data.vort_2_errors.push_back( (double) vort_2_err );
                Helmholtz_data.div_2_errors.push_back(  (double) div_2_err  );

                Helmholtz_data.vel_inf_errors.push_back(  (double) vel_inf_err  );
                Helmholtz_data.vort_inf_errors.push_back( (double) vort_inf_err );
                Helmholtz_data.div_inf_errors.push_back(  (double) div_inf_err  );

                // Save the convergence reference values
                Helmholtz_data.vel_2_norms.push_back(  (double) vel_2_norm  );
                Helmholtz_data.vort_2_norms.push_back( (double) vort_2_norm );
                Helmholtz_data.div_2_norms.push_back(  (double) div_2_norm  );

                Helmholtz_data.vel_inf_norms.push_back(  (double) vel_inf_norm  );
                Helmholtz_data.vort_inf_norms.push_back( (double) vort_inf_norm );
                Helmholtz_data.div_inf_norms.push_back(  (double) div_inf_norm  );


                // If sufficient convergence, stop
                if ( Helmholtz_data.IsConverged() ) { 
                    keep_solving = false;
                    fprintf( stdout,"  Halting at solver cycle %'zu. Solver has converged to desired tolerance.\n", iter_cycle );
                }

                #if DEBUG >= 0
                fprintf( stdout,"  Solver cycle %'zu complete. Relative 2-norm errors of vels, vort, and div are %.2e, %.2e%s, and %.2e%s\n", 
                        iter_cycle,
                        Helmholtz_data.vel_2_errors.back()  / Helmholtz_data.vel_2_norms.back(),
                        Helmholtz_data.vort_2_errors.back() / Helmholtz_data.vort_2_norms.back(), 
                        vort_stagnation_count > 0 ? "(*)" : "",
                        Helmholtz_data.div_2_errors.back()  / Helmholtz_data.div_2_norms.back() ,
                        div_stagnation_count > 0 ? "(*)" : ""
                       );
                #endif
            }
        }

        fprintf( stdout,"\nSolver completed with relative 2-norm errors of vels, vort, and div are %.2e, %.2e, and %.2e\n", 
                Helmholtz_data.vel_2_errors.back()  / Helmholtz_data.vel_2_norms.back(),
                Helmholtz_data.vort_2_errors.back() / Helmholtz_data.vort_2_norms.back(), 
                Helmholtz_data.div_2_errors.back()  / Helmholtz_data.div_2_norms.back() 
               );

        if (constants::DO_TIMING) { clock_on = MPI_Wtime(); }

        coarsened_grid->variables.at("Psi") = solution_grid.variables.at("Psi");
        coarsened_grid->variables.at("Phi") = solution_grid.variables.at("Phi");

        // Compute velocity components
        get_vels_vort_and_div( coarsened_grid, *coarsened_grid, cos_ip, cos_im );

        demean( coarsened_grid->variables.at("proj_vort"), coarsened_grid, mean_vort );
        demean( coarsened_grid->variables.at("proj_div" ), coarsened_grid, mean_div  );
        fprintf( stdout, "Final mean values of vorticity and divergence are %.2e, %.2e\n", mean_vort, mean_div );

        coarsened_grid->variables.at("projected_vort") = coarsened_grid->variables.at("proj_vort");
        coarsened_grid->variables.at("projected_div") = coarsened_grid->variables.at("proj_div");

        #pragma omp parallel default(none) private(Ipt) \
        shared( coarsened_grid ) firstprivate( Npts_coarse )
        {
            #pragma omp for collapse(1) schedule(static)
            for ( size_t Ipt = 0; Ipt < Npts_coarse; Ipt++ ) {
                coarsened_grid->variables.at("projected_u_lon")[Ipt] = 
                      coarsened_grid->variables.at("proj_uo_tor")[Ipt] 
                    + coarsened_grid->variables.at("proj_uo_pot")[Ipt];
                coarsened_grid->variables.at("projected_u_lat")[Ipt] = 
                      coarsened_grid->variables.at("proj_vo_tor")[Ipt] 
                    + coarsened_grid->variables.at("proj_vo_pot")[Ipt];
            }
        }
        if (constants::DO_TIMING) { timing_records.add_to_record(MPI_Wtime() - clock_on, "ProcessingSolution"); }

        // As a sanity check, get the divergence of the toroidal field and vorticity of the potential field
        get_vort_div_of_pot_tor( coarsened_grid, cos_ip, cos_im );
        
        /*
        // Vorticity
        curl_vec_dot_er( 
                coarsened_grid->variables.at("proj_vort_of_pot"),
                coarsened_grid->variables.at("proj_uo_pot"),
                coarsened_grid->variables.at("proj_vo_pot"),
                *coarsened_grid,
                coarsened_grid->mask
                );

        // Divergence
        div_vec( 
                coarsened_grid->variables.at("proj_div_of_tor"),
                coarsened_grid->variables.at("proj_uo_tor"),
                coarsened_grid->variables.at("proj_vo_tor"),
                *coarsened_grid,
                coarsened_grid->mask
               );
        */

        std::ostringstream name_stream;
        name_stream << "projection_R" << (refine_level) << ".nc";
        const std::string output_filename = name_stream.str();

        // Put the mask back
        //coarsened_grid->mask = saved_mask;

        write_Helmholtz_output( output_filename, Helmholtz_data, *coarsened_grid, solution_grid );
        std::fill( coarsened_grid->mask.begin(), coarsened_grid->mask.end(), true );

        // If we're doing timings, then print out and reset values now
        if (constants::DO_TIMING) { 
            timing_records.print();
            timing_records.reset();
            fflush(stdout);
        }

    }
}
