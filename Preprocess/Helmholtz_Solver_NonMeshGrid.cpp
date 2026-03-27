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
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/QR>

// Pragma-magic to allow reduction over vector operator
//      thanks to: https://stackoverflow.com/questions/43168661/openmp-and-reduction-on-stdvector
#pragma omp declare reduction(vec_long_plus : std::vector<long double> : \
                              std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<long double>())) \
                    initializer(omp_priv = decltype(omp_orig)(omp_orig.size()))

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

void filter_PsiPhi_point(
        double & Psi_bar,
        double & Phi_bar,
        const size_t target_index,
        std::vector<double> & Psi,
        std::vector<double> & Phi,
        dataset * coarsened_grid,
        const double filter_scale
        ){

    const size_t num_pts = Psi.size();
    const size_t num_neighbours = coarsened_grid->num_neighbours;

    // Want to use bitset directly, but does not accept runtime-specified size
    // so will use std::vector<bool>, which is a bitset under the hood
    // Note that each element is 1 bit, not byte (because of the 
    // specialization of vector<bool>), so that, while e.g. LLC4320 requires
    // ~1.8GB to store a single field at double precision, these bitsets 
    // only require <30MB. So, having each thread / processor keep their own
    // copy shouldn't be too onerous memory-wise
    std::deque<size_t> points_to_test;
    std::vector<bool> was_rejected(num_pts, false),
                      was_accepted(num_pts, false),
                      planned_for_testing(num_pts, false);

    #ifdef NDEBUG
    const double target_lat = coarsened_grid->latitude[  target_index ];
    const double target_lon = coarsened_grid->longitude[ target_index ];
    double dArea = coarsened_grid->areas[target_index];
    #else
    const double target_lat = coarsened_grid->latitude.at(  target_index );
    const double target_lon = coarsened_grid->longitude.at( target_index );
    double dArea = coarsened_grid->areas.at(target_index);
    #endif

    double kernel_normalization = 0.0;
    kernel_normalization += 1. * dArea;

    // Set up the filter value with the kernel centre
    // assumed to always have weight of 1
    double local_Psi_val = Psi[target_index],
           local_Phi_val = Phi[target_index],
           filter_Psi_value = dArea * local_Psi_val,
           filter_Phi_value = dArea * local_Phi_val;

    // Next, seed the 'points to test' with the adjacent points
    points_to_test.clear();
    for ( size_t neighbour_index = 0; neighbour_index < num_neighbours; neighbour_index++ ) {
        #ifdef NDEBUG
        points_to_test.push_back( 
                coarsened_grid->adjacency_indices[target_index][neighbour_index] 
                );
        planned_for_testing[ 
            coarsened_grid->adjacency_indices[target_index][neighbour_index] 
        ] = true;
        #else
        points_to_test.push_back( 
                coarsened_grid->adjacency_indices.at(target_index)[neighbour_index] 
                );
        planned_for_testing[ 
            coarsened_grid->adjacency_indices.at(target_index)[neighbour_index] 
        ] = true;
        #endif
    }

    // So long as we still have points to test, keep testing!
    // This implements a breadth-first search through the adjacency matrix to build
    //      filtering kernel. If a point is within distance, add it to the kernel,
    //      and then test it's neighbours. If those are in, test their neighbours,
    //      and so on. If a point is too far away, we do not test it's neighbours.
    // This then assumes that for any point X within distance L of Y, that there
    //      is an adjacency path from Y to X strictly using points within distance
    //      L of Y.
    // Along the way, test for double inclusing to make sure that we don't test 
    //      points repeatedly etc.
    while ( points_to_test.size() > 0 ) {

        // Pull out the most-recently-added point, and remove it from the 'to test' list,
        // since we're testing it now.
        // Since we're pulling out the most-recently-added, that effectively makes this a
        // depth-first search to build the kernel.
        size_t search_index = points_to_test.front();
        points_to_test.pop_front();
        planned_for_testing[ search_index ] = false;

        double local_dist = coarsened_grid->distance_from_indices( target_index, search_index );

        double kern_val = kernel( local_dist, filter_scale );

        if ( ( kern_val > 1e-10 ) or ( local_dist <= 0.5 * filter_scale )  ) {
            was_accepted[search_index] = true;

            // Accumulate the coarse values
            #ifdef NDEBUG
            dArea = coarsened_grid->areas[search_index];
            local_Psi_val = Psi[search_index];
            local_Phi_val = Phi[search_index];
            #else
            dArea = coarsened_grid->areas.at(search_index);
            local_Psi_val = Psi.at(search_index);
            local_Phi_val = Phi.at(search_index);
            #endif
            double kA = kern_val * dArea;
            kernel_normalization += kA;
            filter_Psi_value += kA * local_Psi_val;
            filter_Phi_value += kA * local_Phi_val;

            // Since this point is in the kernel
            // add its neighbours to the 'to search' list
            for ( int ii = 0; ii < num_neighbours; ii++ ) {
                #ifdef NDEBUG
                size_t neighbour_index = coarsened_grid->adjacency_indices[search_index][ii];
                #else
                size_t neighbour_index = coarsened_grid->adjacency_indices.at(search_index).at(ii);
                #endif

                // but first check if that neighbour has been rejected already
                if ( was_rejected[neighbour_index] ) { continue; }

                // then check if that neighbour is already accepted
                if ( was_accepted[neighbour_index] ) { continue; }

                // then check if that neighbour is already on the search list
                if ( planned_for_testing[neighbour_index] ) { continue; }

                // if it's not on any of those list already
                // then add it to the 'points to test'
                points_to_test.push_back( neighbour_index );
                planned_for_testing[neighbour_index] = true;

            }
        } else {
            // Otherwise, record this point as rejected, and move on
            was_rejected[search_index] = true;
        }

    }

   Psi_bar = (kernel_normalization == 0) ? 0. : (filter_Psi_value / kernel_normalization);
   Phi_bar = (kernel_normalization == 0) ? 0. : (filter_Phi_value / kernel_normalization);

}

void get_true_vort_div_fast(
        dataset * coarsened_grid,
        dataset & PsiPhi_grid
            ) {

    // Similar to get_vels_vort_and_div, except this manually implements the derivatives
    // By lumping all of the derivative calculations together we can be more efficient, especially
    // with the threading.

    const std::vector<double>   &latitude   = coarsened_grid->latitude,
                                &longitude  = coarsened_grid->longitude;
    const std::vector<int>  &myCounts = coarsened_grid->myCounts;

    const int   Ntime  = myCounts.at(0),
                Ndepth = myCounts.at(1),
                Nlat   = (constants::GRID_TYPE == constants::GridType::MeshGrid) ? myCounts.at(2) : 1,
                Nlon   = (constants::GRID_TYPE == constants::GridType::MeshGrid) ? myCounts.at(3) : 1;
    const size_t Npts  = (constants::GRID_TYPE == constants::GridType::MeshGrid) ? Nlat*Nlon : latitude.size();


    size_t Ipt, In;

    // zero out the fields
    #pragma omp parallel default(none) private( Ipt ) \
    shared( coarsened_grid ) firstprivate( Npts )
    {
        #pragma omp for collapse(1) schedule(static)
        for ( Ipt = 0; Ipt < Npts; Ipt++ ) {

            coarsened_grid->variables.at("vort")[Ipt] = 0;
            coarsened_grid->variables.at("div")[Ipt]  = 0;
        }
    }

    //
    // loop to get vorticity and divergence
    //

    #pragma omp parallel default(none) private( Ipt, In ) \
    firstprivate( Ntime, Ndepth, Nlat, Nlon, Npts ) \
    shared( coarsened_grid, PsiPhi_grid )
    {

        #pragma omp for collapse(1) schedule(static)
        for ( Ipt = 0; Ipt < Npts; Ipt++ ) {

            double vort_tmp = 0, div_tmp = 0;

            const double uo_ref = coarsened_grid->variables.at("u_lon")[Ipt],
                         vo_ref = coarsened_grid->variables.at("u_lat")[Ipt];

            const std::vector<double> grad_weights_lon = coarsened_grid->grad_weights_lon[Ipt],
                                      grad_weights_lat = coarsened_grid->grad_weights_lat[Ipt],
                                      perp_grad_weights_lon = coarsened_grid->perp_grad_weights_lon[Ipt],
                                      perp_grad_weights_lat = coarsened_grid->perp_grad_weights_lat[Ipt];

            for ( In = 0; In < coarsened_grid->num_neighbours; In++ ) {
                /*
                size_t Jpt = coarsened_grid->adjacency_indices[Ipt][In];

                const double uo_val = coarsened_grid->variables.at("u_lon")[Jpt];
                const double vo_val = coarsened_grid->variables.at("u_lat")[Jpt];
                div_tmp +=   grad_weights_lon[In] * ( uo_val + uo_ref ) * 0.5
                           + grad_weights_lat[In] * ( vo_val + vo_ref ) * 0.5;

                vort_tmp +=   perp_grad_weights_lon[In] * ( uo_val + uo_ref ) * 0.5
                            + perp_grad_weights_lat[In] * ( vo_val + vo_ref ) * 0.5;
                */
                div_tmp +=   grad_weights_lon[In] * uo_ref
                           + grad_weights_lat[In] * vo_ref;

                vort_tmp +=   perp_grad_weights_lon[In] uo_ref
                            + perp_grad_weights_lat[In] vo_ref;

            }
            coarsened_grid->variables["div" ][Ipt] = div_tmp;
            coarsened_grid->variables["vort"][Ipt] = vort_tmp;
        }
    }


}

void get_vels_vort_and_div_fast(
        dataset * coarsened_grid,
        dataset & PsiPhi_grid
            ) {

    // Similar to get_vels_vort_and_div, except this manually implements the derivatives
    // By lumping all of the derivative calculations together we can be more efficient, especially
    // with the threading.


    const std::vector<double>   &latitude   = coarsened_grid->latitude,
                                &longitude  = coarsened_grid->longitude;
    const std::vector<int>  &myCounts = coarsened_grid->myCounts;

    const int   Ntime  = myCounts.at(0),
                Ndepth = myCounts.at(1),
                Nlat   = (constants::GRID_TYPE == constants::GridType::MeshGrid) ? myCounts.at(2) : 1,
                Nlon   = (constants::GRID_TYPE == constants::GridType::MeshGrid) ? myCounts.at(3) : 1;
    const size_t Npts  = (constants::GRID_TYPE == constants::GridType::MeshGrid) ? Nlat*Nlon : latitude.size();


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

    //const double R_E_inv = 1. / constants::R_earth;


    //
    // First loop set gets the velocity components
    //

    #pragma omp parallel default(none) private( Ipt, In ) \
    firstprivate( Ntime, Ndepth, Nlat, Nlon, Npts ) \
    shared( coarsened_grid, PsiPhi_grid, Psi, Phi, \
         proj_uo_tor, proj_vo_tor, proj_uo_pot, proj_vo_pot )
    {

        #pragma omp for collapse(1) schedule(static)
        for ( Ipt = 0; Ipt < Npts; Ipt++ ) {

            double uo_tor_tmp = 0, vo_tor_tmp = 0, uo_pot_tmp = 0, vo_pot_tmp = 0;

            const double Psi_ref = Psi[Ipt];
            const double Phi_ref = Phi[Ipt];

            const std::vector<double> grad_weights_lon = coarsened_grid->grad_weights_lon[Ipt],
                                      grad_weights_lat = coarsened_grid->grad_weights_lat[Ipt],
                                      perp_grad_weights_lon = coarsened_grid->perp_grad_weights_lon[Ipt],
                                      perp_grad_weights_lat = coarsened_grid->perp_grad_weights_lat[Ipt];

            for ( In = 0; In < coarsened_grid->num_neighbours; In++ ) {
                size_t Jpt = coarsened_grid->adjacency_indices[Ipt][In];


                #ifdef NDEBUG
                const double Psi_loc = Psi[Jpt],
                             Phi_loc = Phi[Jpt];
                #else
                const double Psi_loc = Psi.at(Jpt),
                             Phi_loc = Phi.at(Jpt);
                #endif

                uo_tor_tmp += perp_grad_weights_lon[In] * ( Psi_loc - Psi_ref );
                vo_tor_tmp += perp_grad_weights_lat[In] * ( Psi_loc - Psi_ref );

                uo_pot_tmp += grad_weights_lon[In] * ( Phi_loc - Phi_ref );
                vo_pot_tmp += grad_weights_lat[In] * ( Phi_loc - Phi_ref );

            }
            #ifdef NDEBUG
            proj_uo_tor[Ipt] = uo_tor_tmp;
            proj_vo_tor[Ipt] = vo_tor_tmp;

            proj_uo_pot[Ipt] = uo_pot_tmp;
            proj_vo_pot[Ipt] = vo_pot_tmp;
            #else
            proj_uo_tor.at(Ipt) = uo_tor_tmp;
            proj_vo_tor.at(Ipt) = vo_tor_tmp;

            proj_uo_pot.at(Ipt) = uo_pot_tmp;
            proj_vo_pot.at(Ipt) = vo_pot_tmp;
            #endif
        }
    }


    //
    // Second loop sets gets vorticity and divergence
    //

    #pragma omp parallel default(none) private( Ipt, In ) \
    firstprivate( Ntime, Ndepth, Nlat, Nlon, Npts ) \
    shared( coarsened_grid, PsiPhi_grid, proj_vort, proj_div, proj_uo_tor, \
           proj_vo_tor, proj_uo_pot, proj_vo_pot )
    {

        #pragma omp for collapse(1) schedule(static)
        for ( Ipt = 0; Ipt < Npts; Ipt++ ) {
            double cos_lat = 
                coarsened_grid->sin_cos_lat_pre_comp ?
                coarsened_grid->cos_lat[Ipt] :
                cos( coarsened_grid->latitude[Ipt] );
            const double R_E_cos_lat_inv = 1. / ( constants::R_earth * cos_lat );

            double vort_tmp = 0, div_tmp = 0;

            const double uo_tor_ref = proj_uo_tor[Ipt],
                         vo_tor_ref = proj_vo_tor[Ipt],
                         uo_pot_ref = proj_uo_pot[Ipt],
                         vo_pot_ref = proj_vo_pot[Ipt];

            const std::vector<double> grad_weights_lon = coarsened_grid->grad_weights_lon[Ipt],
                                      grad_weights_lat = coarsened_grid->grad_weights_lat[Ipt],
                                      perp_grad_weights_lon = coarsened_grid->perp_grad_weights_lon[Ipt],
                                      perp_grad_weights_lat = coarsened_grid->perp_grad_weights_lat[Ipt];

            for ( In = 0; In < coarsened_grid->num_neighbours; In++ ) {
                /*
                size_t Jpt = coarsened_grid->adjacency_indices[Ipt][In];
                div_tmp +=   grad_weights_lon[In] * ( proj_uo_pot[Jpt] + uo_pot_ref ) * 0.5
                           + grad_weights_lat[In] * ( proj_vo_pot[Jpt] + vo_pot_ref ) * 0.5;

                vort_tmp +=   perp_grad_weights_lon[In] * ( proj_uo_tor[Jpt] + uo_tor_ref ) * 0.5
                            + perp_grad_weights_lat[In] * ( proj_vo_tor[Jpt] + vo_tor_ref ) * 0.5;
                */
                div_tmp +=   grad_weights_lon[In] * uo_pot_ref
                           + grad_weights_lat[In] * vo_pot_ref;

                vort_tmp +=   perp_grad_weights_lon[In] * uo_tor_ref
                            + perp_grad_weights_lat[In] * vo_tor_ref;

            }
            proj_vort[Ipt] = vort_tmp;
            proj_div[ Ipt] = div_tmp;
        }
    }

}

double get_mean( 
            std::vector<double> & f,
            dataset * data
             ) {


    size_t Ipt, Npts = f.size();
    long double f_sum = 0., area_sum = 0.;
    #pragma omp parallel default(none) \
    shared( f, data ) private( Ipt ) firstprivate( Npts ) \
    reduction( +:f_sum,area_sum )
    {
        #pragma omp for collapse(1) schedule(static)
        for ( Ipt = 0; Ipt < Npts; Ipt++) {
            long double dA = data->areas[Ipt];
            area_sum += dA;
            f_sum  += dA * f[Ipt];
        }
    }
    double mean_val = (double) (f_sum / area_sum);
    return mean_val;

}

void demean( std::vector<double> & f,
             dataset * data,
             double & mean_val,
             std::vector<short int> * saved_mask
             ) {


    size_t Ipt, Npts = f.size();
    mean_val = get_mean( f, data );
    #pragma omp parallel default(none) \
    shared( f, data, saved_mask ) private( Ipt ) \
    firstprivate( Npts, mean_val )
    {
    #pragma omp for collapse(1) schedule(dynamic,16)
        for ( Ipt = 0; Ipt < Npts; Ipt++) {
            if ( (saved_mask == NULL) or (saved_mask->at(Ipt)) ) {
                f[Ipt] -= mean_val;
            }
        }
    }

}

void get_norms(
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
    bool use_vel = true, use_vort_div = true;
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

    #pragma omp parallel default(none) \
    shared( coarsened_grid, u_lon, u_lat, vort, div, Phi, Psi, \
         proj_uo_tor, proj_vo_tor, proj_uo_pot, proj_vo_pot, proj_vort, proj_div ) \
    private( Ipt ) \
    firstprivate( Npts_coarse, use_vel, use_vort_div ) \
    reduction( +:total_area,vel_2_err_L,vort_2_err_L,div_2_err_L,vel_2_norm_L,vort_2_norm_L,div_2_norm_L ) \
    reduction( +:vort_2_ener_L,div_2_ener_L,Phi_2_norm_L,Psi_2_norm_L )\
    reduction( max:vel_inf_err,vort_inf_err,div_inf_err,vel_inf_norm,vort_inf_norm,div_inf_norm )
    {
        #pragma omp for collapse(1) schedule(static)
        for ( Ipt = 0; Ipt < Npts_coarse; Ipt++ ) {

            long double area = coarsened_grid->areas[Ipt];
            long double incr;

            total_area += area;

            if ( use_vel ) {
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
            }

            if ( use_vort_div ) {

                // Vorticity errors
                #ifdef NDEBUG
                incr = std::fabs( proj_vort[Ipt] - vort[Ipt] );
                #else
                incr = std::fabs(
                          coarsened_grid->variables.at("proj_vort").at(Ipt)
                        - coarsened_grid->variables.at("vort").at(Ipt)
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
            }

            // Reference norms
            if ( use_vel ) {
                #ifdef NDEBUG
                incr = sqrt( pow( u_lon[Ipt], 2) + pow( u_lat[Ipt], 2) );
                #else
                incr = sqrt(   pow(  coarsened_grid->variables.at("u_lon").at(Ipt), 2)
                             + pow(  coarsened_grid->variables.at("u_lat").at(Ipt), 2)
                        );
                #endif
                vel_inf_norm = std::fmax( incr, vel_inf_norm );
                vel_2_norm_L += area * pow( incr, 2 ); 
            }

            if ( use_vort_div ) {
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
            }

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


void Helmholtz_Solver_NonMeshGrid(
        const std::string output_fname,
        dataset & source_data,
        const double rel_tol,
        const unsigned int max_iters,
        const unsigned int iters_per_batch,
        const bool weight_err,
        const bool use_mask,
        const bool use_vel,
        const bool use_vort_div,
        const bool collapse_land,
        const int num_refinements,
        const double CFL,
        const bool seed_provided,
        const std::vector<double> & Psi_seed,
        const std::vector<double> & Phi_seed,
        const MPI_Comm comm
        ) {

    if ( not(use_vort_div and use_vel) ) {
        throw std::runtime_error("Must use both vortdiv and vels for time-step solver.");
    }

    int wRank, wSize;
    MPI_Comm_rank( comm, &wRank );
    MPI_Comm_size( comm, &wSize );

    // If we've passed the DO_TIMING flag, then create some timing vars
    Timing_Records timing_records;
    double clock_on;

    const size_t Npts = source_data.mask.size();

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

    #if DEBUG >= 2
    fprintf(stdout, "Initializing class objects for data storage.\n");
    #endif
    HelmholtzDataClass Helmholtz_data;
    Helmholtz_data.collapse_land = collapse_land;
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
            Nlat_coarse = floor(sqrt( 0.5 * (Npts / pow(4., 2*(refine_level) )) ));
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

            if ( constants::GRID_TYPE == constants::GridType::MeshGrid ) {
                Npts_coarse = (size_t) Nlat_coarse * Nlon_coarse;
            } else {
                Npts_coarse = coarsened_grid->longitude.size(); // The coarsened grid might have a different size
            }
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

            coarsened_grid->variables.insert( std::pair< std::string, std::vector<double> >(
                        "vort", std::vector<double>(Npts_coarse, 0.) ) );
            coarsened_grid->variables.insert( std::pair< std::string, std::vector<double> >(
                        "div", std::vector<double>(Npts_coarse, 0.) ) );

        }

        saved_mask = coarsened_grid->mask; // store a copy

        // In addition to cell areas, we also need cell side lengths
        // this is necessary for cells with wild aspect ratios
        // such as in NEMO tripolar grid near the tripole
        std::vector<double> min_side_length( Npts_coarse, 0. ), diff_coeff( Npts_coarse, 0. );
        std::map< size_t, std::vector<size_t> > local_kernel_indices;
        std::map< size_t, std::vector<double> > local_kernel_kA;
        double min_da = 4 * M_PI * pow( constants::R_earth, 2.);;
        #pragma omp parallel default(none) private( Ipt ) \
        shared(min_side_length, coarsened_grid, diff_coeff, saved_mask ) \
        firstprivate(Npts_coarse, refine_level) \
        reduction( min:min_da )
        {
            #pragma omp for collapse(1) schedule(guided) 
            for ( Ipt = 0; Ipt < Npts_coarse; Ipt++ ) {

                if ( constants::GRID_TYPE == constants::GridType::MeshGrid ) {
                    min_side_length[Ipt] = sqrt( coarsened_grid->areas[Ipt] );
                } else {

                    min_side_length[Ipt] = sqrt( coarsened_grid->areas[Ipt] );

                    //min_da = std::fmin( min_da, coarsened_grid->areas[Ipt]);

                    double sum_wij_mags = 0.;
                    for ( size_t In = 0; In < coarsened_grid->num_neighbours; In++ ) {
                        sum_wij_mags += coarsened_grid->grad_weights_lon[Ipt][In] * coarsened_grid->grad_weights_lon[Ipt][In]
                            + coarsened_grid->grad_weights_lat[Ipt][In] * coarsened_grid->grad_weights_lat[Ipt][In];
                    }
                    min_da = std::fmin( min_da, 1. / sum_wij_mags );

                }
                diff_coeff[Ipt] = pow(min_side_length[Ipt], 2.);
            }
        }
        fprintf( stdout, "min_da = %.3e\n", min_da );

        if ( not(constants::GRID_TYPE == constants::GridType::MeshGrid) ) {
            std::ostringstream adj_stream;
            adj_stream << "adjacency_R" << (refine_level) << ".nc";
            coarsened_grid->write_adjacency( adj_stream.str() );
        }

        // We're done with the land mask, so remove it. This causes vort and
        // div to be computed on land cells to ensure agreement
        std::fill( coarsened_grid->mask.begin(), coarsened_grid->mask.end(), true );

        // Get vorticity and divergence
        #if DEBUG >= 2
        fprintf(stdout, "Computing vorticity and divergence on coarse grid.\n");
        #endif
        if (constants::DO_TIMING) { clock_on = MPI_Wtime(); }

        get_true_vort_div_fast( coarsened_grid, *coarsened_grid );

        double mean_vort, mean_div;
        //demean( coarsened_grid->variables.at("vort"), coarsened_grid, mean_vort, refine_level == 0 ? &saved_mask : NULL );
        //demean( coarsened_grid->variables.at("div" ), coarsened_grid, mean_div,  refine_level == 0 ? &saved_mask : NULL );
        demean( coarsened_grid->variables.at("vort"), coarsened_grid, mean_vort );
        demean( coarsened_grid->variables.at("div" ), coarsened_grid, mean_div  );
        fprintf( stdout, "Mean values of vorticity and divergence are %.2e, %.2e\n", mean_vort, mean_div );

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
        std::fill( solution_grid.mask.begin(), solution_grid.mask.end(), true ); 

        solution_grid.variables.insert( std::pair< std::string, std::vector<double> >(
                    "Psi", std::vector<double>(Npts_coarse, 0.) ) );
        solution_grid.variables.insert( std::pair< std::string, std::vector<double> >(
                    "Phi", std::vector<double>(Npts_coarse, 0.) ) );


        // Reset the Helmholtz data for this resolution
        Helmholtz_data.clear();
        //Helmholtz_data.collapse_land = collapse_land;
        Helmholtz_data.use_vort_div = use_vort_div;
        Helmholtz_data.use_vel      = use_vel;
        Helmholtz_data.weight_err   = weight_err;
        Helmholtz_data.tolerance    = rel_tol;

        Helmholtz_data.iterations_per_cycle = refine_level == 0 ?   5000 : Npts_coarse * 1;
        Helmholtz_data.iteration_max        = refine_level == 0 ? 500000 : Npts_coarse * 1000;

        // Use spherical harmonics to set the seed
        coarsened_grid->variables.insert( std::pair< std::string, std::vector<double> >(
                    "Psi", std::vector<double>(Npts_coarse, 0.) ) );
        coarsened_grid->variables.insert( std::pair< std::string, std::vector<double> >(
                    "Phi", std::vector<double>(Npts_coarse, 0.) ) );
        if ( (num_refinements == 0) and (seed_provided) ) {
            coarsened_grid->variables["Psi"] = Psi_seed;
            coarsened_grid->variables["Phi"] = Phi_seed;
        } else if ( refine_level == num_refinements ) {
            if ( not(constants::GRID_TYPE == constants::GridType::MeshGrid) ) {
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
        get_vels_vort_and_div_fast( coarsened_grid, *coarsened_grid );
        //demean( coarsened_grid->variables.at("proj_vort"), coarsened_grid, mean_vort, refine_level == 0 ? &saved_mask : NULL );
        //demean( coarsened_grid->variables.at("proj_div" ), coarsened_grid, mean_div,  refine_level == 0 ? &saved_mask : NULL );
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

        get_norms(  vel_2_err,  vort_2_err,  div_2_err,
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

        const bool common_dt = true;
        if ( common_dt ) {
            fprintf( stdout, "\n\nCommon dt at all points\n\n\n" );
        } else {
            fprintf( stdout, "\n\ndt is spatially varying\n\n\n" );
        }

        const double Lap_scale = 5. * sqrt( 4 * M_PI * pow(constants::R_earth, 2.) / Npts_coarse );

        while ( keep_solving ) {

            const double dt = CFL;

            if (constants::DO_TIMING) { clock_on = MPI_Wtime(); }

            #pragma omp parallel default(none) private(Ipt) \
            shared( Helmholtz_data, coarsened_grid, solution_grid, dPsidt, dPhidt, \
                    soln_Psi, soln_Phi, cgrid_proj_vort, cgrid_proj_div, cgrid_vort, cgrid_div, \
                    diff_coeff, saved_mask, local_kernel_indices, local_kernel_kA ) \
            firstprivate( Npts_coarse, dt, refine_level, CFL, min_da )
            {
                #pragma omp for collapse(1) schedule(static)
                for ( Ipt = 0; Ipt < Npts_coarse; Ipt++ ) {

                    double Psi_force, Phi_force, diff_coeff_local;
                    diff_coeff_local = common_dt ? min_da : coarsened_grid->areas[Ipt];
                    #ifdef NDEBUG
                    Psi_force = (cgrid_proj_vort[Ipt] - cgrid_vort[Ipt]);
                    Phi_force = (cgrid_proj_div[ Ipt] - cgrid_div[ Ipt]);
                    #else
                    Psi_force = (cgrid_proj_vort.at(Ipt) - cgrid_vort.at(Ipt));
                    Phi_force = (cgrid_proj_div.at( Ipt) - cgrid_div.at( Ipt));
                    #endif

                    // d(Psi)/dt = Lap(Psi) - vorticity
                    dPsidt[Ipt] = diff_coeff_local * Psi_force;

                    // d(Phi)/dt = Lap(Phi) - divergence
                    dPhidt[Ipt] = diff_coeff_local * Phi_force;

                    // Euler for Diffusion
                    #ifdef NDEBUG
                    soln_Psi[Ipt] += dt * dPsidt[Ipt];
                    soln_Phi[Ipt] += dt * dPhidt[Ipt];
                    #else
                    soln_Psi.at(Ipt) += dt * dPsidt.at(Ipt);
                    soln_Phi.at(Ipt) += dt * dPhidt.at(Ipt);
                    #endif
                }
            }

            // Compute velocity components and vort,div
            // To measure errors and use for next iteration
            get_vels_vort_and_div_fast( coarsened_grid, solution_grid );

            double mean_vort, mean_div;
            double mean_Psi, mean_Phi;
            demean( solution_grid.variables.at("Psi"),  coarsened_grid, mean_Psi );
            demean( solution_grid.variables.at("Phi" ), coarsened_grid, mean_Phi );


            iter_cycle++;
            if ( iter_cycle > Helmholtz_data.iteration_max ) { keep_solving = false; }
            if (constants::DO_TIMING) { timing_records.add_to_record(MPI_Wtime() - clock_on, "Solving"); }

            if ( iter_cycle % Helmholtz_data.iterations_per_cycle == 0 ) {

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

                get_norms(  vel_2_err,  vort_2_err,  div_2_err,
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
                } else if ( std::fabs( (vort_2_err / vort_2_norm) - min_vort_err ) < std::fmax(1e-5 * min_vort_err, 1e-14) ) {
                    vort_stagnation_count++;
                } else {
                    vort_stagnation_count = 0;
                }
                min_vort_err = std::fmin( min_vort_err, vort_2_err / vort_2_norm );

                // If diverged or stagnated, halt
                if ( div_2_err / div_2_norm > 500 * min_div_err ) {
                    keep_solving = false;
                } else if ( std::fabs( (div_2_err / div_2_norm) - min_div_err ) < std::fmax(1e-5 * min_div_err, 1e-14) ) {
                    div_stagnation_count++;
                } else {
                    div_stagnation_count = 0;
                }
                min_div_err = std::fmin( min_div_err, div_2_err / div_2_norm );

                // If the solver has stagnated, move on to the next resolution tier
                if ( ( div_stagnation_count > 5 ) and ( vort_stagnation_count > 5 ) ) { keep_solving = false; }

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
                if (use_vort_div and use_vel) {
                    fprintf( stdout,"  Solver cycle %'zu complete. Relative 2-norm errors of vels, vort, and div are %.2e, %.2e, and %.2e\n", 
                            iter_cycle,
                            Helmholtz_data.vel_2_errors.back()  / Helmholtz_data.vel_2_norms.back(),
                            Helmholtz_data.vort_2_errors.back() / Helmholtz_data.vort_2_norms.back(), 
                            Helmholtz_data.div_2_errors.back()  / Helmholtz_data.div_2_norms.back() 
                           );
                } else if (use_vort_div) {
                    fprintf( stdout,"  Solver cycle %'zu complete. Relative 2-norm errors of vort and div are %.2e and %.2e\n", 
                            iter_cycle,
                            Helmholtz_data.vort_2_errors.back() / Helmholtz_data.vort_2_norms.back(), 
                            Helmholtz_data.div_2_errors.back()  / Helmholtz_data.div_2_norms.back() 
                           );
                } else {
                    fprintf( stdout,"  Solver cycle %'zu complete. Relative 2-norm errors of vels is %.2e\n", 
                            iter_cycle,
                            Helmholtz_data.vel_2_errors.back()  / Helmholtz_data.vel_2_norms.back()
                           );
                }
                #endif
            }
        }

        if (use_vort_div and use_vel) {
            fprintf( stdout,"\nSolver completed with relative 2-norm errors of vels, vort, and div are %.2e, %.2e, and %.2e\n", 
                    Helmholtz_data.vel_2_errors.back()  / Helmholtz_data.vel_2_norms.back(),
                    Helmholtz_data.vort_2_errors.back() / Helmholtz_data.vort_2_norms.back(), 
                    Helmholtz_data.div_2_errors.back()  / Helmholtz_data.div_2_norms.back() 
                   );
        } else if (use_vort_div) {
            fprintf( stdout,"\nSolver completed with relative 2-norm errors of vort and div are %.2e and %.2e\n", 
                    Helmholtz_data.vort_2_errors.back() / Helmholtz_data.vort_2_norms.back(), 
                    Helmholtz_data.div_2_errors.back()  / Helmholtz_data.div_2_norms.back() 
                   );
        } else {
            fprintf( stdout,"\nSolver completed with relative 2-norm errors of vels is %.2e\n", 
                    Helmholtz_data.vel_2_errors.back()  / Helmholtz_data.vel_2_norms.back()
                   );
        }

        if (constants::DO_TIMING) { clock_on = MPI_Wtime(); }

        coarsened_grid->variables.at("Psi") = solution_grid.variables.at("Psi");
        coarsened_grid->variables.at("Phi") = solution_grid.variables.at("Phi");

        // Compute velocity components
        get_vels_vort_and_div_fast( coarsened_grid, *coarsened_grid );

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

        std::ostringstream name_stream;
        name_stream << "projection_R" << (refine_level) << ".nc";
        const std::string output_filename = name_stream.str();

        // Put the mask back
        //coarsened_grid->mask = saved_mask;

        write_Helmholtz_output( output_filename, Helmholtz_data, *coarsened_grid, solution_grid );
        add_attr_to_file( "CFL", CFL, output_filename );
        std::fill( coarsened_grid->mask.begin(), coarsened_grid->mask.end(), true );

        // If we're doing timings, then print out and reset values now
        if (constants::DO_TIMING) { 
            timing_records.print();
            timing_records.reset();
            fflush(stdout);
        }

    }
}
