#include <math.h>
#include <algorithm>
#include <vector>
#include <omp.h>
#include <mpi.h>

#include "../functions.hpp"
#include "../constants.hpp"
#include "../differentiation_tools.hpp"

void compute_KE_spectra_and_slopes( 
        std::vector<double> & u_spectrum_tot, 
        std::vector<double> & u_spectrum_tor, 
        std::vector<double> & u_spectrum_pot,
        std::vector<double> & v_spectrum_tot, 
        std::vector<double> & v_spectrum_tor, 
        std::vector<double> & v_spectrum_pot,
        std::vector<double> & spec_slope_tot, 
        std::vector<double> & spec_slope_tor, 
        std::vector<double> & spec_slope_pot,
        std::vector<double> & u_spectrum_tot_HighPass, 
        std::vector<double> & u_spectrum_tor_HighPass, 
        std::vector<double> & u_spectrum_pot_HighPass,
        std::vector<double> & v_spectrum_tot_HighPass, 
        std::vector<double> & v_spectrum_tor_HighPass, 
        std::vector<double> & v_spectrum_pot_HighPass,
        std::vector<double> & spec_slope_tot_HighPass, 
        std::vector<double> & spec_slope_tor_HighPass, 
        std::vector<double> & spec_slope_pot_HighPass,
        const std::vector<double> & u_lon_tot, 
        const std::vector<double> & u_lon_tor, 
        const std::vector<double> & u_lon_pot,
        const std::vector<double> & u_lat_tot, 
        const std::vector<double> & u_lat_tor, 
        const std::vector<double> & u_lat_pot,
        const std::vector<double> & u_lon_tor_unfilt, 
        const std::vector<double> & u_lon_pot_unfilt,
        const std::vector<double> & u_lat_tor_unfilt, 
        const std::vector<double> & u_lat_pot_unfilt,
        const std::vector<double> & dl_coarse_Phi, 
        const std::vector<double> & dl_coarse_Psi,
        const std::vector<double> & dll_coarse_Phi, 
        const std::vector<double> & dll_coarse_Psi,
        const dataset & source_data,
        const double filter_scale
        ) {

    // Create some tidy names for variables
    const std::vector<short int> &mask = source_data.mask;

    const size_t num_pts = mask.size();
    size_t index;

    // Create some (temporary) storage arrays to hold the ell-derivatives
    std::vector<double> 
        u_lon_tor_tmp( num_pts, 0. ),
        u_lat_tor_tmp( num_pts, 0. ),
        u_lon_pot_tmp( num_pts, 0. ),
        u_lat_pot_tmp( num_pts, 0. );

    //
    //// first ell derivatives
    //

    //toroidal_vel_from_F(  u_lon_tor_tmp, u_lat_tor_tmp, dl_coarse_Psi, source_data, mask);
    //potential_vel_from_F( u_lon_pot_tmp, u_lat_pot_tmp, dl_coarse_Phi, source_data, mask);
    perp_grad_scalar(  u_lon_tor_tmp, u_lat_tor_tmp, dl_coarse_Psi, source_data, mask);
    grad_scalar(       u_lon_pot_tmp, u_lat_pot_tmp, dl_coarse_Phi, source_data, mask);

    #pragma omp parallel \
    default( none ) \
    shared( mask, \
            u_spectrum_tor, v_spectrum_tor, spec_slope_tor, \
            u_spectrum_pot, v_spectrum_pot, spec_slope_pot, \
            u_spectrum_tot, v_spectrum_tot, spec_slope_tot, \
            u_spectrum_tor_HighPass, v_spectrum_tor_HighPass, spec_slope_tor_HighPass, \
            u_spectrum_pot_HighPass, v_spectrum_pot_HighPass, spec_slope_pot_HighPass, \
            u_spectrum_tot_HighPass, v_spectrum_tot_HighPass, spec_slope_tot_HighPass, \
            u_lon_tor_unfilt, u_lon_pot_unfilt, \
            u_lat_tor_unfilt, u_lat_pot_unfilt, \
            u_lon_tor, u_lat_tor, u_lon_tor_tmp, u_lat_tor_tmp,\
            u_lon_pot, u_lat_pot, u_lon_pot_tmp, u_lat_pot_tmp,  \
            u_lon_tot, u_lat_tot ) \
    private( index ) \
    firstprivate( num_pts )
    {
    #pragma omp for collapse(1) schedule(guided)
        for (index = 0; index < num_pts; ++index) {
            if ( mask.at(index) ) {

                // Merge toroidal and potential into total
                double u_lon_tot_tmp = u_lon_tor_tmp.at(index) + u_lon_pot_tmp.at(index);
                double u_lat_tot_tmp = u_lat_tor_tmp.at(index) + u_lat_pot_tmp.at(index);

                // Toroidal spectrum
                u_spectrum_tor.at(index) = 2 * u_lon_tor.at(index) * u_lon_tor_tmp.at(index);
                v_spectrum_tor.at(index) = 2 * u_lat_tor.at(index) * u_lat_tor_tmp.at(index);

                // Potential spectrum
                u_spectrum_pot.at(index) = 2 * u_lon_pot.at(index) * u_lon_pot_tmp.at(index);
                v_spectrum_pot.at(index) = 2 * u_lat_pot.at(index) * u_lat_pot_tmp.at(index);

                // Total spectrum
                u_spectrum_tot.at(index) = 2 * u_lon_tot.at(index) * u_lon_tot_tmp;
                v_spectrum_tot.at(index) = 2 * u_lat_tot.at(index) * u_lat_tot_tmp;

                // Contributions of first derivatives to slopes
                spec_slope_tor.at(index) =    my_pow<2>( u_lon_tor_tmp.at(index) ) 
                                            + my_pow<2>( u_lat_tor_tmp.at(index) );
                spec_slope_pot.at(index) =    my_pow<2>( u_lon_pot_tmp.at(index) ) 
                                            + my_pow<2>( u_lat_pot_tmp.at(index) );
                spec_slope_tot.at(index) =    my_pow<2>( u_lon_tot_tmp ) 
                                            + my_pow<2>( u_lat_tot_tmp );


                // 
                //// Also store the high-pass versions
                //

                // Toroidal spectrum
                u_spectrum_tor_HighPass.at(index) =  2 * u_lon_tor_unfilt.at(index) * u_lon_tor_tmp.at(index) - u_spectrum_tor.at(index);
                v_spectrum_tor_HighPass.at(index) =  2 * u_lat_tor_unfilt.at(index) * u_lat_tor_tmp.at(index) - v_spectrum_tor.at(index);

                // Potential spectrum
                u_spectrum_pot_HighPass.at(index) =  2 * u_lon_pot_unfilt.at(index) * u_lon_pot_tmp.at(index) - u_spectrum_pot.at(index);
                v_spectrum_pot_HighPass.at(index) =  2 * u_lat_pot_unfilt.at(index) * u_lat_pot_tmp.at(index) - v_spectrum_pot.at(index);

                // Total spectrum
                double u_lon_tot_unfilt = u_lon_tor_unfilt.at(index) + u_lon_pot_unfilt.at(index);
                double u_lat_tot_unfilt = u_lat_tor_unfilt.at(index) + u_lat_pot_unfilt.at(index);
                u_spectrum_tot_HighPass.at(index) =  2 * u_lon_tot_unfilt * u_lon_tot_tmp - u_spectrum_tot.at(index);
                v_spectrum_tot_HighPass.at(index) =  2 * u_lat_tot_unfilt * u_lat_tot_tmp - v_spectrum_tot.at(index);

                // Contributions of first derivatives to slopes
                spec_slope_tor_HighPass.at(index) =  -(  my_pow<2>( u_lon_tor_tmp.at(index) ) 
                                                       + my_pow<2>( u_lat_tor_tmp.at(index) ) );
                spec_slope_pot_HighPass.at(index) =  -(  my_pow<2>( u_lon_pot_tmp.at(index) ) 
                                                       + my_pow<2>( u_lat_pot_tmp.at(index) ) );
                spec_slope_tot_HighPass.at(index) =  -(  my_pow<2>( u_lon_tot_tmp ) 
                                                       + my_pow<2>( u_lat_tot_tmp ) );
            }
        }
    }

    //
    //// second ell derivatives
    //

    //toroidal_vel_from_F(  u_lon_tor_tmp, u_lat_tor_tmp, dll_coarse_Psi, source_data, mask);
    //potential_vel_from_F( u_lon_pot_tmp, u_lat_pot_tmp, dll_coarse_Phi, source_data, mask);
    perp_grad_scalar(  u_lon_tor_tmp, u_lat_tor_tmp, dll_coarse_Psi, source_data, mask);
    grad_scalar(       u_lon_pot_tmp, u_lat_pot_tmp, dll_coarse_Phi, source_data, mask);

    double numer, denom;

    #pragma omp parallel \
    default( none ) \
    shared( mask,  \
            u_spectrum_tor, v_spectrum_tor, spec_slope_tor, \
            u_spectrum_pot, v_spectrum_pot, spec_slope_pot, \
            u_spectrum_tot, v_spectrum_tot, spec_slope_tot, \
            u_spectrum_tor_HighPass, v_spectrum_tor_HighPass, spec_slope_tor_HighPass, \
            u_spectrum_pot_HighPass, v_spectrum_pot_HighPass, spec_slope_pot_HighPass, \
            u_spectrum_tot_HighPass, v_spectrum_tot_HighPass, spec_slope_tot_HighPass, \
            u_lon_tor, u_lat_tor, u_lon_tor_tmp, u_lat_tor_tmp,\
            u_lon_pot, u_lat_pot, u_lon_pot_tmp, u_lat_pot_tmp,  \
            u_lon_tot, u_lat_tot, \
            u_lon_tor_unfilt, u_lon_pot_unfilt, \
            u_lat_tor_unfilt, u_lat_pot_unfilt ) \
    private( index, numer, denom ) \
    firstprivate( num_pts, filter_scale )
    {
    #pragma omp for collapse(1) schedule(guided)
        for (index = 0; index < num_pts; ++index) {
            if ( mask.at(index) ) {

                double speed;

                // Merge toroidal and potential into total
                double u_lon_tot_tmp = u_lon_tor_tmp.at(index) + u_lon_pot_tmp.at(index);
                double u_lat_tot_tmp = u_lat_tor_tmp.at(index) + u_lat_pot_tmp.at(index);

                // Contributions of first derivatives to slopes
                spec_slope_tor.at(index) +=   u_lon_tor.at(index) * u_lon_tor_tmp.at(index)
                                            + u_lat_tor.at(index) * u_lat_tor_tmp.at(index);
                spec_slope_pot.at(index) +=   u_lon_pot.at(index) * u_lon_pot_tmp.at(index)
                                            + u_lat_pot.at(index) * u_lat_pot_tmp.at(index);
                spec_slope_tot.at(index) +=   u_lon_tot.at(index) * u_lon_tot_tmp
                                            + u_lat_tot.at(index) * u_lat_tot_tmp;


                // And now apply the scaling factors
                // If there's too little energy, then just set the slope to zero
                numer = filter_scale * spec_slope_tor.at(index);
                denom = u_spectrum_tor.at(index) + v_spectrum_tor.at(index);
                speed = sqrt( my_pow<2>(u_lon_tor[index]) + my_pow<2>(u_lat_tor[index]) );
                spec_slope_tor.at(index) = ( (fabs(denom) < 1e-15 ) or (speed < 1e-10) ) ? 0. : - 2 * (1 + (numer / denom));

                numer = filter_scale * spec_slope_pot.at(index);
                denom = u_spectrum_pot.at(index) + v_spectrum_pot.at(index);
                speed = sqrt( my_pow<2>(u_lon_pot[index]) + my_pow<2>(u_lat_pot[index]) );
                spec_slope_pot.at(index) = ( (fabs(denom) < 1e-15 ) or (speed < 1e-10) ) ? 0. : - 2 * (1 + (numer / denom));

                numer = filter_scale * spec_slope_tot.at(index);
                denom = u_spectrum_tot.at(index) + v_spectrum_tot.at(index);
                speed = sqrt( my_pow<2>(u_lon_tor[index] + u_lon_pot[index]) + my_pow<2>(u_lat_tor[index] + u_lat_pot[index]) );
                spec_slope_tot.at(index) = ( (fabs(denom) < 1e-15 ) or (speed < 1e-10) ) ? 0. : - 2 * (1 + (numer / denom));

                //
                //// Repeat for high-pass
                //

                // Merge toroidal and potential into total
                double u_lon_tot_unfilt = u_lon_tor_unfilt.at(index) + u_lon_pot_unfilt.at(index);
                double u_lat_tot_unfilt = u_lat_tor_unfilt.at(index) + u_lat_pot_unfilt.at(index);

                double u_lon_tot_HighPass = u_lon_tot_unfilt - ( u_lon_tor[index] + u_lon_pot[index] );
                double u_lat_tot_HighPass = u_lat_tot_unfilt - ( u_lat_tor[index] + u_lat_pot[index] );

                double u_lon_tor_HighPass = u_lon_tor_unfilt[index] - u_lon_tor[index];
                double u_lat_tor_HighPass = u_lat_tor_unfilt[index] - u_lat_tor[index];

                double u_lon_pot_HighPass = u_lon_pot_unfilt[index] - u_lon_pot[index];
                double u_lat_pot_HighPass = u_lat_pot_unfilt[index] - u_lat_pot[index];

                // Contributions of second derivatives to slopes
                spec_slope_tor_HighPass.at(index) +=   u_lon_tor_HighPass * u_lon_tor_tmp[index]
                                                     + u_lat_tor_HighPass * u_lat_tor_tmp[index];
                spec_slope_pot_HighPass.at(index) +=   u_lon_pot_HighPass * u_lon_pot_tmp[index]
                                                     + u_lat_pot_HighPass * u_lat_pot_tmp[index];
                spec_slope_tot_HighPass.at(index) +=   u_lon_tot_HighPass * u_lon_tot_tmp
                                                     + u_lat_tot_HighPass * u_lat_tot_tmp;


                // And now apply the scaling factors
                // If there's too little energy, then just set the slope to zero
                numer = filter_scale * spec_slope_tor_HighPass.at(index);
                denom = u_spectrum_tor_HighPass.at(index) + v_spectrum_tor_HighPass.at(index);
                speed = sqrt( my_pow<2>(u_lon_tor[index]) + my_pow<2>(u_lat_tor[index]) );
                spec_slope_tor_HighPass.at(index) = ( (fabs(denom) < 1e-15 ) or (speed < 1e-10) ) ? 0. : - 2 * (1 + (numer / denom));

                numer = filter_scale * spec_slope_pot_HighPass.at(index);
                denom = u_spectrum_pot_HighPass.at(index) + v_spectrum_pot_HighPass.at(index);
                speed = sqrt( my_pow<2>(u_lon_pot[index]) + my_pow<2>(u_lat_pot[index]) );
                spec_slope_pot_HighPass.at(index) = ( (fabs(denom) < 1e-15 ) or (speed < 1e-10) ) ? 0. : - 2 * (1 + (numer / denom));

                numer = filter_scale * spec_slope_tot_HighPass.at(index);
                denom = u_spectrum_tot_HighPass.at(index) + v_spectrum_tot_HighPass.at(index);
                speed = sqrt( my_pow<2>(u_lon_tor[index] + u_lon_pot[index]) + my_pow<2>(u_lat_tor[index] + u_lat_pot[index]) );
                spec_slope_tot_HighPass.at(index) = ( (fabs(denom) < 1e-15 ) or (speed < 1e-10) ) ? 0. : - 2 * (1 + (numer / denom));


                //
                ////  And apply the scaling to the power spectrum
                //

                const double half_rho_ell_sqr = 0.5 * constants::rho0 * my_pow<2>(filter_scale);

                u_spectrum_tor.at(index) = - half_rho_ell_sqr * u_spectrum_tor.at(index);
                v_spectrum_tor.at(index) = - half_rho_ell_sqr * v_spectrum_tor.at(index);

                u_spectrum_pot.at(index) = - half_rho_ell_sqr * u_spectrum_pot.at(index);
                v_spectrum_pot.at(index) = - half_rho_ell_sqr * v_spectrum_pot.at(index);

                u_spectrum_tot.at(index) = - half_rho_ell_sqr * u_spectrum_tot.at(index);
                v_spectrum_tot.at(index) = - half_rho_ell_sqr * v_spectrum_tot.at(index);

                // Repeat for high-pass
                u_spectrum_tor_HighPass.at(index) = - half_rho_ell_sqr * u_spectrum_tor_HighPass.at(index);
                v_spectrum_tor_HighPass.at(index) = - half_rho_ell_sqr * v_spectrum_tor_HighPass.at(index);

                u_spectrum_pot_HighPass.at(index) = - half_rho_ell_sqr * u_spectrum_pot_HighPass.at(index);
                v_spectrum_pot_HighPass.at(index) = - half_rho_ell_sqr * v_spectrum_pot_HighPass.at(index);

                u_spectrum_tot_HighPass.at(index) = - half_rho_ell_sqr * u_spectrum_tot_HighPass.at(index);
                v_spectrum_tot_HighPass.at(index) = - half_rho_ell_sqr * v_spectrum_tot_HighPass.at(index);
            }
        }
    }

}
