#include <math.h>
#include <algorithm>
#include <vector>
#include <omp.h>
#include <mpi.h>

#include "../functions.hpp"
#include "../constants.hpp"
#include "../preprocess.hpp"

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
        const std::vector<double> & u_lon_tot, 
        const std::vector<double> & u_lon_tor, 
        const std::vector<double> & u_lon_pot,
        const std::vector<double> & u_lat_tot, 
        const std::vector<double> & u_lat_tor, 
        const std::vector<double> & u_lat_pot,
        const std::vector<double> & dl_coarse_Phi, 
        const std::vector<double> & dl_coarse_Psi,
        const std::vector<double> & dll_coarse_Phi, 
        const std::vector<double> & dll_coarse_Psi,
        const dataset & source_data,
        const double filter_scale
        ) {

    // Get dimension sizes
    const int   Ntime   = source_data.Ntime,    // this is the MPI-local Ntime, not the full Ntime
                Ndepth  = source_data.Ndepth,   // this is the MPI-local Ndepth, not the full Ndepth
                Nlat    = source_data.Nlat,
                Nlon    = source_data.Nlon;

    // Create some tidy names for variables
    const std::vector<double>   &latitude   = source_data.latitude,
                                &longitude  = source_data.longitude;
    const std::vector<bool> &mask = source_data.mask;

    const size_t num_pts = u_lon_tot.size();
    size_t index;

    // Create some (temporary) storage arrays to hold the ell-derivatives
    std::vector<double> u_lon_tor_tmp( num_pts, 0. ),
        u_lat_tor_tmp( num_pts, 0. ),
        u_lon_pot_tmp( num_pts, 0. ),
        u_lat_pot_tmp( num_pts, 0. ),
        u_lon_tot_tmp( num_pts, 0. ),
        u_lat_tot_tmp( num_pts, 0. );

    //
    //// first ell derivatives
    //

    toroidal_vel_from_F( u_lon_tor_tmp, u_lat_tor_tmp, dl_coarse_Psi, 
            longitude, latitude, Ntime, Ndepth, Nlat, Nlon, mask);
    potential_vel_from_F( u_lon_pot_tmp, u_lat_pot_tmp, dl_coarse_Phi, 
            longitude, latitude, Ntime, Ndepth, Nlat, Nlon, mask);

    #pragma omp parallel \
    default( none ) \
    shared( mask, \
            u_spectrum_tor, v_spectrum_tor, spec_slope_tor, \
            u_spectrum_pot, v_spectrum_pot, spec_slope_pot, \
            u_spectrum_tot, v_spectrum_tot, spec_slope_tot, \
            u_lon_tor, u_lat_tor, u_lon_tor_tmp, u_lat_tor_tmp,\
            u_lon_pot, u_lat_pot, u_lon_pot_tmp, u_lat_pot_tmp,  \
            u_lon_tot, u_lat_tot, u_lon_tot_tmp, u_lat_tot_tmp ) \
    private( index ) \
    firstprivate( num_pts )
    {
    #pragma omp for collapse(1) schedule(guided)
        for (index = 0; index < num_pts; ++index) {
            if ( mask.at(index) ) {

                // Merge toroidal and potential into total
                u_lon_tot_tmp.at(index) = u_lon_tor_tmp.at(index) + u_lon_pot_tmp.at(index);
                u_lat_tot_tmp.at(index) = u_lat_tor_tmp.at(index) + u_lat_pot_tmp.at(index);

                // Totoidal spectrum
                u_spectrum_tor.at(index) = 2 * u_lon_tor.at(index) * u_lon_tor_tmp.at(index);
                v_spectrum_tor.at(index) = 2 * u_lat_tor.at(index) * u_lat_tor_tmp.at(index);

                // Potential spectrum
                u_spectrum_pot.at(index) = 2 * u_lon_pot.at(index) * u_lon_pot_tmp.at(index);
                v_spectrum_pot.at(index) = 2 * u_lat_pot.at(index) * u_lat_pot_tmp.at(index);

                // Total spectrum
                u_spectrum_tot.at(index) = 2 * u_lon_tot.at(index) * u_lon_tot_tmp.at(index);
                v_spectrum_tot.at(index) = 2 * u_lat_tot.at(index) * u_lat_tot_tmp.at(index);

                // Contributions of first derivatives to slopes
                spec_slope_tor.at(index) =    pow( u_lon_tor_tmp.at(index), 2 ) 
                                            + pow( u_lat_tor_tmp.at(index), 2 );
                spec_slope_pot.at(index) =    pow( u_lon_pot_tmp.at(index), 2 ) 
                                            + pow( u_lat_pot_tmp.at(index), 2 );
                spec_slope_tot.at(index) =    pow( u_lon_tot_tmp.at(index), 2 ) 
                                            + pow( u_lat_tot_tmp.at(index), 2 );
            }
        }
    }

    //
    //// second ell derivatives
    //

    toroidal_vel_from_F( u_lon_tor_tmp, u_lat_tor_tmp, dll_coarse_Psi, 
            longitude, latitude, Ntime, Ndepth, Nlat, Nlon, mask);
    potential_vel_from_F( u_lon_pot_tmp, u_lat_pot_tmp, dll_coarse_Phi, 
            longitude, latitude, Ntime, Ndepth, Nlat, Nlon, mask);

    double numer, denom;

    #pragma omp parallel \
    default( none ) \
    shared( mask,  \
            u_spectrum_tor, v_spectrum_tor, spec_slope_tor, \
            u_spectrum_pot, v_spectrum_pot, spec_slope_pot, \
            u_spectrum_tot, v_spectrum_tot, spec_slope_tot, \
            u_lon_tor, u_lat_tor, u_lon_tor_tmp, u_lat_tor_tmp,\
            u_lon_pot, u_lat_pot, u_lon_pot_tmp, u_lat_pot_tmp,  \
            u_lon_tot, u_lat_tot, u_lon_tot_tmp, u_lat_tot_tmp ) \
    private( index, numer, denom ) \
    firstprivate( num_pts, filter_scale )
    {
    #pragma omp for collapse(1) schedule(guided)
        for (index = 0; index < num_pts; ++index) {
            if ( mask.at(index) ) {

                // Merge toroidal and potential into total
                u_lon_tot_tmp.at(index) = u_lon_tor_tmp.at(index) + u_lon_pot_tmp.at(index);
                u_lat_tot_tmp.at(index) = u_lat_tor_tmp.at(index) + u_lat_pot_tmp.at(index);

                // Contributions of first derivatives to slopes
                spec_slope_tor.at(index) +=   u_lon_tor.at(index) * u_lon_tor_tmp.at(index)
                                            + u_lat_tor.at(index) * u_lat_tor_tmp.at(index);
                spec_slope_pot.at(index) +=   u_lon_pot.at(index) * u_lon_pot_tmp.at(index)
                                            + u_lat_pot.at(index) * u_lat_pot_tmp.at(index);
                spec_slope_tot.at(index) +=   u_lon_tot.at(index) * u_lon_tot_tmp.at(index)
                                            + u_lat_tot.at(index) * u_lat_tot_tmp.at(index);


                // And now apply the scaling factors
                // If there's too little energy, then just set the slope to zero
                numer = filter_scale * spec_slope_tor.at(index);
                denom = u_spectrum_tor.at(index) + v_spectrum_tor.at(index);
                spec_slope_tor.at(index) = (fabs(denom) < 1e-20 ) ? 0. : - 2 - (numer / denom);

                numer = filter_scale * spec_slope_pot.at(index);
                denom = u_spectrum_pot.at(index) + v_spectrum_pot.at(index);
                spec_slope_pot.at(index) = (fabs(denom) < 1e-20 ) ? 0. : - 2 - (numer / denom);

                numer = filter_scale * spec_slope_tot.at(index);
                denom = u_spectrum_tot.at(index) + v_spectrum_tot.at(index);
                spec_slope_tot.at(index) = (fabs(denom) < 1e-20 ) ? 0. : - 2 - (numer / denom);


                // And apply the scaling to the power spectrum
                u_spectrum_tor.at(index) = - 0.5 * constants::rho0 * pow(filter_scale,2) * u_spectrum_tor.at(index);
                v_spectrum_tor.at(index) = - 0.5 * constants::rho0 * pow(filter_scale,2) * v_spectrum_tor.at(index);

                u_spectrum_pot.at(index) = - 0.5 * constants::rho0 * pow(filter_scale,2) * u_spectrum_pot.at(index);
                v_spectrum_pot.at(index) = - 0.5 * constants::rho0 * pow(filter_scale,2) * v_spectrum_pot.at(index);

                u_spectrum_tot.at(index) = - 0.5 * constants::rho0 * pow(filter_scale,2) * u_spectrum_tot.at(index);
                v_spectrum_tot.at(index) = - 0.5 * constants::rho0 * pow(filter_scale,2) * v_spectrum_tot.at(index);
            }
        }
    }

}
