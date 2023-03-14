#include <vector>
#include <omp.h>
#include <math.h>    
#include "../functions.hpp"
#include "../constants.hpp"
#include "../differentiation_tools.hpp"

/*!
 * \brief Compute the energy transfer through the current filter scale
 *
 * In particular, computes \f$ \rho_0 * ( u_i \tau_{ij,j} - (u_i \tau_{ij})_{,j}  ) \f$
 * 
 * This computation is applied to the Cartesian velocity components
 *
 * @param[in,out]   energy_transfer                 where to store the computed values (array)
 * @param[in]       source_data                     dataset class instance containing data (Psi, Phi, etc)
 * @param[in]       ulon,ulat                       coarse Spherical velocity components
 * @param[in]       ulon_ulon,ulon_ulat,ulat_ulat   coarse velocity products (e.g. bar(u*v) )  
 * @param[in]       comm                            MPI communicator object
 *
 */
void compute_Pi_Helmholtz(
        std::vector<double> & energy_transfer,
        const dataset & source_data,
        const std::vector<double> & ulon,
        const std::vector<double> & ulat,
        const std::vector<double> & ulon_ulon,
        const std::vector<double> & ulon_ulat,
        const std::vector<double> & ulat_ulat,
        const MPI_Comm comm
        ) {

    const std::vector<double>   &latitude   = source_data.latitude,
                                &longitude  = source_data.longitude;

    const std::vector<bool> &mask = source_data.mask;

    const int   Ntime   = source_data.Ntime,
                Ndepth  = source_data.Ndepth,
                Nlat    = source_data.Nlat,
                Nlon    = source_data.Nlon;

    #if DEBUG >= 2
    int wRank, wSize;
    MPI_Comm_rank( comm, &wRank );
    MPI_Comm_size( comm, &wSize );

    if (wRank == 0) { fprintf(stdout, "  Starting Pi computation.\n"); }
    #endif

    double pi_tmp;
    int Itime, Idepth, Ilat, Ilon, ii, jj;
    size_t index;
    const size_t Npts = energy_transfer.size();

    double tau_ij_j, u_i_tau_ij_j, u_i_j, u_j_i;
    std::vector<double> tau_ij;
    std::vector<double> u_i_tau_ij;
    tau_ij.resize(    ulon.size() );
    u_i_tau_ij.resize(ulon.size() );

    // Some convenience handles
    //   note: the pointers aren't constant, but the things
    //         to which they are pointing are
    double ui_loc, uj_loc, uiuj_loc;
    const std::vector<double> *uiuj, *ui, *uj;

    // Set up the derivatives to pass through the differentiation functions
    std::vector<double*> i_deriv_vals, j_deriv_vals;
    std::vector<const std::vector<double>*> i_deriv_fields(1), j_deriv_fields(3);

    j_deriv_fields.at(1) = &tau_ij;
    j_deriv_fields.at(2) = &u_i_tau_ij;

    // Zero out energy transfer before we start
    std::fill( energy_transfer.begin(), energy_transfer.end(), 0.);

    for (ii = 0; ii < 2; ii++) {
        for (jj = 0; jj < 2; jj++) {

            //   Assign the handy pointers: uiuj, ui, uj
            //
            //   0 -> lon
            //   1 -> lat

            // ui
            switch (ii) {
                case 0 : ui = &ulon; j_deriv_fields.at(0) = &ulon; break;
                case 1 : ui = &ulat; j_deriv_fields.at(0) = &ulat; break;
            }

            // uj
            switch (jj) {
                case 0 : uj = &ulon; i_deriv_fields.at(0) = &ulon; break;
                case 1 : uj = &ulat; i_deriv_fields.at(0) = &ulat; break;
            }

            // uiuj (note that they're symmetric i.e. uiuj = ujui)
            switch (ii) {
                case 0 :
                    switch (jj) {
                        case 0 : uiuj = &ulon_ulon; break;
                        case 1 : uiuj = &ulon_ulat; break;
                    }
                    break;
                case 1 :
                    switch (jj) {
                        case 0 : uiuj = &ulon_ulat; break;
                        case 1 : uiuj = &ulat_ulat; break;
                    }
                    break;
            }

            // First, compute the appropriate
            //   tau_ij and u_i * tau_ij
            #pragma omp parallel \
            default(none) \
            shared(tau_ij, u_i_tau_ij, mask, ui, uj, uiuj)\
            private(index, uiuj_loc, ui_loc, uj_loc) \
            firstprivate( Npts )
            {
                #pragma omp for collapse(1) schedule(guided)
                for (index = 0; index < Npts; index++) {

                    if ( mask.at(index) ) {

                        ui_loc   = ui->at(  index);
                        uj_loc   = uj->at(  index);
                        uiuj_loc = uiuj->at(index);

                        tau_ij.at(index) = uiuj_loc - ui_loc * uj_loc;
                        u_i_tau_ij.at(index) = ui_loc * tau_ij.at(index);

                    }
                }
            }

            #pragma omp parallel \
            default(none) \
            shared(energy_transfer, latitude, longitude, mask, source_data, \
                    ii, jj, ui, uj, tau_ij, u_i_tau_ij, i_deriv_fields, j_deriv_fields)\
            private(Itime, Idepth, Ilat, Ilon, index,\
                    pi_tmp, tau_ij_j, u_i_tau_ij_j, u_i_j, u_j_i, i_deriv_vals, j_deriv_vals) \
            firstprivate( Npts, Nlon, Nlat, Ndepth, Ntime )
            {

                // Now set the appropriate derivative pointers
                //   in order to compute
                //     tau_ij,j
                //     (u_i * tau_ij)_,j
                //     u_i,j
                j_deriv_vals.resize(3);
                j_deriv_vals.at(0) = &u_i_j;
                j_deriv_vals.at(1) = &tau_ij_j;
                j_deriv_vals.at(2) = &u_i_tau_ij_j;

                i_deriv_vals.resize(1);
                i_deriv_vals.at(0) = &u_j_i;

                // Now actually compute Pi
                //   in particular, compute
                //           u_i * tau_ij,j - (u_i * tau_ij)_,j
                #pragma omp for collapse(1) schedule(guided)
                for (index = 0; index < Npts; index++) {

                    if ( mask.at(index) ) {

                        Index1to4(index, Itime, Idepth, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon);

                        switch (ii) {
                            case 0 : spher_derivative_at_point( i_deriv_vals, i_deriv_fields, longitude, "lon", source_data,
                                Itime, Idepth, Ilat, Ilon, mask); break;
                            case 1 : spher_derivative_at_point( i_deriv_vals, i_deriv_fields, latitude, "lat", source_data,
                                Itime, Idepth, Ilat, Ilon, mask); break;
                        }
                        switch (jj) {
                            case 0 : spher_derivative_at_point( j_deriv_vals, j_deriv_fields, longitude, "lon", source_data,
                                Itime, Idepth, Ilat, Ilon, mask); break;
                            case 1 : spher_derivative_at_point( j_deriv_vals, j_deriv_fields, latitude, "lat", source_data,
                                Itime, Idepth, Ilat, Ilon, mask); break;
                        }

                        double i_scale_factor = ( (ii == 0) ? 1. / ( cos(latitude.at(Ilat) ) ) : 1. ) / constants::R_earth;
                        u_j_i        *= i_scale_factor;

                        double j_scale_factor = ( (jj == 0) ? 1. / ( cos(latitude.at(Ilat) ) ) : 1. ) / constants::R_earth;
                        tau_ij_j     *= j_scale_factor;
                        u_i_tau_ij_j *= j_scale_factor;
                        u_i_j        *= j_scale_factor;

                        // u_i * tau_ij,j - (u_i * tau_ij)_,j
                        //pi_tmp = ui->at(index) * tau_ij_j  -  u_i_tau_ij_j;
                        // - 0.5 * S_ij * tau_ij
                        pi_tmp = - 0.5 * ( u_i_j + u_j_i ) * tau_ij.at(index);
                        energy_transfer.at(index) += constants::rho0 * pi_tmp;
                    }
                }
            }
        }
    }
    #if DEBUG >= 2
    if (wRank == 0) { fprintf(stdout, "     ... done.\n"); }
    #endif
}

