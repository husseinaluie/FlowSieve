#include <vector>
#include <omp.h>
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
 * @param[in]       ux,uy,uz                        coarse Cartesian velocity components
 * @param[in]       uxux,uxuy,uxuz,uyuy,uyuz,uzuz   coarse velocity products (e.g. bar(u*v) )  
 * @param[in]       comm                            MPI communicator object
 *
 */
void compute_Pi(
        std::vector<double> & energy_transfer,
        const dataset & source_data,
        const std::vector<double> & ux,
        const std::vector<double> & uy,
        const std::vector<double> & uz,
        const std::vector<double> & uxux,
        const std::vector<double> & uxuy,
        const std::vector<double> & uxuz,
        const std::vector<double> & uyuy,
        const std::vector<double> & uyuz,
        const std::vector<double> & uzuz,
        const std::vector<double> * ux_in_tau,
        const std::vector<double> * uy_in_tau,
        const std::vector<double> * uz_in_tau,
        const MPI_Comm comm
        ) {

    const std::vector<bool> &mask = source_data.mask;

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

    double ui_j, uj_i;
    std::vector<double> tau_ij( Npts );

    // Some convenience handles
    //   note: the pointers aren't constant, but the things
    //         to which they are pointing are
    double ui_loc, uj_loc, uiuj_loc;
    const std::vector<double> *uiuj, *ui, *uj, *ui_in_tau, *uj_in_tau;

    //const std::vector<double> *ux_in_tau = ( ux_for_tau == NULL ) ? &ux[0] : ux_for_tau;
    //const std::vector<double> *uy_in_tau = ( uy_for_tau == NULL ) ? &uy[0] : uy_for_tau;
    //const std::vector<double> *uz_in_tau = ( uz_for_tau == NULL ) ? &uz[0] : uz_for_tau;
    if ( ux_in_tau == NULL ) { ux_in_tau = &ux; }
    if ( uy_in_tau == NULL ) { uy_in_tau = &uy; }
    if ( uz_in_tau == NULL ) { uz_in_tau = &uz; }

    // Set up the derivatives to pass through the differentiation functions
    std::vector<double*> x_deriv_vals, y_deriv_vals, z_deriv_vals;
    std::vector<const std::vector<double>*> deriv_fields(2);

    // Zero out energy transfer before we start
    std::fill( energy_transfer.begin(), energy_transfer.end(), 0.);

    for (ii = 0; ii < 3; ii++) {
        for (jj = 0; jj < 3; jj++) {

            //   Assign some handy pointers: uiuj, ui, uj
            //
            //   0 -> x
            //   1 -> y
            //   2 -> z

            // ui
            ui = ( ii == 0 ) ? &ux : ( ii == 1 ) ? &uy : &uz;
            uj = ( jj == 0 ) ? &ux : ( jj == 1 ) ? &uy : &uz;

            ui_in_tau = ( ii == 0 ) ? ux_in_tau : ( ii == 1 ) ? uy_in_tau : uz_in_tau;
            uj_in_tau = ( jj == 0 ) ? ux_in_tau : ( jj == 1 ) ? uy_in_tau : uz_in_tau;

            deriv_fields[0] = ui;
            deriv_fields[1] = uj;

            // uiuj (note that they're symmetric i.e. uiuj = ujui)
            uiuj =  ( ii == 0 ) ? ( ( jj == 0 ) ? &uxux : ( jj == 1 ) ? &uxuy : &uxuz ) :
                    ( ii == 1 ) ? ( ( jj == 0 ) ? &uxuy : ( jj == 1 ) ? &uyuy : &uyuz ) :
                                  ( ( jj == 0 ) ? &uxuz : ( jj == 1 ) ? &uyuz : &uzuz ) ;

            // Compute tau_ij 
            #pragma omp parallel \
            default(none) \
            shared(tau_ij, mask, ui_in_tau, uj_in_tau, uiuj, source_data) \
            private(index, uiuj_loc, ui_loc, uj_loc) \
            firstprivate( Npts )
            {
                #pragma omp for collapse(1) schedule(guided)
                for (index = 0; index < Npts; index++) {

                    if ( constants::FILTER_OVER_LAND or mask.at(index) ) {
                        ui_loc   = ui_in_tau->at(  index);
                        uj_loc   = uj_in_tau->at(  index);
                        uiuj_loc = uiuj->at(index);

                        tau_ij.at(index) = uiuj_loc - ui_loc * uj_loc;
                    }
                }
            }

            #pragma omp parallel default(none) \
            shared( source_data, energy_transfer, mask, \
                    ii, jj, ui, uj, tau_ij, deriv_fields)\
            private( Itime, Idepth, Ilat, Ilon, index, pi_tmp, ui_j, uj_i, \
                     x_deriv_vals, y_deriv_vals, z_deriv_vals) \
            firstprivate( Npts )
            {

                x_deriv_vals.resize(2);
                y_deriv_vals.resize(2);
                z_deriv_vals.resize(2);

                // Now set the appropriate derivative pointers
                //   in order to compute
                //   ui,j and uj,i
                x_deriv_vals.at(0) = (jj == 0) ? &ui_j : NULL;
                x_deriv_vals.at(1) = (ii == 0) ? &uj_i : NULL;

                y_deriv_vals.at(0) = (jj == 1) ? &ui_j : NULL;
                y_deriv_vals.at(1) = (ii == 1) ? &uj_i : NULL;

                z_deriv_vals.at(0) = (jj == 2) ? &ui_j : NULL;
                z_deriv_vals.at(1) = (ii == 2) ? &uj_i : NULL;

                // Now actually compute Pi
                //   in particular, compute S_ij * tau_ij
                #pragma omp for collapse(1) schedule(guided)
                for (index = 0; index < Npts; index++) {

                    if ( constants::FILTER_OVER_LAND or mask.at(index) ) {

                        source_data.index1to4_local( index, Itime, Idepth, Ilat, Ilon);

                        // Compute the desired derivatives
                        Cart_derivatives_at_point(
                                x_deriv_vals, y_deriv_vals, z_deriv_vals, deriv_fields, 
                                source_data, Itime, Idepth, Ilat, Ilon,
                                1, constants::DiffOrd);

                        double Sij = 0.5 * ( ui_j + uj_i );
                        pi_tmp = - constants::rho0 * Sij * tau_ij.at(index);
                        energy_transfer.at(index) += pi_tmp;
                    }
                }
            }
        }
    }
    #if DEBUG >= 2
    if (wRank == 0) { fprintf(stdout, "     ... done.\n"); }
    #endif
}
