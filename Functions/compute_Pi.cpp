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
        const MPI_Comm comm
        ) {

    const std::vector<double>   &latitude   = source_data.latitude,
                                &longitude  = source_data.longitude;

    const std::vector<bool> &mask = source_data.mask;

    const int   Ntime   = source_data.Ntime,
                Ndepth  = source_data.Ndepth,
                Nlat    = source_data.Nlat,
                Nlon    = source_data.Nlon;

    const int OMP_chunksize = get_omp_chunksize(Nlat,Nlon);

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
    std::vector<double> tau_ij( ux.size() );

    // Some convenience handles
    //   note: the pointers aren't constant, but the things
    //         to which they are pointing are
    double ui_loc, uj_loc, uiuj_loc;
    const std::vector<double> *uiuj, *ui, *uj;

    // Set up the derivatives to pass through the differentiation functions
    std::vector<double*> x_deriv_vals, y_deriv_vals, z_deriv_vals;
    std::vector<const std::vector<double>*> deriv_fields(2);

    // Zero out energy transfer before we start
    std::fill( energy_transfer.begin(), energy_transfer.end(), 0.);

    for (ii = 0; ii < 3; ii++) {
        for (jj = 0; jj < 3; jj++) {

            //   Assign the handy pointers: uiuj, ui, uj
            //
            //   0 -> x
            //   1 -> y
            //   2 -> z

            // ui
            switch (ii) {
                case 0 : ui = &ux; break;
                case 1 : ui = &uy; break;
                case 2 : ui = &uz; break;
            }

            // uj
            switch (jj) {
                case 0 : uj = &ux; break;
                case 1 : uj = &uy; break;
                case 2 : uj = &uz; break;
            }

            deriv_fields[0] = ui;
            deriv_fields[1] = uj;

            // uiuj (note that they're symmetric i.e. uiuj = ujui)
            switch (ii) {
                case 0 :
                    switch (jj) {
                        case 0 : uiuj = &uxux; break;
                        case 1 : uiuj = &uxuy; break;
                        case 2 : uiuj = &uxuz; break;
                    }
                    break;
                case 1 :
                    switch (jj) {
                        case 0 : uiuj = &uxuy; break;
                        case 1 : uiuj = &uyuy; break;
                        case 2 : uiuj = &uyuz; break;
                    }
                    break;
                case 2 :
                    switch (jj) {
                        case 0 : uiuj = &uxuz; break;
                        case 1 : uiuj = &uyuz; break;
                        case 2 : uiuj = &uzuz; break;
                    }
                    break;
            }

            // Compute tau_ij 
            #pragma omp parallel default(none) shared(tau_ij, mask, ui, uj, uiuj) private(index, uiuj_loc, ui_loc, uj_loc)
            {
                #pragma omp for collapse(1) schedule(dynamic, OMP_chunksize)
                for (index = 0; index < Npts; index++) {

                    if ( mask.at(index) ) {

                        ui_loc   = ui->at(  index);
                        uj_loc   = uj->at(  index);
                        uiuj_loc = uiuj->at(index);

                        tau_ij.at(index) = uiuj_loc - ui_loc * uj_loc;
                    }
                }
            }

            #pragma omp parallel default(none) \
            shared(energy_transfer, latitude, longitude, mask, ii, jj, ui, uj, tau_ij, deriv_fields)\
            private(Itime, Idepth, Ilat, Ilon, index, pi_tmp, ui_j, uj_i, x_deriv_vals, y_deriv_vals, z_deriv_vals)
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
                #pragma omp for collapse(1) schedule(dynamic, OMP_chunksize)
                for (index = 0; index < Npts; index++) {

                    if ( mask.at(index) ) {

                        Index1to4(index, Itime, Idepth, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon);

                        // Compute the desired derivatives
                        Cart_derivatives_at_point(
                                x_deriv_vals, y_deriv_vals, z_deriv_vals, deriv_fields,
                                latitude, longitude, Itime, Idepth, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon,
                                mask);

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
