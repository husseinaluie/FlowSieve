#include <vector>
#include <stdio.h>
#include <iostream>
#include <omp.h>
#include "../functions.hpp"
#include "../constants.hpp"
#include "../differentiation_tools.hpp"

/*!
 * \brief Compute the energy transfer through the current filter scale
 *
 * In particular, computes \$ \rho_0 * ( u_i \tau_{ij,j} - (u_i \tau_{ij})_{,j}  ) \$
 * 
 * This computation is applied to the Cartesian velocity components
 *
 * @param[in,out]   enstrophy_transfer              where to store the computed values (array)
 * @param[in]       ux,uy,uz                        coarse Cartesian velocity components
 * @param[in]       coarse_vort_r                   coarse radial vorticity
 * @param[in]       vort_ux,vort_uy,vort_uz         coarse vort-velocity products (e.g. bar(omega * u_x) )  
 * @param[in]       Ntime,Ndepth,Nlat,Nlon          Size of dimensions (MPI-local sizes)
 * @param[in]       longitude,latitude              1D dimension vectors
 * @param[in]       mask                            2D array to distinguish land from water
 * @param[in]       comm                            MPI communicator object
 *
 */
void compute_Z(
        std::vector<double> & enstrophy_transfer,
        const std::vector<double> & ux,
        const std::vector<double> & uy,
        const std::vector<double> & uz,
        const std::vector<double> & coarse_vort_r,
        const std::vector<double> & vort_ux,
        const std::vector<double> & vort_uy,
        const std::vector<double> & vort_uz,
        const int Ntime,
        const int Ndepth,
        const int Nlat,
        const int Nlon,
        const std::vector<double> & longitude,
        const std::vector<double> & latitude,
        const std::vector<bool> & mask,
        const MPI_Comm comm
        ) {

    const int OMP_chunksize = get_omp_chunksize(Nlat,Nlon);

    #if DEBUG >= 2
    int wRank, wSize;
    MPI_Comm_rank( comm, &wRank );
    MPI_Comm_size( comm, &wSize );

    if (wRank == 0) { fprintf(stdout, "  Starting Z computation.\n"); }
    #endif

    double Z_tmp;
    const int ii = 0;
    int Itime, Idepth, Ilat, Ilon, jj;
    size_t index;
    const size_t Npts = enstrophy_transfer.size();

    double tau_ij_j, u_i_tau_ij_j;
    std::vector<double> tau_ij;
    std::vector<double> u_i_tau_ij;
    tau_ij.resize(ux.size());
    u_i_tau_ij.resize(ux.size());

    // Some convenience handles
    //   note: the pointers aren't constant, but the things
    //         to which they are pointing are
    double ui_loc, uj_loc, uiuj_loc;
    const std::vector<double> *uiuj, *ui, *uj;

    // Set up the derivatives to pass through the differentiation functions
    std::vector<double*> x_deriv_vals, y_deriv_vals, z_deriv_vals;
    std::vector<const std::vector<double>*> deriv_fields;

    deriv_fields.push_back(&tau_ij);
    deriv_fields.push_back(&u_i_tau_ij);

    // Zero out enstrophy transfer before we start
    #pragma omp parallel \
    default(none) shared(mask, enstrophy_transfer) \
    private(index)
    {
        #pragma omp for collapse(1) schedule(static)
        for (index = 0; index < Npts; index++) {
            enstrophy_transfer.at(index) = mask.at(index) ? 0. : constants::fill_value;
        }
    }

    for (jj = 0; jj < 3; jj++) {

        //   Assign the handy pointers: uiuj, ui, uj
        //
        //   0 -> x
        //   1 -> y
        //   2 -> z

        // ui
        ui = &coarse_vort_r;

        // uj
        switch (jj) {
            case 0 : uj = &ux; break;
            case 1 : uj = &uy; break;
            case 2 : uj = &uz; break;
        }

        // uiuj (note that they're symmetric i.e. uiuj = ujui)
        switch (jj) {
            case 0 : uiuj = &vort_ux; break;
            case 1 : uiuj = &vort_uy; break;
            case 2 : uiuj = &vort_uz; break;
        }
        
            
        // First, compute the appropriate
        //   tau_ij and u_i * tau_ij
        #pragma omp parallel \
        default(none) \
        shared(tau_ij, u_i_tau_ij, mask, ui, uj, uiuj)\
        private(index, uiuj_loc, ui_loc, uj_loc)
        {
            #pragma omp for collapse(1) schedule(dynamic, OMP_chunksize)
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
        shared(enstrophy_transfer, latitude, longitude, mask, stdout, \
                jj, ui, tau_ij, u_i_tau_ij, deriv_fields,std::cout)\
        private(Itime, Idepth, Ilat, Ilon, index,\
                Z_tmp, tau_ij_j, u_i_tau_ij_j,\
                x_deriv_vals, y_deriv_vals, z_deriv_vals)
        {

            x_deriv_vals.resize(2);
            y_deriv_vals.resize(2);
            z_deriv_vals.resize(2);

            // Now set the appropriate derivative pointers in order to compute
            //     tau_ij,j
            //     (u_i * tau_ij)_,j
            x_deriv_vals.at(0) = (jj == 0) ? &tau_ij_j     : NULL;
            x_deriv_vals.at(1) = (jj == 0) ? &u_i_tau_ij_j : NULL;

            y_deriv_vals.at(0) = (jj == 1) ? &tau_ij_j     : NULL;
            y_deriv_vals.at(1) = (jj == 1) ? &u_i_tau_ij_j : NULL;

            z_deriv_vals.at(0) = (jj == 2) ? &tau_ij_j     : NULL;
            z_deriv_vals.at(1) = (jj == 2) ? &u_i_tau_ij_j : NULL;

            // Now actually compute Z -  in particular, compute
            //           u_i * tau_ij,j - (u_i * tau_ij)_,j               
            #pragma omp for collapse(1) schedule(dynamic, OMP_chunksize)
            for (index = 0; index < Npts; index++) {

                if ( mask.at(index) ) {

                    Index1to4(index, Itime, Idepth, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon);

                    // Compute the desired derivatives
                    Cart_derivatives_at_point(
                            x_deriv_vals, y_deriv_vals,
                            z_deriv_vals, deriv_fields,
                            latitude, longitude,
                            Itime, Idepth, Ilat, Ilon,
                            Ntime, Ndepth, Nlat, Nlon,
                            mask);

                    // u_i * tau_ij,j - (u_i * tau_ij)_,j
                    Z_tmp = ui->at(index) * tau_ij_j  -  u_i_tau_ij_j;
                    enstrophy_transfer.at(index) += Z_tmp;
                        
                }
            }
        }
    }
    #if DEBUG >= 2
    if (wRank == 0) { fprintf(stdout, "     ... done.\n"); }
    #endif
}
