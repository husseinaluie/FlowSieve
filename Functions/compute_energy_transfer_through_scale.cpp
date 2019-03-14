#include <vector>
#include <omp.h>
#include "../functions.hpp"
#include "../constants.hpp"

void compute_energy_transfer_through_scale(
        std::vector<double> & energy_transfer,  /**< [in] where to store the energy transfer */
        const std::vector<double> & ux,         /**< [in] filtered u_x */
        const std::vector<double> & uy,         /**< [in] filtered u_y */
        const std::vector<double> & uz,         /**< [in] filtered u_z */
        const std::vector<double> & uxux,       /**< [in] filtered u_x * u_x */
        const std::vector<double> & uxuy,       /**< [in] filtered u_x * u_y */
        const std::vector<double> & uxuz,       /**< [in] filtered u_x * u_z */
        const std::vector<double> & uyuy,       /**< [in] filtered u_y * u_y */
        const std::vector<double> & uyuz,       /**< [in] filtered u_y * u_z */
        const std::vector<double> & uzuz,       /**< [in] filtered u_z * u_z */
        const int Ntime,                        /**< [in] Length of time dimension */
        const int Ndepth,                       /**< [in] Length of depth dimension */
        const int Nlat,                         /**< [in] Length of latitude dimension */
        const int Nlon,                         /**< [in] Length of longitude dimension */
        const std::vector<double> & longitude,  /**< [in] Longitude dimension (1D) */
        const std::vector<double> & latitude,   /**< [in] Latitude dimension (1D) */
        const std::vector<double> & mask        /**< [in] Mask array (2D) to distinguish land from water */
        ) {

    double tau_xx, tau_xy, tau_xz, tau_yy, tau_yz, tau_zz;
    double S_xx,   S_xy,   S_xz,   S_yy,   S_yz,   S_zz;
    double pi_tmp;
    int index, mask_index, Ilat, Ilon;

    for (int Itime = 0; Itime < Ntime; Itime++) {
         for (int Idepth = 0; Idepth < Ndepth; Idepth++) {
            #pragma omp parallel \
            default(none) \
            shared(Itime, Idepth, mask,\
                    uxux, uxuy, uxuz, uyuy, uyuz, uzuz, ux, uy, uz,\
                    energy_transfer, latitude, longitude)\
            private(Ilat, Ilon, index, mask_index,\
                    tau_xx, tau_xy, tau_xz, tau_yy, tau_yz, tau_zz,\
                    S_xx, S_xy, S_xz, S_yy, S_yz, S_zz,\
                    pi_tmp)
            {
                #pragma omp for collapse(2) schedule(dynamic)
                for (Ilat = 0; Ilat < Nlat; Ilat++) {
                    for (Ilon = 0; Ilon < Nlon; Ilon++) {

                        // Convert our four-index to a one-index
                        index = Index(Itime, Idepth, Ilat, Ilon,
                                      Ntime, Ndepth, Nlat, Nlon);
                        mask_index = Index(0,     0,      Ilat, Ilon,
                                           Ntime, Ndepth, Nlat, Nlon);

                        if (mask.at(mask_index) == 1) { // Skip land areas

                            // Compute subfilter-scale stress
                            tau_xx = uxux.at(index) - ux.at(index)*ux.at(index);
                            tau_xy = uxuy.at(index) - ux.at(index)*uy.at(index);
                            tau_xz = uxuz.at(index) - ux.at(index)*uz.at(index);
                            tau_yy = uyuy.at(index) - uy.at(index)*uy.at(index);
                            tau_yz = uyuz.at(index) - uy.at(index)*uz.at(index);
                            tau_zz = uzuz.at(index) - uz.at(index)*uz.at(index);

                            // Compute large-scale strain
                            compute_largescale_strain(
                                    S_xx, S_xy, S_xz, S_yy, S_yz, S_zz, 
                                    ux, uy, uz,
                                    Itime, Idepth, Ilat, Ilon,
                                    Ntime, Ndepth, Nlat, Nlon,
                                    longitude, latitude, mask);

                            pi_tmp = -constants::rho0 * 
                                (         S_xx * tau_xx + S_yy * tau_yy + S_zz * tau_zz
                                          +  2 * ( S_xy * tau_xy + S_xz * tau_xz + S_yz * tau_yz  )
                                );

                        } else {
                            pi_tmp = 0.;
                        }

                        energy_transfer.at(index) = pi_tmp;
                    }
                }
            }
         }
    }
}
