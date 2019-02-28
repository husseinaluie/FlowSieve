#include "../functions.hpp"
#include "../constants.hpp"

void compute_energy_transfer_through_scale(
        double * energy_transfer,  /**< [in] where to store the energy transfer */
        const double * ux,         /**< [in] filtered u_x */
        const double * uy,         /**< [in] filtered u_y */
        const double * uz,         /**< [in] filtered u_z */
        const double * uxux,       /**< [in] filtered u_x * u_x */
        const double * uxuy,       /**< [in] filtered u_x * u_y */
        const double * uxuz,       /**< [in] filtered u_x * u_z */
        const double * uyuy,       /**< [in] filtered u_y * u_y */
        const double * uyuz,       /**< [in] filtered u_y * u_z */
        const double * uzuz,       /**< [in] filtered u_z * u_z */
        const int Ntime,           /**< [in] Length of time dimension */
        const int Ndepth,          /**< [in] Length of depth dimension */
        const int Nlat,            /**< [in] Length of latitude dimension */
        const int Nlon,            /**< [in] Length of longitude dimension */
        const double * longitude,  /**< [in] Longitude dimension (1D) */
        const double * latitude,   /**< [in] Latitude dimension (1D) */
        const double * mask        /**< [in] Mask array (2D) to distinguish land from water */
        ) {

    double tau_xx, tau_xy, tau_xz, tau_yy, tau_yz, tau_zz;
    double S_xx,   S_xy,   S_xz,   S_yy,   S_yz,   S_zz;
    double pi_tmp;
    int index, mask_index;

    for (int Itime = 0; Itime < Ntime; Itime++) {
        for (int Idepth = 0; Idepth < Ndepth; Idepth++) {
            for (int Ilat = 0; Ilat < Nlat; Ilat++) {
                for (int Ilon = 0; Ilon < Nlon; Ilon++) {

                    // Convert our four-index to a one-index
                    index = Index(Itime, Idepth, Ilat, Ilon,
                                  Ntime, Ndepth, Nlat, Nlon);
                    mask_index = Index(0,     0,      Ilat, Ilon,
                                       Ntime, Ndepth, Nlat, Nlon);

                    if (mask[mask_index] == 1) { // Skip land areas

                        // Compute subfilter-scale stress
                        tau_xx = uxux[index] - ux[index]*ux[index];
                        tau_xy = uxuy[index] - ux[index]*uy[index];
                        tau_xz = uxuz[index] - ux[index]*uz[index];
                        tau_yy = uyuy[index] - uy[index]*uy[index];
                        tau_yz = uyuz[index] - uy[index]*uz[index];
                        tau_zz = uzuz[index] - uz[index]*uz[index];

                        // Compute large-scale strain
                        compute_largescale_strain(
                                S_xx, S_xy, S_xz, S_yy, S_yz, S_zz, 
                                ux, uy, uz,
                                Ntime, Ndepth, Nlat, Nlon,
                                Itime, Idepth, Ilat, Ilon,
                                longitude, latitude, mask);

                        pi_tmp = -constants::rho0 * 
                            (         S_xx * tau_xx + S_yy * tau_yy + S_zz * tau_zz
                             +  2 * ( S_xy * tau_xy + S_xz * tau_xz + S_yz * tau_yz  )
                            );

                    } else {
                        pi_tmp = 0.;
                    }

                    energy_transfer[index] = pi_tmp;
                }
            }
        }
    }
}
