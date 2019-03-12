#include <math.h>
#include <vector>
#include "../functions.hpp"
#include "../constants.hpp"
#include "../differentiation_tools.hpp"

void  compute_baroclinic_transfer(
    std::vector<double> & baroclinic_transfer,  /**< [in] Where to store the computed values*/
    const std::vector<double> & coarse_vort_r,  /**< [in] Full vorticity (r   component) */ 
    const std::vector<double> & coarse_vort_lon,/**< [in] Full vorticity (lon component) */ 
    const std::vector<double> & coarse_vort_lat,/**< [in] Full vorticity (lat component) */
    const std::vector<double> & coarse_rho,     /**< [in] Coarse density field */ 
    const std::vector<double> & coarse_p,       /**< [in] Coarse pressure field */
    const int Ntime,                            /**< [in] Length of time dimension */
    const int Ndepth,                           /**< [in] Length of depth dimension */
    const int Nlat,                             /**< [in] Length of latitude dimension */
    const int Nlon,                             /**< [in] Length of longitude dimension */
    const std::vector<double> & longitude,      /**< [in] Longitude dimension (1D) */
    const std::vector<double> & latitude,       /**< [in] Latitude dimension (1D) */
    const std::vector<double> & mask            /**< [in] Mask array (2D) to distinguish land from water */
    ) {

    // For the moment, only use vort_r
    double drhodlat, drhodlon, dpdlat, dpdlon;
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

                    if (mask.at(mask_index) == 1) { // Skip land areas

                        // We need a few derivatives
                        drhodlat = spher_derivative_at_point(
                                coarse_rho, latitude, "lat",
                                Itime, Idepth, Ilat, Ilon,
                                Ntime, Ndepth, Nlat, Nlon,
                                mask);
                        
                        dpdlat = spher_derivative_at_point(
                                coarse_p, latitude, "lat",
                                Itime, Idepth, Ilat, Ilon,
                                Ntime, Ndepth, Nlat, Nlon,
                                mask);
                        
                        drhodlon = spher_derivative_at_point(
                                coarse_rho, longitude, "lon",
                                Itime, Idepth, Ilat, Ilon,
                                Ntime, Ndepth, Nlat, Nlon,
                                mask);

                        dpdlon = spher_derivative_at_point(
                                coarse_p, longitude, "lon",
                                Itime, Idepth, Ilat, Ilon,
                                Ntime, Ndepth, Nlat, Nlon,
                                mask);

                        baroclinic_transfer.at(index) = 
                              coarse_vort_r.at(index) * ( drhodlon * dpdlat  -  drhodlat * dpdlon ) 
                            / ( coarse_rho.at(index) * pow(constants::R_earth,2) * cos(latitude.at(Ilat)) );

                    }
                }
            }
        }
    }

}
