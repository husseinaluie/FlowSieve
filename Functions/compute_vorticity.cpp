#include "../functions.hpp"

void compute_vorticity(
        double * vort_r,           /**< [in] where to store vort_r */
        double * vort_lon,         /**< [in] where to store vort_lon */
        double * vort_lat,         /**< [in] where to store vort_lat */
        const double * u_r,        /**< [in] full u_r for calculation */
        const double * u_lon,      /**< [in] full u_lon for calculation */
        const double * u_lat,      /**< [in] full u_lat for calculation */
        const int Ntime,           /**< [in] Length of time dimension */
        const int Ndepth,          /**< [in] Length of depth dimension */
        const int Nlat,            /**< [in] Length of latitude dimension */
        const int Nlon,            /**< [in] Length of longitude dimension */
        const double * longitude,  /**< [in] Longitude dimension (1D) */
        const double * latitude,   /**< [in] Latitude dimension (1D) */
        const double * mask        /**< [in] Mask array (2D) to distinguish land from water */
        ) {

    double vort_r_tmp, vort_lon_tmp, vort_lat_tmp;
    int index, mask_index;

    for (int Itime = 0; Itime < Ntime; Itime++) {
        for (int Idepth = 0; Idepth < Ndepth; Idepth++) {
            for (int Ilat = 0; Ilat < Nlat; Ilat++) {
                for (int Ilon = 0; Ilon < Nlon; Ilon++) {

                    vort_r_tmp = 0.;
                    vort_lon_tmp = 0.;
                    vort_lat_tmp = 0.;

                    // Convert our four-index to a one-index
                    index = Index(Itime, Idepth, Ilat, Ilon,
                                  Ntime, Ndepth, Nlat, Nlon);
                    mask_index = Index(0,     0,      Ilat, Ilon,
                                       Ntime, Ndepth, Nlat, Nlon);

                    if (mask[mask_index] == 1) { // Skip land areas

                         compute_vorticity_at_point(
                            vort_r_tmp, vort_lon_tmp, vort_lat_tmp,
                            u_r,        u_lon,        u_lat,
                            Ntime, Ndepth, Nlat, Nlon,
                            Itime, Idepth, Ilat, Ilon,
                            longitude, latitude,
                            mask);

                    }
                    vort_r[index]   = vort_r_tmp;
                    vort_lon[index] = vort_lon_tmp;
                    vort_lat[index] = vort_lat_tmp;
                }
            }
        }
    }
}
