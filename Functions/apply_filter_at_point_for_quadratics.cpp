#include <math.h>
#include <algorithm>
#include <vector>
#include "../functions.hpp"
#include "../constants.hpp"

void apply_filter_at_point_for_quadratics(
        double & uxux_tmp,                      /**< [in] where to store filtered (u_x)*(u_x) */
        double & uxuy_tmp,                      /**< [in] where to store filtered (u_x)*(u_y) */
        double & uxuz_tmp,                      /**< [in] where to store filtered (u_x)*(u_z) */
        double & uyuy_tmp,                      /**< [in] where to store filtered (u_y)*(u_y) */
        double & uyuz_tmp,                      /**< [in] where to store filtered (u_y)*(u_z) */
        double & uzuz_tmp,                      /**< [in] where to store filtered (u_z)*(u_z) */
        const std::vector<double> & u_x,        /**< [in] (full) u_x to filter */
        const std::vector<double> & u_y,        /**< [in] (full) u_y to filter */
        const std::vector<double> & u_z,        /**< [in] (full) u_z to filter */
        const int Ntime,                        /**< [in] Length of time dimension */
        const int Ndepth,                       /**< [in] Length of depth dimension */
        const int Nlat,                         /**< [in] Length of latitude dimension */
        const int Nlon,                         /**< [in] Length of longitude dimension */
        const int Itime,                        /**< [in] Current position in time dimension */
        const int Idepth,                       /**< [in] Current position in depth dimension */
        const int Ilat,                         /**< [in] Current position in latitude dimension */
        const int Ilon,                         /**< [in] Current position in longitude dimension */
        const std::vector<double> & longitude,  /**< [in] Longitude dimension (1D) */
        const std::vector<double> & latitude,   /**< [in] Latitude dimension (1D) */
        const int LAT_lb,
        const int LAT_ub,
        const std::vector<double> & dAreas,     /**< [in] Array of cell areas (2D) (compute_areas())*/
        const double scale,                     /**< [in] The filtering scale */
        const std::vector<double> & mask,       /**< [in] Array to distinguish between land and water cells (2D) */
        const std::vector<double> * local_kernel    /**< [in] Array of local kernel (if not NULL) */
        ) {


    double kA_sum, dist, kern, area, mask_val;
    double u_x_loc, u_y_loc, u_z_loc;
    int index, mask_index;
    int curr_lon, curr_lat;

    kA_sum  = 0.;
    uxux_tmp = 0.;
    uxuy_tmp = 0.;
    uxuz_tmp = 0.;
    uyuy_tmp = 0.;
    uyuz_tmp = 0.;
    uzuz_tmp = 0.;

    int    LON_lb, LON_ub;
    double lat_at_curr, lat_at_ilat;
    lat_at_ilat = latitude.at(Ilat);

    for (int LAT = LAT_lb; LAT < LAT_ub; LAT++) {

        // Handle periodicity
        if (constants::PERIODIC_Y) {
            if      (LAT <  0   ) { curr_lat = LAT + Nlat; } 
            else if (LAT >= Nlat) { curr_lat = LAT - Nlat; }
            else                  { curr_lat = LAT; }
        } else {
            curr_lat = LAT;
        }
        lat_at_curr = latitude.at(curr_lat);

        get_lon_bounds(LON_lb, LON_ub, longitude, Ilon, 
                lat_at_ilat, lat_at_curr, scale);

        for (int LON = LON_lb; LON < LON_ub; LON++) {

            // Handle periodicity
            if (constants::PERIODIC_X) {
                if      (LON <  0   ) { curr_lon = LON + Nlon; }
                else if (LON >= Nlon) { curr_lon = LON - Nlon; }
                else                  { curr_lon = LON; }
            } else {
                curr_lon = LON;
            }

            index = Index(Itime, Idepth, curr_lat, curr_lon,
                          Ntime, Ndepth, Nlat,     Nlon);

            mask_index = Index(0,     0,      curr_lat, curr_lon,
                               Ntime, Ndepth, Nlat,     Nlon);

            if (local_kernel == NULL) {
                if (constants::CARTESIAN) {
                    double dlat_m = latitude.at( 1) - latitude.at( 0);
                    double dlon_m = longitude.at(1) - longitude.at(0);
                    dist = distance(longitude.at(Ilon),     lat_at_ilat,
                                    longitude.at(curr_lon), lat_at_curr,
                                    dlon_m * Nlon, dlat_m * Nlat);
                } else {
                    dist = distance(longitude.at(Ilon),     lat_at_ilat,
                                    longitude.at(curr_lon), lat_at_curr);
                }
                kern = kernel(dist, scale);
            } else {
                kern = local_kernel->at(mask_index);
            }

            area     = dAreas.at(mask_index);
            kA_sum  += kern * area;
            mask_val = mask.at(mask_index);

            u_x_loc = u_x.at(index);
            u_y_loc = u_y.at(index);
            u_z_loc = u_z.at(index);

            uxux_tmp += u_x_loc * u_x_loc * kern * area * mask_val;
            uxuy_tmp += u_x_loc * u_y_loc * kern * area * mask_val;
            uxuz_tmp += u_x_loc * u_z_loc * kern * area * mask_val;
            uyuy_tmp += u_y_loc * u_y_loc * kern * area * mask_val;
            uyuz_tmp += u_y_loc * u_z_loc * kern * area * mask_val;
            uzuz_tmp += u_z_loc * u_z_loc * kern * area * mask_val;

        }
    }

    uxux_tmp *= 1. / kA_sum;
    uxuy_tmp *= 1. / kA_sum;
    uxuz_tmp *= 1. / kA_sum;
    uyuy_tmp *= 1. / kA_sum;
    uyuz_tmp *= 1. / kA_sum;
    uzuz_tmp *= 1. / kA_sum;
}

