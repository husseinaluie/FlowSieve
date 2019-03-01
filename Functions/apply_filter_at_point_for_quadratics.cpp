#include <algorithm>
#include "../functions.hpp"

#ifndef DEBUG
    #define DEBUG 0
#endif

void apply_filter_at_point_for_quadratics(
        double & uxux_tmp,         /**< [in] where to store filtered (u_x)*(u_x) */
        double & uxuy_tmp,         /**< [in] where to store filtered (u_x)*(u_y) */
        double & uxuz_tmp,         /**< [in] where to store filtered (u_x)*(u_z) */
        double & uyuy_tmp,         /**< [in] where to store filtered (u_y)*(u_y) */
        double & uyuz_tmp,         /**< [in] where to store filtered (u_y)*(u_z) */
        double & uzuz_tmp,         /**< [in] where to store filtered (u_z)*(u_z) */
        const double * u_x,       /**< [in] (full) u_x to filter */
        const double * u_y,       /**< [in] (full) u_y to filter */
        const double * u_z,       /**< [in] (full) u_z to filter */
        const int dlon_N,         /**< [in] Maximum longitudinal width of kernel area in cells */
        const int dlat_N,         /**< [in] Maximum latitudinal width of kernel area in cells */
        const int Ntime,          /**< [in] Length of time dimension */
        const int Ndepth,         /**< [in] Length of depth dimension */
        const int Nlat,           /**< [in] Length of latitude dimension */
        const int Nlon,           /**< [in] Length of longitude dimension */
        const int Itime,          /**< [in] Current position in time dimension */
        const int Idepth,         /**< [in] Current position in depth dimension */
        const int Ilat,           /**< [in] Current position in latitude dimension */
        const int Ilon,           /**< [in] Current position in longitude dimension */
        const double * longitude, /**< [in] Longitude dimension (1D) */
        const double * latitude,  /**< [in] Latitude dimension (1D) */
        const double * dAreas,    /**< [in] Array of cell areas (2D) (compute_areas())*/
        const double scale,       /**< [in] The filtering scale */
        const double * mask       /**< [in] Array to distinguish between land and water cells (2D) */
        ) {


    double kA_sum, dist, kern, area;
    int index, mask_index;
    int curr_lat, curr_lon;

    kA_sum  = 0.;
    uxux_tmp = 0.;
    uxuy_tmp = 0.;
    uxuz_tmp = 0.;
    uyuy_tmp = 0.;
    uyuz_tmp = 0.;
    uzuz_tmp = 0.;

    int LAT_lb = std::max( -Nlat,   Ilat - dlat_N);
    int LAT_ub = std::min(2*Nlat-1, Ilat + dlat_N);

    int LON_lb = std::max( -Nlon,   Ilon - dlon_N);
    int LON_ub = std::min(2*Nlon-1, Ilon + dlon_N);

    for (int LAT = LAT_lb; LAT < LAT_ub; LAT++) {

        // Handle periodicity
        if (LAT < 0) {
            curr_lat = LAT + Nlat;
        } else if (LAT > Nlat) {
            curr_lat = LAT - Nlat;
        } else {
            curr_lat = LAT;
        }

        for (int LON = LON_lb; LON < LON_ub; LON++) {

            // Handle periodicity
            if (LON < 0) {
                curr_lon = LON + Nlon;
            } else if (LON > Nlon) {
                curr_lon = LON - Nlon;
            } else {
                curr_lon = LON;
            }

            dist = distance(longitude[Ilon],     latitude[Ilat],
                            longitude[curr_lon], latitude[curr_lat]);

            kern = kernel(dist, scale);

            index = Index(Itime, Idepth, curr_lat, curr_lon,
                          Ntime, Ndepth, Nlat,     Nlon);

            mask_index = Index(0,     0,      curr_lat, curr_lon,
                               Ntime, Ndepth, Nlat,     Nlon);

            area    = dAreas[index];
            kA_sum += kern * area;// * mask[mask_index];

            uxux_tmp += u_x[index] * u_x[index] * kern * area * mask[mask_index];
            uxuy_tmp += u_x[index] * u_y[index] * kern * area * mask[mask_index];
            uxuz_tmp += u_x[index] * u_z[index] * kern * area * mask[mask_index];
            uyuy_tmp += u_y[index] * u_y[index] * kern * area * mask[mask_index];
            uyuz_tmp += u_y[index] * u_z[index] * kern * area * mask[mask_index];
            uzuz_tmp += u_z[index] * u_z[index] * kern * area * mask[mask_index];

        }
    }

    uxux_tmp *= 1. / kA_sum;
    uxuy_tmp *= 1. / kA_sum;
    uxuz_tmp *= 1. / kA_sum;
    uyuy_tmp *= 1. / kA_sum;
    uyuz_tmp *= 1. / kA_sum;
    uzuz_tmp *= 1. / kA_sum;
}

