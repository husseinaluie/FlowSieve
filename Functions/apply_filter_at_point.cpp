#include <algorithm>
#include "../functions.hpp"

#ifndef DEBUG
    #define DEBUG 0
#endif

void apply_filter_at_point(
        double & u_x_tmp,         /**< [in] where to store filtered u_x */
        double & u_y_tmp,         /**< [in] where to store filtered u_y */
        double & u_z_tmp,         /**< [in] where to store filtered u_z */
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
    u_x_tmp = 0.;
    u_y_tmp = 0.;
    u_z_tmp = 0.;

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

            u_x_tmp += u_x[index] * kern * area * mask[mask_index];
            u_y_tmp += u_y[index] * kern * area * mask[mask_index];
            u_z_tmp += u_z[index] * kern * area * mask[mask_index];

        }
    }

    u_x_tmp = u_x_tmp / kA_sum;
    u_y_tmp = u_y_tmp / kA_sum;
    u_z_tmp = u_z_tmp / kA_sum;
}

