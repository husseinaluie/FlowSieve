#include <math.h>
#include <algorithm>
#include <vector>
#include "../functions.hpp"
#include "../constants.hpp"

void apply_filter_at_point(
        double & coarse_val,                    /**< [in] where to store filtered value */
        const std::vector<double> & field,      /**< [in] field to filter */
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
        const bool use_mask,                    /**< [in] Whether or not to mask (zero) out land cells when integrating */
        const std::vector<double> * local_kernel    /**< [in] Array of local_kernel (if not NULL) */
        ) {

    double dist, kern, area;
    int index, mask_index;
    int curr_lon, curr_lat;

    double kA_sum = 0.;
    double coarse_val_tmp = 0.;
    double mask_val = 0.;

    // Grid spacing: assume uniform grid
    const double dlat = latitude.at( 1) - latitude.at( 0);
    const double dlon = longitude.at(1) - longitude.at(0);

    double dlat_m, dlon_m; 
    int LON_lb, LON_ub;

    // The spacing (in metres and points) betwee latitude gridpoints
    //   The factor of 2 is diameter->radius 
    if (constants::CARTESIAN) { dlat_m = dlat; } 
    else                      { dlat_m = dlat * constants::R_earth; }

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

        // Get the longitude grid spacing (in m) at this latitude
        if (constants::CARTESIAN) { dlon_m = dlon; } 
        else { dlon_m = dlon * constants::R_earth * cos(lat_at_curr); }

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

            area      = dAreas.at(mask_index);
            kA_sum   += kern * area;
            mask_val  = use_mask ? mask.at(mask_index) : 1.;

            coarse_val_tmp += field.at(index) * kern * area * mask_val;

        }
    }
    if (kA_sum != 0) { coarse_val = coarse_val_tmp / kA_sum; }
    else { coarse_val = 0.; }
}
