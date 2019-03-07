#include <math.h>
#include <algorithm>
#include <vector>
#include "../functions.hpp"
#include "../constants.hpp"

#ifndef DEBUG
    #define DEBUG 0
#endif

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
        const std::vector<double> & dAreas,     /**< [in] Array of cell areas (2D) (compute_areas())*/
        const double scale,                     /**< [in] The filtering scale */
        const std::vector<double> & mask        /**< [in] Array to distinguish between land and water cells (2D) */
        ) {


    double kA_sum, dist, kern, area;
    int index, mask_index;
    int curr_lon;

    kA_sum     = 0.;
    double coarse_val_tmp = 0.;

    // Grid spacing: assume uniform grid
    double dlat = latitude.at( 1) - latitude.at( 0);
    double dlon = longitude.at(1) - longitude.at(0);

    double dlat_m, dlon_m; 
    int    dlat_N, dlon_N;
    // The spacing (in metres and points) betwee latitude gridpoints
    //   The factor of 2 is diameter->radius 
    dlat_m = dlat * constants::R_earth;
    dlat_N = ceil( (1.1*scale / dlat_m) / 2 );

    int LAT_lb, LAT_ub, LON_lb, LON_ub;
    // Latitude periodicity is a little different / awkward
    //   for the moment, hope it doesn't become an issue
    LAT_lb = std::max(0,      Ilat - dlat_N);
    LAT_ub = std::min(Nlat-1, Ilat + dlat_N);

    double local_scale, delta_lat;

    for (int curr_lat = LAT_lb; curr_lat < LAT_ub; curr_lat++) {

        // Get the longitude grid spacing (in m) at this latitude
        dlon_m = dlon * constants::R_earth * cos(latitude.at(curr_lat));

        // Next determine how far we need to go (since we're already 
        //   at a finite distance in latitude, we don't need to go
        //   the full scale distance in longitude).
        // Essentially, use circular integration regions, not square
        //    this will further reduce the number of cells
        //    required (by up to a factor of 4), which should
        //    improve performance.
        //  The abs in local_scale is to handle the 'comfort zone'
        //    where delta_lat > scale (from the 1.1 factor in dlat_N)
        delta_lat   = constants::R_earth * abs(   latitude.at(Ilat) 
                                                - latitude.at(curr_lat));
        local_scale = sqrt( abs( pow(scale, 2) - pow(delta_lat, 2) ));

        // Now find the appropriate integration region
        //   The factor of 2 is diameter->radius 
        dlon_N = ceil( ( 1.1 * local_scale / dlon_m) / 2 );
        LON_lb = std::max( -Nlon,   Ilon - dlon_N);
        LON_ub = std::min(2*Nlon-1, Ilon + dlon_N);

        for (int LON = LON_lb; LON < LON_ub; LON++) {

            // Handle periodicity
            if (LON < 0) {
                curr_lon = LON + Nlon;
            } else if (LON >= Nlon) {
                curr_lon = LON - Nlon;
            } else {
                curr_lon = LON;
            }

            dist = distance(longitude.at(Ilon),     latitude.at(Ilat),
                            longitude.at(curr_lon), latitude.at(curr_lat));

            kern = kernel(dist, scale);

            index = Index(Itime, Idepth, curr_lat, curr_lon,
                          Ntime, Ndepth, Nlat,     Nlon);

            mask_index = Index(0,     0,      curr_lat, curr_lon,
                               Ntime, Ndepth, Nlat,     Nlon);

            area    = dAreas.at(index);
            kA_sum += kern * area;

            coarse_val_tmp += field.at(index) * kern * area * mask.at(mask_index);

        }
    }

    coarse_val = coarse_val_tmp / kA_sum;
}

