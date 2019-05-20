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
        const std::vector<double> & dAreas,     /**< [in] Array of cell areas (2D) (compute_areas())*/
        const double scale,                     /**< [in] The filtering scale */
        const std::vector<double> & mask,       /**< [in] Array to distinguish between land and water cells (2D) */
        const bool use_mask,                    /**< [in] Whether or not to mask (zero) out land cells when integrating */
        const std::vector<double> * distances   /**< [in] Array of distances (if not NULL) */
        ) {


    double dist, kern, area;
    int index, mask_index;
    int curr_lon, curr_lat;

    double kA_sum = 0.;
    double coarse_val_tmp = 0.;
    double mask_val = 0.;
    const double KernPad = constants::KernPad;

    // Grid spacing: assume uniform grid
    const double dlat = latitude.at( 1) - latitude.at( 0);
    const double dlon = longitude.at(1) - longitude.at(0);

    double dlat_m, dlon_m; 
    int    dlat_N, dlon_N;
    // The spacing (in metres and points) betwee latitude gridpoints
    //   The factor of 2 is diameter->radius 
    if (constants::CARTESIAN) { dlat_m = dlat; } 
    else                      { dlat_m = dlat * constants::R_earth; }

    if (KernPad < 0) { dlat_N = Nlat; } 
    else             { dlat_N = ceil( ( KernPad * scale / dlat_m ) / 2.); }

    dlat_N = std::min(Nlat, dlat_N);

    int LAT_lb, LAT_ub, LON_lb, LON_ub;
    // Latitude periodicity is a little different / awkward
    //   for the moment, hope it doesn't become an issue
    if (constants::PERIODIC_Y) {
        LAT_lb = Ilat - dlat_N;
        LAT_ub = Ilat + dlat_N;
        if (LAT_lb + Nlat < LAT_ub) { LAT_ub = LAT_lb + Nlat; }
    } else {
        LAT_lb = std::max(0,    Ilat - dlat_N);
        LAT_ub = std::min(Nlat, Ilat + dlat_N);
    }

    double local_scale, delta_lat;
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

        // Next determine how far we need to go (since we're already 
        //   at a finite distance in latitude, we don't need to go
        //   the full scale distance in longitude).
        // Essentially, use 'circular' integration regions, not square
        //    this will further reduce the number of cells
        //    required (by up to a factor of 4), which should
        //    improve performance.
        //  The abs in local_scale is to handle the 'comfort zone'
        //    where delta_lat > scale (from the KernPad factor in dlat_N)
        if (constants::CARTESIAN) { delta_lat = lat_at_ilat - lat_at_curr; }
        else { delta_lat = constants::R_earth * ( lat_at_ilat - lat_at_curr ); }
        local_scale = sqrt( fabs( scale*scale - delta_lat*delta_lat ));

        // Now find the appropriate integration region
        //   The factor of 2 is diameter->radius 
        if (KernPad < 0) { dlon_N = Nlon; }
        else { dlon_N = ceil( ( KernPad * local_scale / dlon_m ) / 2.); }
        dlon_N = std::min(Nlon, dlon_N);

        if (constants::PERIODIC_X) {
            LON_lb = Ilon - dlon_N;
            LON_ub = Ilon + dlon_N;
            if (LON_ub - LON_lb > Nlon) { LON_ub = LON_lb + Nlon; }
        } else {
            LON_lb = std::max(0,    Ilon - dlon_N);
            LON_ub = std::min(Nlon, Ilon + dlon_N);
        }

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

            if (distances == NULL) {
                if (constants::CARTESIAN) {
                    dist = distance(longitude.at(Ilon),     lat_at_ilat,
                                    longitude.at(curr_lon), lat_at_curr,
                                    dlon_m * Nlon, dlat_m * Nlat);
                } else {
                    dist = distance(longitude.at(Ilon),     lat_at_ilat,
                                    longitude.at(curr_lon), lat_at_curr);
                }
            } else {
                dist = distances->at(mask_index);
            }

            kern = kernel(dist, scale);

            area    = dAreas.at(mask_index);
            kA_sum += kern * area;

            if (use_mask) { mask_val = mask.at(mask_index); }
            else          { mask_val = 1.; }

            coarse_val_tmp += field.at(index) * kern * area * mask_val;

        }
    }
    coarse_val = coarse_val_tmp / kA_sum;
}
