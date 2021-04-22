#include <math.h>
#include <algorithm>
#include <vector>
#include <cassert>
#include "../functions.hpp"
#include "../constants.hpp"

/*!
 * \brief Compute filtered field at a single point
 *
 * Computes the integral of the provided field with the
 * kernel().
 *
 * @param[in,out]   coarse_val              where to store filtered value
 * @param[in]       fields                  fields to filter
 * @param[in]       Ntime,Ndepth,Nlat,Nlon  length of time dimension
 * @param[in]       Itime,Idepth,Ilat,Ilon  current position in time dimension
 * @param[in]       longitude,latitude      grid vectors (lon,lat)
 * @param[in]       LAT_lb,LAT_ub           lower/upper boundd on latitude for kernel
 * @param[in]       dAreas                  array of cell areas (2D - lat,lon)
 * @param[in]       scale                   filtering scale
 * @param[in]       mask                    array to distinguish land from water
 * @param[in]       use_mask                array of booleans indicating whether or not to use mask (i.e. zero out land) or to use the array value
 * @param[in]       local_kernel            pointer to pre-computed kernel (NULL indicates not provided)
 * @param[in]       weight                  pointer to spatial weight (i.e. rho) (NULL indicates not provided)
 *
 */
void apply_filter_at_point(
        std::vector<double*> & coarse_vals,
        const std::vector<const std::vector<double>*> & fields,
        const int Ntime,
        const int Ndepth,
        const int Nlat,
        const int Nlon,
        const int Itime,
        const int Idepth,
        const int Ilat,
        const int Ilon,
        const std::vector<double> & longitude,
        const std::vector<double> & latitude,
        const int LAT_lb,
        const int LAT_ub,
        const std::vector<double> & dAreas,
        const double scale,
        const std::vector<bool> & mask,
        const std::vector<bool> & use_mask,
        const std::vector<double> * local_kernel,
        const std::vector<double> * weight
        ) {

    assert(coarse_vals.size() == fields.size());
    const size_t Nfields = fields.size();

    double dist, kern, area, loc_val;
    size_t index, area_index;
    int curr_lon, curr_lat;

    double kA_sum = 0.;
    std::vector<double> tmp_vals(Nfields);
    double mask_val = 0.;

    double dlat_m, dlon_m; 
    int LON_lb, LON_ub;

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
            area_index = Index(0,     0,      curr_lat, curr_lon,
                               Ntime, Ndepth, Nlat,     Nlon);

            if (local_kernel == NULL) {
                if (constants::CARTESIAN) {
                    dlat_m = latitude.at( 1) - latitude.at( 0);
                    dlon_m = longitude.at(1) - longitude.at(0);
                    dist = distance(longitude.at(Ilon),     lat_at_ilat,
                                    longitude.at(curr_lon), lat_at_curr,
                                    dlon_m * Nlon, dlat_m * Nlat);
                } else {
                    dist = distance(longitude.at(Ilon),     lat_at_ilat,
                                    longitude.at(curr_lon), lat_at_curr);
                }
                kern = kernel(dist, scale);
            } else {
                kern = local_kernel->at(area_index);
            }

            // If cell is water, or if we're not deforming around land, then include the cell area in the integral
            mask_val = ( mask.at(index) or not(constants::DEFORM_AROUND_LAND) ) ? 1. : 0.;

            area    = dAreas.at(area_index);
            kA_sum += kern * area * mask_val;

            for (size_t II = 0; II < Nfields; ++II) {
                // If we are not using the mask, or if we are on a water cell, include the value in the integrak
                if ( mask.at(index) ) {
                    mask_val == 1.;
                } else if ( constants::DEFORM_AROUND_LAND or use_mask.at(II) ) {
                    mask_val == 0;
                }
                //mask_val = ( not(use_mask.at(II)) or mask.at(index) ) ? 1. : 0.;

                loc_val = fields.at(II)->at(index);
                if (weight != NULL) { loc_val *= weight->at(index); }

                tmp_vals.at(II) += loc_val * kern * area * mask_val;
            }

        }
    }

    // On the off chance that the kernel was null (size zero), just return zero
    for (size_t II = 0; II < Nfields; ++II) {
        if (coarse_vals.at(II) != NULL) {
            *(coarse_vals.at(II)) = (kA_sum == 0) ? 0. : tmp_vals.at(II) / kA_sum;
        }
    }
}
