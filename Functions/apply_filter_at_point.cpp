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

    double  kA_sum   = 0.,
            mask_val = 0.;
    std::vector<double> tmp_vals(Nfields);

    int curr_lon, curr_lat, LON_lb, LON_ub;

    double lat_at_curr;
    const double lat_at_ilat = latitude.at(Ilat);

    for (int LAT = LAT_lb; LAT < LAT_ub; LAT++) {

        // Handle periodicity if necessary
        if (constants::PERIODIC_Y) { curr_lat = ( LAT % Nlat + Nlat ) % Nlat; }
        else                       { curr_lat = LAT; }
        lat_at_curr = latitude.at(curr_lat);

        get_lon_bounds(LON_lb, LON_ub, longitude, Ilon, lat_at_ilat, lat_at_curr, scale);
        for (int LON = LON_lb; LON < LON_ub; LON++ ) {

            // Handle periodicity if necessary
            if (constants::PERIODIC_X) { curr_lon = ( LON % Nlon + Nlon ) % Nlon; }
            else                       { curr_lon = LON; }

            index = Index(Itime, Idepth, curr_lat, curr_lon,
                          Ntime, Ndepth, Nlat,     Nlon);
            area_index = Index(0,     0,      curr_lat, curr_lon,
                               Ntime, Ndepth, Nlat,     Nlon);
            area = dAreas.at(area_index);

            if (local_kernel == NULL) {
                // If no pre-computed kernel was provided, then compute it now.
                //  NOTE: This is generally very inefficient. Better to compute
                //        ahead of time
                if (constants::CARTESIAN) {
                    dist = distance(longitude.at(Ilon),     lat_at_ilat,
                                    longitude.at(curr_lon), lat_at_curr,
                                    ( longitude.at(1) - longitude.at(0) ) * Nlon, 
                                    ( latitude.at( 1) - latitude.at( 0) ) * Nlat);
                } else {
                    dist = distance(longitude.at(Ilon),     lat_at_ilat,
                                    longitude.at(curr_lon), lat_at_curr);
                }
                kern = kernel(dist, scale);
            } else {
                size_t kernel_index;
                // If the kernel was provided, then just grab the appropraite point.
                if ( (constants::UNIFORM_LON_GRID) and (constants::FULL_LON_SPAN) and (constants::PERIODIC_X ) ) {
                    // In this case, we can re-use the kernel from a previous Ilon value by just shifting our indices
                    //  This cuts back on the most computation-heavy part of the code (computing kernels / distances)
                    kernel_index = Index(0,     0,      curr_lat, ( (LON - Ilon) % Nlon + Nlon ) % Nlon,
                                         Ntime, Ndepth, Nlat,     Nlon);
                } else {
                    kernel_index = area_index;
                }
                kern = local_kernel->at(kernel_index);
            }

            // If cell is water, or if we're not deforming around land, then include the cell area in the denominator
            if ( mask.at(index) or not(constants::DEFORM_AROUND_LAND) ) {
                kA_sum += kern * area;
            }

            for (size_t II = 0; II < Nfields; ++II) {
                // If we are not using the mask, or if we are on a water cell, include the value in the numerator
                if ( mask.at(index) ) {
                    loc_val = fields.at(II)->at(index);
                    if (weight != NULL) { loc_val *= weight->at(index); }
                    tmp_vals.at(II) += loc_val * kern * area;
                }
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
