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
 * @param[in]       source_data             dataset class instance containing data (Psi, Phi, etc)
 * @param[in]       Itime,Idepth,Ilat,Ilon  current position in time dimension
 * @param[in]       LAT_lb,LAT_ub           lower/upper boundd on latitude for kernel
 * @param[in]       scale                   filtering scale
 * @param[in]       use_mask                array of booleans indicating whether or not to use mask (i.e. zero out land) or to use the array value
 * @param[in]       local_kernel            pre-computed kernel (NULL indicates not provided)
 * @param[in]       weight                  pointer to spatial weight (i.e. rho) (NULL indicates not provided)
 *
 */
void apply_filter_at_point(
        std::vector<double*> & coarse_vals,
        const std::vector<const std::vector<double>*> & fields,
        const dataset & source_data,
        const int Itime,
        const int Idepth,
        const int Ilat,
        const int Ilon,
        const int LAT_lb,
        const int LAT_ub,
        const double scale,
        const std::vector<bool> & use_mask,
        const std::vector<double> & local_kernel,
        const std::vector<double> * weight
        ) {

    assert(coarse_vals.size() == fields.size());
    const size_t Nfields = fields.size();

    const std::vector<double>   &latitude   = source_data.latitude,
                                &longitude  = source_data.longitude,
                                &dAreas     = source_data.areas;

    const std::vector<bool> &mask = source_data.mask;

    const int   Ntime   = source_data.Ntime,
                Ndepth  = source_data.Ndepth,
                Nlat    = source_data.Nlat,
                Nlon    = source_data.Nlon;

    double kern, area, loc_val, loc_weight;
    size_t index, kernel_index;

    double  kA_sum   = 0.;
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

            index = Index(Itime, Idepth, curr_lat, curr_lon, Ntime, Ndepth, Nlat, Nlon);

            if ( (constants::UNIFORM_LON_GRID) and (constants::FULL_LON_SPAN) and (constants::PERIODIC_X) ) {
                // In this case, we can re-use the kernel from a previous Ilon value by just shifting our indices
                //  This cuts back on the most computation-heavy part of the code (computing kernels / distances)
                kernel_index = Index(0, 0, curr_lat, ( (LON - Ilon) % Nlon + Nlon ) % Nlon, Ntime, Ndepth, Nlat, Nlon);
            } else {
                kernel_index = Index(0, 0, curr_lat, curr_lon, Ntime, Ndepth, Nlat, Nlon);
            }
            #if DEBUG >= 1
            kern = local_kernel.at(kernel_index);
            area = dAreas.at(kernel_index);
            bool is_water = mask.at(index);
            #else
            kern = local_kernel[kernel_index];
            area = dAreas[kernel_index];
            bool is_water = mask[index];
            #endif
            loc_weight = kern * area;

            // If cell is water, or if we're not deforming around land, then include the cell area in the denominator
            if ( not(constants::DEFORM_AROUND_LAND) or is_water ) { kA_sum += loc_weight; }

            // If we are not using the mask, or if we are on a water cell, include the value in the numerator
            if (weight != NULL) { loc_weight *= weight->at(index); }
            if ( is_water ) {
                for (size_t II = 0; II < Nfields; ++II) {
                    #if DEBUG >= 1
                    loc_val = fields.at(II)->at(index);
                    tmp_vals.at(II) += loc_val * loc_weight;
                    #else
                    loc_val = fields[II]->at(index);
                    tmp_vals[II] += loc_val * loc_weight;
                    #endif
                }
            }
        }
    }

    // On the off chance that the kernel was null (size zero), just return zero
    for (size_t II = 0; II < Nfields; ++II) {
        if (coarse_vals.at(II) != NULL) { *(coarse_vals.at(II)) = (kA_sum == 0) ? 0. : tmp_vals.at(II) / kA_sum; }
    }
}
