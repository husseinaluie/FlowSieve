#include <math.h>
#include <algorithm>
#include <vector>
#include <cassert>
#include "../functions.hpp"
#include "../constants.hpp"

/*!
 * \brief Compute filtered quadratic velocities at a single point
 *
 * Computes the integral of each quadratic Cartesian velocity
 * with the kernel().
 *
 * In particular, the quadratic terms being filtered are:
 * \$u_xu_x\$, \$u_xu_y\$, \$u_xu_z\$, \$u_yu_y\$, \$u_yu_z\$, \$u_zu_z\$
 *
 * @param[in,out]   uxux_tmp                where to store filtered (u_x)*(u_x)
 * @param[in,out]   uxuy_tmp                where to store filtered (u_x)*(u_y)
 * @param[in,out]   uxuz_tmp                where to store filtered (u_x)*(u_z)
 * @param[in,out]   uyuy_tmp                where to store filtered (u_y)*(u_y)
 * @param[in,out]   uyuz_tmp                where to store filtered (u_y)*(u_z)
 * @param[in,out]   uzuz_tmp                where to store filtered (u_z)*(u_z)
 * @param[in,out]   vort_ux_tmp             where to store filtered (vort_r)*(u_x)
 * @param[in,out]   vort_uy_tmp             where to store filtered (vort_r)*(u_y)
 * @param[in,out]   vort_uz_tmp             where to store filtered (vort_r)*(u_z)
 * @param[in]       u_x,u_y,u_z             fields to filter
 * @param[in]       vort_r                  vorticity field to filter
 * @param[in]       source_data             dataset class instance containing data (Psi, Phi, etc)
 * @param[in]       Itime,Idepth,Ilat,Ilon  current position in time dimension
 * @param[in]       LAT_lb,LAT_ub           lower/upper boundd on latitude for kernel
 * @param[in]       scale                   filtering scale
 * @param[in]       local_kernel            pre-computed kernel (NULL indicates not provided)
 */
void apply_filter_at_point_for_quadratics(
        double & uxux_tmp,
        double & uxuy_tmp,
        double & uxuz_tmp,
        double & uyuy_tmp,
        double & uyuz_tmp,
        double & uzuz_tmp,
        double & vort_ux_tmp,
        double & vort_uy_tmp,
        double & vort_uz_tmp,
        const std::vector<double> & u_x,
        const std::vector<double> & u_y,
        const std::vector<double> & u_z,
        const std::vector<double> & vort_r,
        const dataset & source_data,
        const int Itime,
        const int Idepth,
        const int Ilat,
        const int Ilon,
        const int LAT_lb,
        const int LAT_ub,
        const double scale,
        const std::vector<double> & local_kernel
        ) {



    double  kern = 0, area = 0, kA_sum = 0, local_weight = 0,
            u_x_loc = 0, u_y_loc = 0, u_z_loc = 0, vort_r_loc = 0;
    size_t index, kernel_index;

    const std::vector<double>   &latitude   = source_data.latitude,
                                &longitude  = source_data.longitude,
                                &dAreas     = source_data.areas;

    const std::vector<bool> &mask = source_data.mask;

    const int   Ntime   = source_data.Ntime,
                Ndepth  = source_data.Ndepth,
                Nlat    = source_data.Nlat,
                Nlon    = source_data.Nlon;

    // Zero out the coarse values before we start accumulating (integrating) over space
    uxux_tmp = 0.;
    uxuy_tmp = 0.;
    uxuz_tmp = 0.;
    uyuy_tmp = 0.;
    uyuz_tmp = 0.;
    uzuz_tmp = 0.;

    vort_ux_tmp = 0.;
    vort_uy_tmp = 0.;
    vort_uz_tmp = 0.;


    int    curr_lon, curr_lat, LON_lb, LON_ub;
    double lat_at_curr;
    const double lat_at_ilat = latitude.at(Ilat);

    for (int LAT = LAT_lb; LAT < LAT_ub; LAT++) {

        // Handle periodicity if necessary
        if (constants::PERIODIC_Y) { curr_lat = ( LAT % Nlat + Nlat ) % Nlat; }
        else                       { curr_lat = LAT; }
        lat_at_curr = latitude.at(curr_lat);

        get_lon_bounds(LON_lb, LON_ub, longitude, Ilon, lat_at_ilat, lat_at_curr, scale);

        for (int LON = LON_lb; LON < LON_ub; LON++) {

            // Handle periodicity if necessary
            if (constants::PERIODIC_X) { curr_lon = ( LON % Nlon + Nlon ) % Nlon; }
            else                       { curr_lon = LON; }

            index = Index(Itime, Idepth, curr_lat, curr_lon, Ntime, Ndepth, Nlat, Nlon);

            if ( (constants::UNIFORM_LON_GRID) and (constants::FULL_LON_SPAN) and (constants::PERIODIC_X ) ) {
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
            local_weight = kern * area;

            // If cell is water, or if we're not deforming around land, then include the cell area in the denominator
            //      i.e. treat land cells as zero velocity, unless we're deforming around land
            if ( not(constants::DEFORM_AROUND_LAND) or is_water ) { kA_sum += local_weight; }

            // If the cell is water, add to the numerator
            if ( is_water ) {
                #if DEBUG >= 1
                u_x_loc     = u_x.at(index);
                u_y_loc     = u_y.at(index);
                u_z_loc     = u_z.at(index);
                vort_r_loc  = vort_r.at(index);
                #else
                u_x_loc     = u_x[index];
                u_y_loc     = u_y[index];
                u_z_loc     = u_z[index];
                vort_r_loc  = vort_r[index];
                #endif

                uxux_tmp += u_x_loc * u_x_loc * local_weight;
                uxuy_tmp += u_x_loc * u_y_loc * local_weight;
                uxuz_tmp += u_x_loc * u_z_loc * local_weight;
                uyuy_tmp += u_y_loc * u_y_loc * local_weight;
                uyuz_tmp += u_y_loc * u_z_loc * local_weight;
                uzuz_tmp += u_z_loc * u_z_loc * local_weight;

                vort_ux_tmp += vort_r_loc * u_x_loc * local_weight;
                vort_uy_tmp += vort_r_loc * u_y_loc * local_weight;
                vort_uz_tmp += vort_r_loc * u_z_loc * local_weight;
            }
        }
    }

    uxux_tmp = (kA_sum == 0) ? 0. : uxux_tmp / kA_sum;
    uxuy_tmp = (kA_sum == 0) ? 0. : uxuy_tmp / kA_sum;
    uxuz_tmp = (kA_sum == 0) ? 0. : uxuz_tmp / kA_sum;
    uyuy_tmp = (kA_sum == 0) ? 0. : uyuy_tmp / kA_sum;
    uyuz_tmp = (kA_sum == 0) ? 0. : uyuz_tmp / kA_sum;
    uzuz_tmp = (kA_sum == 0) ? 0. : uzuz_tmp / kA_sum;

    vort_ux_tmp *= (kA_sum == 0) ? 0. : vort_ux_tmp / kA_sum;
    vort_uy_tmp *= (kA_sum == 0) ? 0. : vort_uy_tmp / kA_sum;
    vort_uz_tmp *= (kA_sum == 0) ? 0. : vort_uz_tmp / kA_sum;
}

