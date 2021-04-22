#include <math.h>
#include <algorithm>
#include <vector>
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
 * @param[in]       u_x,u_y,u_z             fields to filter
 * @param[in]       Ntime,Ndepth,Nlat,Nlon  length of time dimension
 * @param[in]       Itime,Idepth,Ilat,Ilon  current position in time dimension
 * @param[in]       longitude,latitude      grid vectors (lon,lat)
 * @param[in]       LAT_lb,LAT_ub           lower/upper boundd on latitude for kernel
 * @param[in]       dAreas                  array of cell areas (2D - lat,lon)
 * @param[in]       scale                   filtering scale
 * @param[in]       mask                    array to distinguish land from water
 * @param[in]       local_kernel            pointer to pre-computed kernel (NULL indicates not provided)
 */
void apply_filter_at_point_for_quadratics(
        double & uxux_tmp,
        double & uxuy_tmp,
        double & uxuz_tmp,
        double & uyuy_tmp,
        double & uyuz_tmp,
        double & uzuz_tmp,
        const std::vector<double> & u_x,
        const std::vector<double> & u_y,
        const std::vector<double> & u_z,
        const int Ntime,
        const int Ndepth,
        const int Nlat,
        const int Nlon,
        const int Itime,
        const int Idepth,
        const int Ilat
        const int Ilon,
        const std::vector<double> & longitude,
        const std::vector<double> & latitude,
        const int LAT_lb,
        const int LAT_ub,
        const std::vector<double> & dAreas,
        const double scale,
        const std::vector<bool>   & mask,
        const std::vector<double> * local_kernel
        ) {


    double kA_sum, dist, kern, area, mask_val;
    double u_x_loc, u_y_loc, u_z_loc;
    int index, area_index;
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

            area_index = Index(0,     0,      curr_lat, curr_lon,
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
                kern = local_kernel->at(area_index);
            }


            // If cell is water, or if we're not deforming around land, then include the cell area in the integral
            mask_val = ( mask.at(index) or not(constants::DEFORM_AROUND_LAND) ) ? 1. : 0.;
            area     = dAreas.at(area_index);
            kA_sum  += kern * area * mask_val;

            // If the cell is water, keep the value, otherwise zero it out
            mask_val = mask.at(index) ? 1. : 0.;

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

