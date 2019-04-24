#include <math.h>
#include <algorithm>
#include <vector>
#include "../functions.hpp"
#include "../constants.hpp"

void apply_filter_at_point_for_quadratics(
        double & uxux_tmp,                      /**< [in] where to store filtered (u_x)*(u_x) */
        double & uxuy_tmp,                      /**< [in] where to store filtered (u_x)*(u_y) */
        double & uxuz_tmp,                      /**< [in] where to store filtered (u_x)*(u_z) */
        double & uyuy_tmp,                      /**< [in] where to store filtered (u_y)*(u_y) */
        double & uyuz_tmp,                      /**< [in] where to store filtered (u_y)*(u_z) */
        double & uzuz_tmp,                      /**< [in] where to store filtered (u_z)*(u_z) */
        const std::vector<double> & u_x,        /**< [in] (full) u_x to filter */
        const std::vector<double> & u_y,        /**< [in] (full) u_y to filter */
        const std::vector<double> & u_z,        /**< [in] (full) u_z to filter */
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
        const std::vector<double> * distances   /**< [in] Array of distances (if not NULL) */
        ) {


    double kA_sum, dist, kern, area, mask_val;
    double u_x_loc, u_y_loc, u_z_loc;
    int index, mask_index;
    int curr_lon, curr_lat;
    const double KernPad = constants::KernPad;

    kA_sum  = 0.;
    uxux_tmp = 0.;
    uxuy_tmp = 0.;
    uxuz_tmp = 0.;
    uyuy_tmp = 0.;
    uyuz_tmp = 0.;
    uzuz_tmp = 0.;

    // Grid spacing: assume uniform grid
    double dlat = latitude.at( 1) - latitude.at( 0);
    double dlon = longitude.at(1) - longitude.at(0);

    double dlat_m, dlon_m; 
    int    dlat_N, dlon_N;
    // The spacing (in metres and points) betwee latitude gridpoints
    //   The factor of 2 is diameter->radius 
    #if CARTESIAN
    dlat_m = dlat;
    #elif not(CARTESIAN)
    dlat_m = dlat * constants::R_earth;
    #endif
    if (KernPad < 0) {
        dlat_N = Nlat;
    } else {
        dlat_N = ceil( ( KernPad * scale / dlat_m ) / 2.);
    }
    dlat_N = std::min(Nlat, dlat_N);

    int LAT_lb, LAT_ub, LON_lb, LON_ub;
    // Latitude periodicity is a little different / awkward
    //   for the moment, hope it doesn't become an issue
    #if PERIODIC_Y
    LAT_lb = Ilat - dlat_N;
    LAT_ub = Ilat + dlat_N;
    if (LAT_lb + Nlat < LAT_ub) { LAT_ub = LAT_lb + Nlat; }
    #else
    LAT_lb = std::max(0,    Ilat - dlat_N);
    LAT_ub = std::min(Nlat, Ilat + dlat_N);
    #endif

    double local_scale, delta_lat;
    double lat_at_curr, lat_at_ilat;
    lat_at_ilat = latitude.at(Ilat);

    for (int LAT = LAT_lb; LAT < LAT_ub; LAT++) {

        #if PERIODIC_Y
        // Handle periodicity
        if      (LAT <  0   ) { curr_lat = LAT + Nlat; } 
        else if (LAT >= Nlat) { curr_lat = LAT - Nlat; }
        else                  { curr_lat = LAT; }
        #else
        curr_lat = LAT;
        #endif
        lat_at_curr = latitude.at(curr_lat);

        // Get the longitude grid spacing (in m) at this latitude
        #if CARTESIAN
        dlon_m = dlon;
        #elif not(CARTESIAN)
        dlon_m = dlon * constants::R_earth * cos(lat_at_curr);
        #endif 

        // Next determine how far we need to go (since we're already 
        //   at a finite distance in latitude, we don't need to go
        //   the full scale distance in longitude).
        // Essentially, use circular integration regions, not square
        //    this will further reduce the number of cells
        //    required (by up to a factor of 4), which should
        //    improve performance.
        //  The abs in local_scale is to handle the 'comfort zone'
        //    where delta_lat > scale (from the KernPad factor in dlat_N)
        #if CARTESIAN
        delta_lat   = lat_at_ilat - lat_at_curr;
        #elif not(CARTESIAN)
        delta_lat   = constants::R_earth * ( lat_at_ilat - lat_at_curr );
        #endif
        local_scale = sqrt( fabs( scale*scale - delta_lat*delta_lat ));

        // Now find the appropriate integration region
        //   The factor of 2 is diameter->radius 
        if (KernPad < 0) {
            dlon_N = Nlon;
        } else {
            dlon_N = ceil( ( KernPad * local_scale / dlon_m ) / 2.);
        }
        dlon_N = std::min(Nlon, dlon_N);
        #if PERIODIC_X
        LON_lb = Ilon - dlon_N;
        LON_ub = Ilon + dlon_N;
        if (LON_lb + Nlon < LON_ub) { LON_ub = LON_lb + Nlon; }
        #else
        LON_lb = std::max(0,    Ilon - dlon_N);
        LON_ub = std::min(Nlon, Ilon + dlon_N);
        #endif

        for (int LON = LON_lb; LON < LON_ub; LON++) {

            #if PERIODIC_X
            // Handle periodicity
            if      (LON <  0   ) { curr_lon = LON + Nlon; }
            else if (LON >= Nlon) { curr_lon = LON - Nlon; }
            else                  { curr_lon = LON; }
            #else
            curr_lon = LON;
            #endif

            index = Index(Itime, Idepth, curr_lat, curr_lon,
                          Ntime, Ndepth, Nlat,     Nlon);

            mask_index = Index(0,     0,      curr_lat, curr_lon,
                               Ntime, Ndepth, Nlat,     Nlon);

            if (distances == NULL) {
                #if CARTESIAN
                dist = distance(longitude.at(Ilon),     lat_at_ilat,
                                longitude.at(curr_lon), lat_at_curr,
                                dlon_m * Nlon, dlat_m * Nlat);
                #elif not(CARTESIAN)
                dist = distance(longitude.at(Ilon),     lat_at_ilat,
                                longitude.at(curr_lon), lat_at_curr);
                #endif
            } else {
                dist = distances->at(mask_index);
            }

            kern = kernel(dist, scale);

            area     = dAreas.at(mask_index);
            kA_sum  += kern * area;
            mask_val = mask.at(mask_index);

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

