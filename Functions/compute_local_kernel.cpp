#include <math.h>
#include <vector>
#include "../functions.hpp"
#include "../constants.hpp"

/*!
 * \brief Compute an array of the kernel values from a 
 * given reference point to every other point in the domain
 *
 * (ref_ilat, ref_ilon) is the reference point from which 
 *   the kernel values are computed.
 *
 * LAT_lb and LAT_ub are the (pre-computed) latitudinal bounds for the kernel.
 *
 * @param[in,out]   local_kernel        where to store the local kernel
 * @param[in]       scale               Filtering scale
 * @param[in]       source_data         dataset class instance containing data (Psi, Phi, etc)
 * @param[in]       Ilat,Ilon           reference coordinate (kernel centre)
 * @param[in]       LAT_lb,LAT_ub       upper and lower latitudinal bounds for kernel
 *
 */
void compute_local_kernel(
        std::vector<double> & local_kernel,
        std::vector<double> & local_dl_kernel,
        std::vector<double> & local_dll_kernel,
        const double scale,
        const dataset & source_data,
        const int Ilat,
        const int Ilon,
        const int LAT_lb,
        const int LAT_ub
        ){

    const std::vector<double>   &latitude   = source_data.latitude,
                                &longitude  = source_data.longitude;

    const int   Ntime   = source_data.Ntime,
                Ndepth  = source_data.Ndepth,
                Nlat    = source_data.Nlat,
                Nlon    = source_data.Nlon;

    double dist, kern, dlat_m, dlon_m;
    size_t index;
    int Jlon, Jlat, LON_lb, LON_ub;

    const double    lat_at_ilat = latitude.at(Ilat),
                    lon_at_ilon = longitude.at(Ilon);
    double lat_at_jlat, lon_at_jlon;

    const bool do_dl  = (local_dl_kernel.size() > 0),
               do_dll = (local_dll_kernel.size() > 0);

    for (int LAT = LAT_lb; LAT < LAT_ub; LAT++) {

        // Handle periodicity
        //if (constants::PERIODIC_Y) { Jlat = ( LAT % Nlat + Nlat ) % Nlat; }
        if (constants::PERIODIC_Y) { Jlat = apply_periodicity( LAT, Nlat ); }
        else                       { Jlat = LAT; }
        #ifdef NDEBUG
        lat_at_jlat = latitude[Jlat];
        #else
        lat_at_jlat = latitude.at(Jlat);
        #endif

        // Get lon bounds at the latitude
        get_lon_bounds(LON_lb, LON_ub, longitude, Ilon, lat_at_ilat, lat_at_jlat, scale);

        for (int LON = LON_lb; LON < LON_ub; LON++) {

            // Handle periodicity
            //if (constants::PERIODIC_X) { Jlon = ( LON % Nlon + Nlon ) % Nlon; }
            if (constants::PERIODIC_X) { Jlon = apply_periodicity( LON, Nlon ); }
            else                       { Jlon = LON; }

            //index = Index(0, 0, Jlat, Jlon, Ntime, Ndepth, Nlat, Nlon);
            index = source_data.local_index(0, 0, Jlat, Jlon);

            #ifdef NDEBUG
            lon_at_jlon = longitude[Jlon];
            #else
            lon_at_jlon = longitude.at(Jlon);
            #endif

            if (constants::CARTESIAN) {
                #ifdef NDEBUG
                dlat_m = latitude[ 1] - latitude[ 0];
                dlon_m = longitude[1] - longitude[0];
                #else
                dlat_m = latitude.at( 1) - latitude.at( 0);
                dlon_m = longitude.at(1) - longitude.at(0);
                #endif
                dist = distance(lon_at_ilon, lat_at_ilat,
                                lon_at_jlon, lat_at_jlat,
                                dlon_m * Nlon, dlat_m * Nlat);
            } else {
                /*
                dist = distance(lon_at_ilon, lat_at_ilat,
                                lon_at_jlon, lat_at_jlat);
                */
                dist = source_data.distance_from_IlatIlon( Ilat, Ilon, Jlat, Jlon );
            }
            kern = kernel(dist, scale);
            #ifdef NDEBUG
            local_kernel[index] = kern;
            #else
            local_kernel.at(index) = kern;
            #endif

            // Also get the first and second ell-derivatives of the kernel
            if (do_dl) { local_dl_kernel.at(index)  = kernel(dist, scale, 1); }
            if (do_dll) { local_dll_kernel.at(index) = kernel(dist, scale, 2); }

        }
    }
}
