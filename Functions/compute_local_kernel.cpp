#include <math.h>
#include <vector>
#include "../functions.hpp"
#include "../constants.hpp"

void compute_local_kernel(
        std::vector<double> & local_kernel,
        const double scale,
        const std::vector<double> & longitude,
        const std::vector<double> & latitude,
        const int Ilat,
        const int Ilon,
        const int Ntime,
        const int Ndepth,
        const int Nlat,
        const int Nlon,
        const int LAT_lb,
        const int LAT_ub
        ){

    double dist, kern, dlat_m, dlon_m;
    size_t index;
    int curr_lon, curr_lat, LON_lb, LON_ub;

    const double lat_at_ilat = latitude.at(Ilat);
    const double lon_at_ilon = longitude.at(Ilon);
    double lat_at_curr;

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

        // Get lon bounds at the latitude
        get_lon_bounds(LON_lb, LON_ub, longitude, Ilon, lat_at_ilat, lat_at_curr, scale);

        for (int LON = LON_lb; LON < LON_ub; LON++) {

            // Handle periodicity
            if (constants::PERIODIC_X) {
                if      (LON <  0   ) { curr_lon = LON + Nlon; }
                else if (LON >= Nlon) { curr_lon = LON - Nlon; }
                else                  { curr_lon = LON; }
            } else {
                curr_lon = LON;
            }

            index = Index(0,     0,      curr_lat, curr_lon,
                          Ntime, Ndepth, Nlat,     Nlon);

            if (constants::CARTESIAN) {
                dlat_m = latitude.at( 1) - latitude.at( 0);
                dlon_m = longitude.at(1) - longitude.at(0);
                dist = distance(lon_at_ilon,     lat_at_ilat,
                                longitude.at(curr_lon), lat_at_curr,
                                dlon_m * Nlon, dlat_m * Nlat);
            } else {
                dist = distance(lon_at_ilon,            lat_at_ilat,
                                longitude.at(curr_lon), lat_at_curr);
            }
            kern = kernel(dist, scale);

            local_kernel.at(index) = kern;

        }
    }
}
