
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <algorithm>
#include <math.h>
#include <vector>
#include "../functions.hpp"
#include "../constants.hpp"

double vel_lat(const double lat, const double lon) {
    double ret_val = cos(8 * lon + 2 * lat) * exp( - pow( lat / (M_PI / 9), 2)) * pow(lat, 2);
    return ret_val;
}

double vel_lon(const double lat, const double lon) {
    double ret_val = sin(8 * lon) * sin(2 * lat) * exp( - pow( lat / (M_PI / 9), 2));
    return ret_val;
}

double true_vort(const double lat, const double lon) {
    double ret_val = 0.;

    // d/dlon (vlat)
    ret_val += - 8 * sin(8 * lon + 2 * lat) * exp( - pow( lat / (M_PI / 9), 2)) * pow(lat, 2);

    // sin(lat) * vlon
    ret_val += sin(lat) * vel_lon(lat, lon);

    // -cos(lat) * d/dlat (vlon)
    ret_val += -cos(lat) * sin(8 * lon) * ( 2 * cos(2 * lat) ) * exp( - pow( lat / (M_PI / 9), 2));
    ret_val += -cos(lat) * sin(8 * lon) * (     sin(2 * lat) ) * exp( - pow( lat / (M_PI / 9), 2)) * ( - 2 * lat / pow(M_PI/9, 2));

    //
    ret_val *= 1./(constants::R_earth * cos(lat));

    return ret_val;
}

double mask_func(const double lat, const double lon) {
    // 1 indicates water, 0 indicates land
    double ret_val = 1.;
    
    // Make a circular island, radius pi/6
    if ( sqrt( lat*lat + lon*lon ) < M_PI/6 ) {
        ret_val *= 0.;
    }

    // Add a square island poking out in the corners
    // Essentially, just don't make the island too smooth
    if ( (abs(lat) < M_PI/7) and (abs(lon) < M_PI/7) ) {
        ret_val *= 0.;
    }

    return ret_val;
}

int main(int argc, char *argv[]) {

    fprintf(stdout, "Beginning vorticity tests.\n");

    const int Nlat = 1024;
    const int Nlon = 512;

    const double lon_min = -M_PI;
    const double lon_max =  M_PI;

    const double lat_min = -M_PI / 3;
    const double lat_max =  M_PI / 3;

    const double dlat = (lat_max - lat_min) / Nlat;
    const double dlon = (lon_max - lon_min) / Nlon;

    // Create the grid
    std::vector<double> longitude(Nlon);
    std::vector<double> latitude(Nlat);
    std::vector<double> dArea(Nlon * Nlat);
    std::vector<double> mask(Nlon * Nlat);
    std::vector<double> u_r(Nlon * Nlat);
    std::vector<double> u_lon(Nlon * Nlat);
    std::vector<double> u_lat(Nlon * Nlat);
    std::vector<double> vort_r(Nlon * Nlat);
    std::vector<double> vort_lon(Nlon * Nlat);
    std::vector<double> vort_lat(Nlon * Nlat);

    for (int II = 0; II < Nlat; II++) { latitude.at( II) = lat_min + (II+0.5) * dlat; }
    for (int II = 0; II < Nlon; II++) { longitude.at(II) = lon_min + (II+0.5) * dlon; }

    compute_areas(dArea, longitude, latitude);
    int index;

    // Initialize velocities
    for (int Ilat = 0; Ilat < Nlat; Ilat++) {
        for (int Ilon = 0; Ilon < Nlon; Ilon++) {
            index = Index(0, 0, Ilat, Ilon,
                    1, 1, Nlat, Nlon);

            mask.at(index) = mask_func(latitude.at(Ilat), longitude.at(Ilon));

            u_r.at(index) = 0.;
            u_lon.at(index) = vel_lon(latitude.at(Ilat), longitude.at(Ilon));
            u_lat.at(index) = vel_lat(latitude.at(Ilat), longitude.at(Ilon));
        }
    }

    // Compute the vorticity (numerical)
    compute_vorticity(vort_r, vort_lon, vort_lat,
            u_r, u_lon, u_lat,
            1, 1, Nlat, Nlon,
            longitude, latitude, mask);

    // Now compute the error
    double err_2     = 0.;
    double err_inf   = 0.;

    double tot_area  = 0.;

    double denom_2_t   = 0.;
    double denom_2_n   = 0.;

    double denom_inf_t = 0.;
    double denom_inf_n = 0.;

    double loc_true_vort;
    for (int Ilat = 0; Ilat < Nlat; Ilat++) {
        for (int Ilon = 0; Ilon < Nlon; Ilon++) {
            index = Index(0, 0, Ilat, Ilon,
                    1, 1, Nlat, Nlon);

            loc_true_vort = true_vort(latitude.at(Ilat), longitude.at(Ilon));

            if (mask.at(index) == 1) {
                err_2   += dArea.at(index) * pow(vort_r.at(index) - loc_true_vort, 2);
                err_inf  = std::max(err_inf, abs(vort_r.at(index) - loc_true_vort));

                denom_2_t += dArea.at(index) * pow(loc_true_vort, 2); 
                denom_2_n += dArea.at(index) * pow(vort_r.at(index), 2); 

                denom_inf_t = std::max(denom_inf_t, abs(loc_true_vort));
                denom_inf_n = std::max(denom_inf_n, abs(vort_r.at(index)));

                tot_area += dArea.at(index);
            }

        }
    }

    err_2   = sqrt(err_2) / (sqrt(denom_2_t) + sqrt(denom_2_n));
    err_inf = err_inf / (denom_inf_t + denom_inf_n);

    fprintf(stdout, "  (normalized) 2-norm   error: %.4g\n", err_2);
    fprintf(stdout, "  (normalized) inf-norm error: %.4g\n", err_inf);

    return 0;
}
