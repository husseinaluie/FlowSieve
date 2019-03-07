
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <algorithm>
#include <math.h>
#include <vector>
#include "../differentiation_tools.hpp"
#include "../functions.hpp"


double field_func(const double lat, const double lon) {
    double ret_val =       cos(8 * lon + 2 * lat) * exp( - pow( lat / (M_PI / 9), 2));
    return ret_val;
}

double true_deriv_lon(const double lat, const double lon) {
    double ret_val = - 8 * sin(8 * lon + 2 * lat) * exp( - pow( lat / (M_PI / 9), 2));
    return ret_val;
}

double true_deriv_lat(const double lat, const double lon) {
    double ret_val =       cos(8 * lon + 2 * lat) * exp( - pow( lat / (M_PI / 9), 2)) * (- 2 * lat / pow(M_PI/9, 2))
                     - 2 * sin(8 * lon + 2 * lat) * exp( - pow( lat / (M_PI / 9), 2));
    return ret_val;
}

void apply_test(double & err_lon, double & err_lat, 
        const std::vector<double> &longitude,
        const std::vector<double> &latitude,
        const std::vector<double> &field,
        const std::vector<double> &dArea,
        const std::vector<double> &mask,
        const int Nlat, const int Nlon) {

    // Now compute the derivatives
    std::vector<double> numer_lon_deriv( Nlat * Nlon );
    std::vector<double> numer_lat_deriv( Nlat * Nlon );
    double tmp;
    int index;

    for (int Ilat = 0; Ilat < Nlat; Ilat++) {
        for (int Ilon = 0; Ilon < Nlon; Ilon++) {
            index = Ilat * Nlon + Ilon;

            // Compute longitudinal derivative
            tmp = longitude_derivative_at_point(field, longitude, 0, 0, Ilat, Ilon, 1, 1, Nlat, Nlon, mask);
            numer_lon_deriv.at(index) = tmp;

            // Compute latitudinal derivative
            tmp = latitude_derivative_at_point(field, latitude, 0, 0, Ilat, Ilon, 1, 1, Nlat, Nlon, mask);
            numer_lat_deriv.at(index) = tmp;
        }
    }

    // Now measure the error
    err_lon = 0.;
    err_lat = 0.;
    double denom_lon = 0.;
    double denom_lat = 0.;
    double lon_true, lon_numer, lat_true, lat_numer;

    for (int Ilat = 0; Ilat < Nlat; Ilat++) {
        for (int Ilon = 0; Ilon < Nlon; Ilon++) {
            index = Ilat * Nlon + Ilon;

            lon_true  = true_deriv_lon(latitude.at(Ilat), longitude.at(Ilon));
            lon_numer = numer_lon_deriv.at(index);

            lat_true  = true_deriv_lat(latitude.at(Ilat), longitude.at(Ilon));
            lat_numer = numer_lat_deriv.at(index);

            denom_lon += dArea.at(index);
            denom_lat += dArea.at(index);

            err_lon += dArea.at(index) * pow( lon_true - lon_numer, 2 );
            err_lat += dArea.at(index) * pow( lat_true - lat_numer, 2 );
        }
    }
    err_lon *= 1./denom_lon;
    err_lat *= 1./denom_lat;

}

int main(int argc, char *argv[]) {

    fprintf(stdout, "Beginning tests for differentiation routines.\n");

    const int num_tests = 8;
    const int base_nlon = 32;
    const int base_nlat = 32;

    std::vector<double> lon_errors(num_tests);
    std::vector<double> lat_errors(num_tests);
    std::vector<int>    lon_points(num_tests);
    std::vector<int>    lat_points(num_tests);

    int Nlat, Nlon;
    double lon_min, lon_max, lat_min, lat_max, dlat, dlon;
    double lon_err, lat_err;

    std::vector<double> longitude, latitude, dArea, 
        mask, field;

    for (int test_ind = 0; test_ind < num_tests; test_ind++) {

        // Specify grid information
        Nlat = base_nlat * pow(2, test_ind);
        Nlon = base_nlon * pow(2, test_ind);

        lon_points.at(test_ind) = Nlon;
        lat_points.at(test_ind) = Nlat;

        lon_min = -M_PI;
        lon_max =  M_PI;

        lat_min = -M_PI / 3;
        lat_max =  M_PI / 3;

        dlat = (lat_max - lat_min) / Nlat;
        dlon = (lon_max - lon_min) / Nlon;

        // Create the grid
        longitude.resize(Nlon);
        latitude.resize(Nlat);
        dArea.resize(Nlon * Nlat);
        mask.resize(Nlon * Nlat);
        field.resize(Nlon * Nlat);

        for (int II = 0; II < Nlat; II++) { latitude.at( II) = lat_min + (II+0.5) * dlat; }
        for (int II = 0; II < Nlon; II++) { longitude.at(II) = lon_min + (II+0.5) * dlon; }

        compute_areas(dArea, longitude, latitude);

        // Create the field to differentiate
        int index;

        // cos(long) * exp( - lat**2 )
        for (int Ilat = 0; Ilat < Nlat; Ilat++) {
            for (int Ilon = 0; Ilon < Nlon; Ilon++) {
                index = Ilat * Nlon + Ilon;
                mask.at(index) = 1.;
                field.at(index) = field_func(latitude.at(Ilat), longitude.at(Ilon));
            }
        }

        apply_test(lon_err, lat_err, longitude, latitude, field, dArea, mask, Nlat, Nlon); 

        lon_errors.at(test_ind) = lon_err;
        lat_errors.at(test_ind) = lat_err;

    }

    // Now that we're done applying the test, confirm that the
    //    order of convergence is as expected
     
    double cov_xy_lon = 0.;
    double cov_xy_lat = 0.;
    double var_x_lon  = 0.;
    double var_x_lat  = 0.;

    double mean_x_lon = 0.;
    double mean_x_lat = 0.;
    double mean_y_lon = 0.;
    double mean_y_lat = 0.;

    // First, get means
    for (int test_ind = 0; test_ind < num_tests; test_ind++) {
        mean_x_lon += log2(lon_points.at(test_ind));
        mean_x_lat += log2(lat_points.at(test_ind));

        mean_y_lon += log2(lon_errors.at(test_ind));
        mean_y_lat += log2(lat_errors.at(test_ind));
    }
    mean_x_lon *= 1./num_tests;
    mean_x_lat *= 1./num_tests;
    mean_y_lon *= 1./num_tests;
    mean_y_lat *= 1./num_tests;

    // Now compute variance / covariances
    for (int test_ind = 0; test_ind < num_tests; test_ind++) {
        var_x_lon += pow(log2(lon_points.at(test_ind)) - mean_x_lon, 2);
        var_x_lat += pow(log2(lat_points.at(test_ind)) - mean_x_lat, 2);

        cov_xy_lon +=   ( log2(lon_points.at(test_ind)) - mean_x_lon )
                      * ( log2(lon_errors.at(test_ind)) - mean_y_lon );
        cov_xy_lat +=   ( log2(lat_points.at(test_ind)) - mean_x_lat )
                      * ( log2(lat_errors.at(test_ind)) - mean_y_lat );
    }

    double slope_lon = cov_xy_lon / var_x_lon;
    double slope_lat = cov_xy_lat / var_x_lat;

    fprintf(stdout, "Mean convergence rates: (lon, lat) = (%.3g, %.3g)\n", slope_lon, slope_lat);
    fprintf(stdout, "  Lon ( Points  Error )  :  Lat ( Points  Error )   (log2)\n");
    for (int test_ind = 0; test_ind < num_tests; test_ind++) {
        fprintf(stdout, "      ( %03d    %.04g )  :      ( %03d    %.04g )\n", 
            (int)log2(lon_points.at(test_ind)), log2(lon_errors.at(test_ind)),
            (int)log2(lat_points.at(test_ind)), log2(lat_errors.at(test_ind)));
    }

    return 0;
}
