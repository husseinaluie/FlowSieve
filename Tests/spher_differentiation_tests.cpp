
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <algorithm>
#include <math.h>
#include <vector>
#include "../differentiation_tools.hpp"
#include "../functions.hpp"


double field_func(const double lat, const double lon) {
    double ret_val =       cos(8 * lon + 10 * lat) * exp( - pow( lat / (M_PI / 9), 2));
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
    if ( (fabs(lat) < M_PI/7) and (fabs(lon) < M_PI/7) ) {
        ret_val *= 0.;
    }

    return ret_val;
}

double true_deriv_lon(const double lat, const double lon) {
    double ret_val = - 8 * sin(8 * lon + 10 * lat) * exp( - pow( lat / (M_PI / 9), 2));
    return ret_val;
}

double true_deriv_lat(const double lat, const double lon) {
    double ret_val =        cos(8 * lon + 10 * lat) * exp( - pow( lat / (M_PI / 9), 2)) * (- 2 * lat / pow(M_PI/9, 2))
                     - 10 * sin(8 * lon + 10 * lat) * exp( - pow( lat / (M_PI / 9), 2));
    return ret_val;
}

void apply_test(double & err2_lon, double & err2_lat, 
        double & errinf_lon, double & errinf_lat,
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
            tmp = spher_derivative_at_point(field, longitude, "lon", 0, 0, Ilat, Ilon, 1, 1, Nlat, Nlon, mask);
            numer_lon_deriv.at(index) = tmp;

            // Compute latitudinal derivative
            tmp = spher_derivative_at_point(field, latitude, "lat", 0, 0, Ilat, Ilon, 1, 1, Nlat, Nlon, mask);
            numer_lat_deriv.at(index) = tmp;
        }
    }

    // Now measure the error
    err2_lon = 0.;
    err2_lat = 0.;
    errinf_lon = 0.;
    errinf_lat = 0.;
    double denom_lon = 0.;
    double denom_lat = 0.;
    double lon_true, lon_numer, lat_true, lat_numer;

    for (int Ilat = 0; Ilat < Nlat; Ilat++) {
        for (int Ilon = 0; Ilon < Nlon; Ilon++) {
            index = Ilat * Nlon + Ilon;

            if (mask.at(index) == 1) {
                lon_true  = true_deriv_lon(latitude.at(Ilat), longitude.at(Ilon));
                lon_numer = numer_lon_deriv.at(index);

                lat_true  = true_deriv_lat(latitude.at(Ilat), longitude.at(Ilon));
                lat_numer = numer_lat_deriv.at(index);

                denom_lon += dArea.at(index);
                denom_lat += dArea.at(index);

                err2_lon += dArea.at(index) * pow( lon_true - lon_numer, 2 );
                err2_lat += dArea.at(index) * pow( lat_true - lat_numer, 2 );

                errinf_lon = std::max(errinf_lon, fabs( lon_true - lon_numer ) );
                errinf_lat = std::max(errinf_lat, fabs( lat_true - lat_numer ) );
            }
        }
    }
    err2_lon = sqrt( err2_lon / denom_lon );
    err2_lat = sqrt( err2_lat / denom_lat );


}

int main(int argc, char *argv[]) {

    fprintf(stdout, "Beginning tests for differentiation routines.\n");

    const int num_tests = 10;
    const int base_nlon = 32;
    const int base_nlat = 32;

    std::vector<double> lon2_errors(num_tests);
    std::vector<double> lat2_errors(num_tests);
    std::vector<double> loninf_errors(num_tests);
    std::vector<double> latinf_errors(num_tests);
    std::vector<int>    lon_points(num_tests);
    std::vector<int>    lat_points(num_tests);

    int Nlat, Nlon;
    double lon_min, lon_max, lat_min, lat_max, dlat, dlon;
    double lon2_err, lat2_err, loninf_err, latinf_err;

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
        longitude.resize( Nlon );
        latitude.resize(  Nlat );
        dArea.resize( Nlon * Nlat );
        mask.resize(  Nlon * Nlat );
        field.resize( Nlon * Nlat );

        for (int II = 0; II < Nlat; II++) { latitude.at( II) = lat_min + (II+0.5) * dlat; }
        for (int II = 0; II < Nlon; II++) { longitude.at(II) = lon_min + (II+0.5) * dlon; }

        compute_areas(dArea, longitude, latitude);

        // Create the field to differentiate
        int index;
        int num_land = 0;

        for (int Ilat = 0; Ilat < Nlat; Ilat++) {
            for (int Ilon = 0; Ilon < Nlon; Ilon++) {
                index = Ilat * Nlon + Ilon;
                mask.at(index)  = mask_func( latitude.at(Ilat), longitude.at(Ilon));
                field.at(index) = field_func(latitude.at(Ilat), longitude.at(Ilon));

                num_land += 1 - mask.at(index);
            }
        }
        if (test_ind == num_tests - 1) {
            fprintf(stdout, "At highest resolution, %.3g%% of tiles were land.\n", 100 * ( (double) num_land) / (Nlat * Nlon));
        }

        apply_test(lon2_err, lat2_err, loninf_err, latinf_err, longitude, latitude, field, dArea, mask, Nlat, Nlon); 

        lon2_errors.at(test_ind) = lon2_err;
        lat2_errors.at(test_ind) = lat2_err;
        loninf_errors.at(test_ind) = loninf_err;
        latinf_errors.at(test_ind) = latinf_err;

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

    //// First for the l2

    // First, get means
    for (int test_ind = 0; test_ind < num_tests; test_ind++) {
        mean_x_lon += log2(lon_points.at(test_ind));
        mean_x_lat += log2(lat_points.at(test_ind));

        mean_y_lon += log2(lon2_errors.at(test_ind));
        mean_y_lat += log2(lat2_errors.at(test_ind));
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
                      * ( log2(lon2_errors.at(test_ind)) - mean_y_lon );
        cov_xy_lat +=   ( log2(lat_points.at(test_ind)) - mean_x_lat )
                      * ( log2(lat2_errors.at(test_ind)) - mean_y_lat );
    }

    double slope_lon = cov_xy_lon / var_x_lon;
    double slope_lat = cov_xy_lat / var_x_lat;

    fprintf(stdout, "2-norm\n");
    fprintf(stdout, "Mean convergence rates: (lon, lat) = (%.3g, %.3g)\n", slope_lon, slope_lat);
    fprintf(stdout, "  Lon ( Points  Error  Ord )  :  Lat ( Points  Error  Ord )   ( log2 of pts and err )\n");
    double ord_lon, ord_lat;
    for (int test_ind = 0; test_ind < num_tests; test_ind++) {
        if (test_ind < num_tests - 1) {
            ord_lon =   ( log2(lon2_errors.at(test_ind + 1)) - log2(lon2_errors.at(test_ind)) ) 
                      / ( log2(lon_points.at( test_ind + 1)) - log2(lon_points.at( test_ind)) );
            ord_lat =   ( log2(lat2_errors.at(test_ind + 1)) - log2(lat2_errors.at(test_ind)) ) 
                      / ( log2(lat_points.at( test_ind + 1)) - log2(lat_points.at( test_ind)) );
        }
        fprintf(stdout, "      ( %03d    %.04g  %.3g )  :      ( %03d    %.04g  %.3g )\n", 
            (int)log2(lon_points.at(test_ind)), log2(lon2_errors.at(test_ind)), ord_lon,
            (int)log2(lat_points.at(test_ind)), log2(lat2_errors.at(test_ind)), ord_lat);
    }

    //// Second for the linf

    mean_x_lon = 0.;
    mean_x_lat = 0.;
    mean_y_lon = 0.;
    mean_y_lat = 0.;
    var_x_lon  = 0.;
    var_x_lat  = 0.;
    cov_xy_lon = 0.;
    cov_xy_lat = 0.;

    // First, get means
    for (int test_ind = 0; test_ind < num_tests; test_ind++) {
        mean_x_lon += log2(lon_points.at(test_ind));
        mean_x_lat += log2(lat_points.at(test_ind));

        mean_y_lon += log2(loninf_errors.at(test_ind));
        mean_y_lat += log2(latinf_errors.at(test_ind));
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
                      * ( log2(loninf_errors.at(test_ind)) - mean_y_lon );
        cov_xy_lat +=   ( log2(lat_points.at(test_ind)) - mean_x_lat )
                      * ( log2(latinf_errors.at(test_ind)) - mean_y_lat );
    }

    slope_lon = cov_xy_lon / var_x_lon;
    slope_lat = cov_xy_lat / var_x_lat;

    fprintf(stdout, "inf-norm\n");
    fprintf(stdout, "Mean convergence rates: (lon, lat) = (%.3g, %.3g)\n", slope_lon, slope_lat);
    fprintf(stdout, "  Lon ( Points  Error  Ord )  :  Lat ( Points  Error  Ord )   ( log2 of pts and err )\n");
    for (int test_ind = 0; test_ind < num_tests; test_ind++) {
        if (test_ind < num_tests - 1) {
            ord_lon =   ( log2(loninf_errors.at(test_ind + 1)) - log2(loninf_errors.at(test_ind)) ) 
                      / ( log2(lon_points.at(   test_ind + 1)) - log2(lon_points.at(   test_ind)) );
            ord_lat =   ( log2(latinf_errors.at(test_ind + 1)) - log2(latinf_errors.at(test_ind)) ) 
                      / ( log2(lat_points.at(   test_ind + 1)) - log2(lat_points.at(   test_ind)) );
        }
        fprintf(stdout, "      ( %03d    %.04g  %.3g )  :      ( %03d    %.04g  %.3g )\n", 
            (int)log2(lon_points.at(test_ind)), log2(loninf_errors.at(test_ind)), ord_lon,
            (int)log2(lat_points.at(test_ind)), log2(latinf_errors.at(test_ind)), ord_lat);
    }

    return 0;
}
