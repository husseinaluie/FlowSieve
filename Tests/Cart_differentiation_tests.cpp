#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <algorithm>
#include <math.h>
#include <vector>
#include <assert.h>
#include "../differentiation_tools.hpp"
#include "../functions.hpp"
#include "../constants.hpp"
#include "../netcdf_io.hpp"

#ifndef SAVE_TO_FILE
    #define SAVE_TO_FILE true
#endif

double field_func(const double lat, const double lon) {
    double ret_val = cos(8 * lon + 10 * lat) * exp( - pow( lat / (M_PI / 9), 2));
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

double true_deriv_x(const double lat, const double lon) {
    double ret_val;
    if (constants::CARTESIAN) {
        ret_val = true_deriv_lon(lat, lon);
    } else {
        ret_val = 
            - true_deriv_lon(lat, lon) * sin(lon) / (cos(lat) * constants::R_earth)
            - true_deriv_lat(lat, lon) * cos(lon) *  sin(lat) / constants::R_earth;
    }
    return ret_val;
}

double true_deriv_y(const double lat, const double lon) {
    double ret_val;
    if (constants::CARTESIAN) {
        ret_val = true_deriv_lat(lat, lon);
    } else {
        ret_val = 
              true_deriv_lon(lat, lon) * cos(lon) / (cos(lat) * constants::R_earth)
            - true_deriv_lat(lat, lon) * sin(lon) *  sin(lat) / constants::R_earth;
    }
    return ret_val;
}

double true_deriv_z(const double lat, const double lon) {
    double ret_val;
    if (constants::CARTESIAN) {
        ret_val = 0.;
    } else {
        ret_val = true_deriv_lat(lat, lon) * cos(lat) / constants::R_earth;
    }
    return ret_val;
}

void apply_test(
        double & err2_x,   double & err2_y,   double & err2_z, 
        double & errinf_x, double & errinf_y, double & errinf_z,
        const std::vector<double> &longitude,
        const std::vector<double> &latitude,
        const std::vector<double> &field,
        const std::vector<double> &dArea,
        const std::vector<double> &mask,
        const int Nlat, const int Nlon) {

    // Now compute the derivatives
    std::vector<double> numer_x_deriv( Nlat * Nlon );
    std::vector<double> numer_y_deriv( Nlat * Nlon );
    std::vector<double> numer_z_deriv( Nlat * Nlon );

    #if SAVE_TO_FILE
    std::vector<double> true_x_deriv( Nlat * Nlon );
    std::vector<double> true_y_deriv( Nlat * Nlon );
    std::vector<double> true_z_deriv( Nlat * Nlon );
    #endif

    double x_deriv_val, y_deriv_val, z_deriv_val;
    std::vector<double*> x_deriv_vals, y_deriv_vals, z_deriv_vals;
    std::vector<const std::vector<double>*> deriv_fields;

    deriv_fields.push_back(&field);
    x_deriv_vals.push_back(&x_deriv_val);
    y_deriv_vals.push_back(&y_deriv_val);
    z_deriv_vals.push_back(&z_deriv_val);

    int index;
    for (int Ilat = 0; Ilat < Nlat; Ilat++) {
        for (int Ilon = 0; Ilon < Nlon; Ilon++) {
            index = Ilat * Nlon + Ilon;

            // Compute derivatives
            Cart_derivatives_at_point(
                    x_deriv_vals, y_deriv_vals,
                    z_deriv_vals, deriv_fields,
                    latitude, longitude,
                    0, 0, Ilat, Ilon,
                    1, 1, Nlat, Nlon,
                    mask);

            // Store them
            if ( mask.at(index) == 0 ) {
                numer_x_deriv.at(index) = constants::fill_value;
                numer_y_deriv.at(index) = constants::fill_value;
                numer_z_deriv.at(index) = constants::fill_value;
            } else {
                numer_x_deriv.at(index) = x_deriv_val;
                numer_y_deriv.at(index) = y_deriv_val;
                numer_z_deriv.at(index) = z_deriv_val;
            }

            #if SAVE_TO_FILE
            if ( mask.at(index) == 0 ) {
                true_x_deriv.at(index) = constants::fill_value;
                true_y_deriv.at(index) = constants::fill_value;
                true_z_deriv.at(index) = constants::fill_value;
            } else {
                true_x_deriv.at(index) = true_deriv_x(latitude.at(Ilat), longitude.at(Ilon));
                true_y_deriv.at(index) = true_deriv_y(latitude.at(Ilat), longitude.at(Ilon));
                true_z_deriv.at(index) = true_deriv_z(latitude.at(Ilat), longitude.at(Ilon));
            }
            #endif
        }
    }

    #if SAVE_TO_FILE
    char fname[50];
    snprintf(fname, 50, "test_%d_%d.nc", Nlat, Nlon);
    size_t starts[4] = {0, 0, 0, 0};
    size_t counts[4] = {1, 1, size_t(Nlat), size_t(Nlon)};

    write_field_to_output(field, "field", starts, counts, fname);

    write_field_to_output(numer_x_deriv, "numer_x_deriv", starts, counts, fname);
    write_field_to_output(numer_y_deriv, "numer_y_deriv", starts, counts, fname);
    write_field_to_output(numer_z_deriv, "numer_z_deriv", starts, counts, fname);

    write_field_to_output(true_x_deriv, "true_x_deriv", starts, counts, fname);
    write_field_to_output(true_y_deriv, "true_y_deriv", starts, counts, fname);
    write_field_to_output(true_z_deriv, "true_z_deriv", starts, counts, fname);
    #endif

    // Now measure the error
    err2_x   = 0.;
    err2_y   = 0.;
    err2_z   = 0.;
    errinf_x = 0.;
    errinf_y = 0.;
    errinf_z = 0.;
    double denom_x = 0.;
    double denom_y = 0.;
    double denom_z = 0.;
    double x_true, x_numer, y_true, y_numer, z_true, z_numer;

    for (int Ilat = 0; Ilat < Nlat; Ilat++) {
        for (int Ilon = 0; Ilon < Nlon; Ilon++) {
            index = Ilat * Nlon + Ilon;

            if (mask.at(index) == 1) {
                x_true  = true_deriv_x(latitude.at(Ilat), longitude.at(Ilon));
                x_numer = numer_x_deriv.at(index);

                y_true  = true_deriv_y(latitude.at(Ilat), longitude.at(Ilon));
                y_numer = numer_y_deriv.at(index);

                z_true  = true_deriv_z(latitude.at(Ilat), longitude.at(Ilon));
                z_numer = numer_z_deriv.at(index);

                denom_x += dArea.at(index);
                denom_y += dArea.at(index);
                denom_z += dArea.at(index);

                err2_x += dArea.at(index) * pow( x_true - x_numer, 2 );
                err2_y += dArea.at(index) * pow( y_true - y_numer, 2 );
                err2_z += dArea.at(index) * pow( z_true - z_numer, 2 );

                errinf_x = std::max(errinf_x, fabs( x_true - x_numer ) );
                errinf_y = std::max(errinf_y, fabs( y_true - y_numer ) );
                errinf_z = std::max(errinf_z, fabs( z_true - z_numer ) );
            }
        }
    }
    err2_x = sqrt( err2_x / denom_x );
    err2_y = sqrt( err2_y / denom_y );
    err2_z = sqrt( err2_z / denom_z );

}

int main(int argc, char *argv[]) {

    fprintf(stdout, "Beginning tests for differentiation routines.\n");

    int thread_safety_provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &thread_safety_provided);
    MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI::ERRORS_THROW_EXCEPTIONS);

    int wRank=-1, wSize=-1;
    MPI_Comm_rank( MPI_COMM_WORLD, &wRank );
    MPI_Comm_size( MPI_COMM_WORLD, &wSize );

    assert(wSize==1);
    assert( constants::CARTESIAN);
    assert( constants::PERIODIC_X);
    assert(!constants::PERIODIC_Y);

    const int num_tests = 6;
    const int base_nlon = 32;
    const int base_nlat = 32;

    std::vector<double> x2_errors(num_tests);
    std::vector<double> y2_errors(num_tests);
    std::vector<double> z2_errors(num_tests);
    std::vector<double> xinf_errors(num_tests);
    std::vector<double> yinf_errors(num_tests);
    std::vector<double> zinf_errors(num_tests);
    std::vector<int>    lon_points(num_tests);
    std::vector<int>    lat_points(num_tests);

    #if SAVE_TO_FILE
    std::vector<double> time{ 0. };
    std::vector<double> depth{ 0. };
    #endif

    int Nlat, Nlon;
    double lon_min, lon_max, lat_min, lat_max, dlat, dlon;
    double x2_err, y2_err, z2_err, xinf_err, yinf_err, zinf_err;

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
        int num_land = 0;

        for (int Ilat = 0; Ilat < Nlat; Ilat++) {
            for (int Ilon = 0; Ilon < Nlon; Ilon++) {
                index = Ilat * Nlon + Ilon;
                mask.at(index)  = mask_func( latitude.at(Ilat), longitude.at(Ilon));
                field.at(index) = field_func(latitude.at(Ilat), longitude.at(Ilon));

                if ( mask.at(index) == 0 ) { field.at(index) = constants::fill_value; }

                num_land += 1 - mask.at(index);
            }
        }
        if (test_ind == num_tests - 1) {
            fprintf(stdout, "At highest resolution, %.3g%% of tiles were land.\n", 
                    100 * ( (double) num_land) / (Nlat * Nlon));
        }

        #if SAVE_TO_FILE
        std::vector<std::string> vars_to_write;
        vars_to_write.push_back("field");
        vars_to_write.push_back("true_x_deriv");
        vars_to_write.push_back("true_y_deriv");
        vars_to_write.push_back("true_z_deriv");
        vars_to_write.push_back("numer_x_deriv");
        vars_to_write.push_back("numer_y_deriv");
        vars_to_write.push_back("numer_z_deriv");
        char fname[50];
        snprintf(fname, 50, "test_%d_%d.nc", Nlat, Nlon);
        initialize_output_file(time, depth, longitude, latitude, mask,
                vars_to_write, fname, 0.);
        #endif

        apply_test(x2_err, y2_err, z2_err, xinf_err, yinf_err, zinf_err, 
                longitude, latitude, field, dArea, mask, Nlat, Nlon); 

        x2_errors.at(test_ind) = x2_err;
        y2_errors.at(test_ind) = y2_err;
        z2_errors.at(test_ind) = z2_err;
        xinf_errors.at(test_ind) = xinf_err;
        yinf_errors.at(test_ind) = yinf_err;
        zinf_errors.at(test_ind) = zinf_err;

    }

    // Now that we're done applying the test, confirm that the
    //    order of convergence is as expected
     
    double cov_xy_x = 0.;
    double cov_xy_y = 0.;
    double cov_xy_z = 0.;
    double var_x_x  = 0.;
    double var_x_y  = 0.;
    double var_x_z  = 0.;

    double mean_x_x = 0.;
    double mean_x_y = 0.;
    double mean_x_z = 0.;
    double mean_y_x = 0.;
    double mean_y_y = 0.;
    double mean_y_z = 0.;

    //// First for the l2 - by lon

    // First, get means
    for (int test_ind = 0; test_ind < num_tests; test_ind++) {
        mean_x_x += log2(lon_points.at(test_ind));
        mean_x_y += log2(lon_points.at(test_ind));
        mean_x_z += log2(lon_points.at(test_ind));

        mean_y_x += log2(x2_errors.at(test_ind));
        mean_y_y += log2(y2_errors.at(test_ind));
        mean_y_z += log2(z2_errors.at(test_ind));
    }
    mean_x_x *= 1./num_tests;
    mean_x_y *= 1./num_tests;
    mean_x_z *= 1./num_tests;
    mean_y_x *= 1./num_tests;
    mean_y_y *= 1./num_tests;
    mean_y_z *= 1./num_tests;

    // Now compute variance / covariances
    for (int test_ind = 0; test_ind < num_tests; test_ind++) {
        var_x_x += pow(log2(lon_points.at(test_ind)) - mean_x_x, 2);
        var_x_y += pow(log2(lon_points.at(test_ind)) - mean_x_y, 2);
        var_x_z += pow(log2(lon_points.at(test_ind)) - mean_x_z, 2);

        cov_xy_x +=   ( log2(lon_points.at(test_ind)) - mean_x_x )
                    * ( log2( x2_errors.at(test_ind)) - mean_y_x );
        cov_xy_y +=   ( log2(lon_points.at(test_ind)) - mean_x_y )
                    * ( log2( y2_errors.at(test_ind)) - mean_y_y );
        cov_xy_z +=   ( log2(lon_points.at(test_ind)) - mean_x_z )
                    * ( log2( z2_errors.at(test_ind)) - mean_y_z );
    }

    double slope_x = cov_xy_x / var_x_x;
    double slope_y = cov_xy_y / var_x_y;
    double slope_z = cov_xy_z / var_x_z;

    fprintf(stdout, "2-norm\n");
    fprintf(stdout, "Mean convergence rates with longitude: (x, y, z) = (%.3g, %.3g, %.3g)\n", slope_x, slope_y, slope_z);
    fprintf(stdout, "  x ( Points  Error  Ord )  :  y ( Points  Error  Ord )  :  z ( Points  Error  Ord )   ( log2 of pts and err )\n");
    double ord_x=0, ord_y=0, ord_z=0;
    for (int test_ind = 0; test_ind < num_tests; test_ind++) {
        if (test_ind < num_tests - 1) {
            ord_x =   ( log2( x2_errors.at(test_ind + 1)) - log2( x2_errors.at(test_ind)) ) 
                    / ( log2(lon_points.at(test_ind + 1)) - log2(lon_points.at(test_ind)) );
            ord_y =   ( log2( y2_errors.at(test_ind + 1)) - log2( y2_errors.at(test_ind)) ) 
                    / ( log2(lon_points.at(test_ind + 1)) - log2(lon_points.at(test_ind)) );
            ord_z =   ( log2( z2_errors.at(test_ind + 1)) - log2( z2_errors.at(test_ind)) ) 
                    / ( log2(lon_points.at(test_ind + 1)) - log2(lon_points.at(test_ind)) );
        }
        fprintf(stdout, "      ( %03d    %.04g  %.3g )  :      ( %03d    %.04g  %.3g )  :      ( %03d    %.04g  %.3g )\n", 
            (int)log2(lon_points.at(test_ind)), log2(x2_errors.at(test_ind)), ord_x,
            (int)log2(lon_points.at(test_ind)), log2(y2_errors.at(test_ind)), ord_y,
            (int)log2(lon_points.at(test_ind)), log2(z2_errors.at(test_ind)), ord_z);
    }

    //// Second for the l2 - by lat

    mean_x_x = 0.;
    mean_x_y = 0.;
    mean_x_z = 0.;
    mean_y_x = 0.;
    mean_y_y = 0.;
    mean_y_z = 0.;
    // First, get means
    for (int test_ind = 0; test_ind < num_tests; test_ind++) {
        mean_x_x += log2(lat_points.at(test_ind));
        mean_x_y += log2(lat_points.at(test_ind));
        mean_x_z += log2(lat_points.at(test_ind));

        mean_y_x += log2(x2_errors.at(test_ind));
        mean_y_y += log2(y2_errors.at(test_ind));
        mean_y_z += log2(z2_errors.at(test_ind));
    }
    mean_x_x *= 1./num_tests;
    mean_x_y *= 1./num_tests;
    mean_x_z *= 1./num_tests;
    mean_y_x *= 1./num_tests;
    mean_y_y *= 1./num_tests;
    mean_y_z *= 1./num_tests;

    var_x_x = 0.;
    var_x_y = 0.;
    var_x_z = 0.;
    cov_xy_x = 0.;
    cov_xy_y = 0.;
    cov_xy_z = 0.;
    // Now compute variance / covariances
    for (int test_ind = 0; test_ind < num_tests; test_ind++) {
        var_x_x += pow(log2(lat_points.at(test_ind)) - mean_x_x, 2);
        var_x_y += pow(log2(lat_points.at(test_ind)) - mean_x_y, 2);
        var_x_z += pow(log2(lat_points.at(test_ind)) - mean_x_z, 2);

        cov_xy_x +=   ( log2(lat_points.at(test_ind)) - mean_x_x )
                    * ( log2( x2_errors.at(test_ind)) - mean_y_x );
        cov_xy_y +=   ( log2(lat_points.at(test_ind)) - mean_x_y )
                    * ( log2( y2_errors.at(test_ind)) - mean_y_y );
        cov_xy_z +=   ( log2(lat_points.at(test_ind)) - mean_x_z )
                    * ( log2( z2_errors.at(test_ind)) - mean_y_z );
    }

    slope_x = cov_xy_x / var_x_x;
    slope_y = cov_xy_y / var_x_y;
    slope_z = cov_xy_z / var_x_z;

    fprintf(stdout, "2-norm\n");
    fprintf(stdout, "Mean convergence rates with latitude: (x, y, z) = (%.3g, %.3g, %.3g)\n", slope_x, slope_y, slope_z);
    fprintf(stdout, "  x ( Points  Error  Ord )  :  y ( Points  Error  Ord )  :  z ( Points  Error  Ord )   ( log2 of pts and err )\n");
    for (int test_ind = 0; test_ind < num_tests; test_ind++) {
        if (test_ind < num_tests - 1) {
            ord_x =   ( log2( x2_errors.at(test_ind + 1)) - log2( x2_errors.at(test_ind)) ) 
                    / ( log2(lat_points.at(test_ind + 1)) - log2(lat_points.at(test_ind)) );
            ord_y =   ( log2( y2_errors.at(test_ind + 1)) - log2( y2_errors.at(test_ind)) ) 
                    / ( log2(lat_points.at(test_ind + 1)) - log2(lat_points.at(test_ind)) );
            ord_z =   ( log2( z2_errors.at(test_ind + 1)) - log2( z2_errors.at(test_ind)) ) 
                    / ( log2(lat_points.at(test_ind + 1)) - log2(lat_points.at(test_ind)) );
        }
        fprintf(stdout, "      ( %03d    %.04g  %.3g )  :      ( %03d    %.04g  %.3g )  :      ( %03d    %.04g  %.3g )\n", 
            (int)log2(lat_points.at(test_ind)), log2(x2_errors.at(test_ind)), ord_x,
            (int)log2(lat_points.at(test_ind)), log2(y2_errors.at(test_ind)), ord_y,
            (int)log2(lat_points.at(test_ind)), log2(z2_errors.at(test_ind)), ord_z);
    }

    return 0;
}
