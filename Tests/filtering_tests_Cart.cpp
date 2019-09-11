
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <algorithm>
#include <math.h>
#include <vector>
#include <cassert>
#include "../functions.hpp"
#include "../constants.hpp"
#include "../netcdf_io.hpp"

double constant_field(const double lat, const double lon) {
    double ret_val = 100.;
    return ret_val;
}

double linear_lon_field(const double lat, const double lon) {
    double ret_val = lon;
    return ret_val;
}

double linear_lat_field(const double lat, const double lon) {
    double ret_val = lat;
    return ret_val;
}

double mask_func(const double lat, const double lon) {
    // 1 indicates water, 0 indicates land
    double ret_val = 1.;
    return ret_val;
}

int main(int argc, char *argv[]) {

    fprintf(stdout, "Beginning filtering tests.\n");

    static_assert(constants::CARTESIAN);

    int thread_safety_provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &thread_safety_provided);

    const int Nlat = 256;
    const int Nlon = 256;
    const int Npts = Nlat * Nlon;

    const double scale = 250e3;

    const double lon_min = -500e3;
    const double lon_max =  500e3;

    const double lat_min = -500e3;
    const double lat_max =  500e3;

    const double dlat = (lat_max - lat_min) / Nlat;
    const double dlon = (lon_max - lon_min) / Nlon;

    // Create the grid
    std::vector<double> time  = { 0. };
    std::vector<double> depth = { 0. };
    std::vector<double> longitude(Nlon);
    std::vector<double> latitude( Nlat);
    std::vector<double> dAreas(Npts);
    std::vector<double> mask( Npts);

    for (int II = 0; II < Nlat; II++) { latitude.at( II) = lat_min + (II+0.5) * dlat; }
    for (int II = 0; II < Nlon; II++) { longitude.at(II) = lon_min + (II+0.5) * dlon; }

    compute_areas(dAreas, longitude, latitude);
    int index, LAT_lb, LAT_ub, Ilat, Ilon;
    std::vector<double> field(Npts);
    std::vector<double> filtered(Npts);
    std::vector<double> local_kernel(Npts);
    double filter_tmp;

    // Initialize mask
    for (int Ilat = 0; Ilat < Nlat; Ilat++) {
        for (int Ilon = 0; Ilon < Nlon; Ilon++) {
            index = Index(0, 0, Ilat, Ilon, 1, 1, Nlat, Nlon);
            mask.at( index) = mask_func(latitude.at(Ilat) / lat_max, longitude.at(Ilon) / lon_max);
        }
    }

    // Prepare the output
    std::vector<std::string> vars_to_write;

    vars_to_write.push_back("kernel_sample");

    vars_to_write.push_back("constant_full");
    vars_to_write.push_back("constant_filtered");

    vars_to_write.push_back("linear_lon_full");
    vars_to_write.push_back("linear_lon_filtered");

    vars_to_write.push_back("linear_lat_full");
    vars_to_write.push_back("linear_lat_filtered");

    char fname [50];
    snprintf(fname, 50, "filtering_test_%.6gkm.nc", scale/1e3);
    initialize_output_file(time, depth, longitude, latitude, 
            mask, vars_to_write, fname, scale);

    const int ndims = 4;
    size_t starts[ndims] = { 0, 0, 0,            0 };
    size_t counts[ndims] = { 1, 1, size_t(Nlat), size_t(Nlon) };

    //
    //// Kernel sample
    //
    fprintf(stdout, "   kernel sample\n");

    Ilat = 5*Nlat / 6;
    Ilon = Nlon / 2;
    get_lat_bounds(LAT_lb, LAT_ub, latitude, Ilat, scale); 
    fprintf(stdout, "    Ilat = %d : LAT_lb, LAT_ub = %d, %d\n", Ilat, LAT_lb, LAT_ub);
    compute_local_kernel(local_kernel, scale, longitude, latitude,
            Ilat, Ilon, 1, 1, Nlat, Nlon);

    write_field_to_output(local_kernel, "kernel_sample", starts, counts, fname, &mask);

    //
    //// Constant Test
    //
    fprintf(stdout, "   filtering constants\n");

    // Initialize velocities
    for (int Ilat = 0; Ilat < Nlat; Ilat++) {
        for (int Ilon = 0; Ilon < Nlon; Ilon++) {
            index = Index(0, 0, Ilat, Ilon, 1, 1, Nlat, Nlon);
            field.at(index) = constant_field(latitude.at(Ilat) / lat_max, longitude.at(Ilon) / lon_max);
        }
    }

    for (int Ilat = 0; Ilat < Nlat; Ilat++) {
        get_lat_bounds(LAT_lb, LAT_ub, latitude, Ilat, scale); 
        for (int Ilon = 0; Ilon < Nlon; Ilon++) {
            compute_local_kernel(local_kernel, scale, longitude, latitude,
                    Ilat, Ilon, 1, 1, Nlat, Nlon);
            apply_filter_at_point( filter_tmp, field,     
                    1, 1, Nlat, Nlon, 0, 0, Ilat, Ilon,
                    longitude, latitude, LAT_lb, LAT_ub,
                    dAreas, scale, mask, true, &local_kernel );
            index = Index(0, 0, Ilat, Ilon, 1, 1, Nlat, Nlon);
            filtered.at(index) = filter_tmp;
        }
    }

    write_field_to_output(field,    "constant_full",     starts, counts, fname, &mask);
    write_field_to_output(filtered, "constant_filtered", starts, counts, fname, &mask);

    //
    //// Linear lon test
    //
    fprintf(stdout, "   filtering linear lon\n");

    // Initialize velocities
    for (int Ilat = 0; Ilat < Nlat; Ilat++) {
        for (int Ilon = 0; Ilon < Nlon; Ilon++) {
            index = Index(0, 0, Ilat, Ilon, 1, 1, Nlat, Nlon);
            field.at(index) = linear_lon_field(latitude.at(Ilat) / lat_max, longitude.at(Ilon) / lon_max);
        }
    }

    for (int Ilat = 0; Ilat < Nlat; Ilat++) {
        get_lat_bounds(LAT_lb, LAT_ub, latitude, Ilat, scale); 
        for (int Ilon = 0; Ilon < Nlon; Ilon++) {
            compute_local_kernel(local_kernel, scale, longitude, latitude,
                    Ilat, Ilon, 1, 1, Nlat, Nlon);
            apply_filter_at_point( filter_tmp, field,     
                    1, 1, Nlat, Nlon, 0, 0, Ilat, Ilon,
                    longitude, latitude, LAT_lb, LAT_ub,
                    dAreas, scale, mask, true, &local_kernel );
            index = Index(0, 0, Ilat, Ilon, 1, 1, Nlat, Nlon);
            filtered.at(index) = filter_tmp;
        }
    }

    write_field_to_output(field,    "linear_lon_full",     starts, counts, fname, &mask);
    write_field_to_output(filtered, "linear_lon_filtered", starts, counts, fname, &mask);


    //
    //// Linear lat test
    //
    fprintf(stdout, "   filtering linear lat\n");

    // Initialize velocities
    for (int Ilat = 0; Ilat < Nlat; Ilat++) {
        for (int Ilon = 0; Ilon < Nlon; Ilon++) {
            index = Index(0, 0, Ilat, Ilon, 1, 1, Nlat, Nlon);
            field.at(index) = linear_lat_field(latitude.at(Ilat) / lat_max, longitude.at(Ilon) / lon_max);
        }
    }

    for (int Ilat = 0; Ilat < Nlat; Ilat++) {
        get_lat_bounds(LAT_lb, LAT_ub, latitude, Ilat, scale); 
        for (int Ilon = 0; Ilon < Nlon; Ilon++) {
            compute_local_kernel(local_kernel, scale, longitude, latitude,
                    Ilat, Ilon, 1, 1, Nlat, Nlon);
            apply_filter_at_point( filter_tmp, field,     
                    1, 1, Nlat, Nlon, 0, 0, Ilat, Ilon,
                    longitude, latitude, LAT_lb, LAT_ub,
                    dAreas, scale, mask, true, &local_kernel );
            index = Index(0, 0, Ilat, Ilon, 1, 1, Nlat, Nlon);
            filtered.at(index) = filter_tmp;
        }
    }

    write_field_to_output(field,    "linear_lat_full",     starts, counts, fname, &mask);
    write_field_to_output(filtered, "linear_lat_filtered", starts, counts, fname, &mask);

    return 0;
}
