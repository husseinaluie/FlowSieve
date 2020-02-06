
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <algorithm>
#include <math.h>
#include <vector>
#include "../functions.hpp"
#include "../constants.hpp"
#include "../netcdf_io.hpp"

const double L_env = M_PI / 9;
// Irrotational component
//    f = ( sin(lon) * cos(lat)  +  0.5 * cos(10 * lon) * sin(9 * lat) ) 
//        * np.exp( - (lat / L_env)^2 )
const double irrot_coef = 1.;



const double toroi_coef = 0.;



double u_r_field(const double lat, const double lon) {
    double ret_val = 0.;
    return ret_val;
}

double u_lon_field(const double lat, const double lon) {
    double ret_val = 0.;
    double env = exp( - pow(lat / L_env, 2) );

    // Irrotational component
    //    ddlon(f) / cos(lat)
    ret_val += irrot_coef 
        * ( cos(lon) * cos(lat)  -  0.5 * 10 * sin(10*lon) * sin(9*lat) ) 
        * env / cos(lat);

    // Toroidal component
    //    

    return ret_val;
}

double u_lat_field(const double lat, const double lon) {
    double ret_val = 0.;
    double env = exp( - pow(lat / L_env, 2) );

    // Irrotational component
    //    ddlat(f)
    ret_val += irrot_coef
        * ( (-2*lat/pow(L_env,2)) * ( sin(lon) * cos(lat)  +  0.5 *   cos(10*lon) * sin(9*lat) ) 
            +                       (-sin(lon) * sin(lat)  +  0.5 * 9*cos(10*lon) * cos(9*lat) ) )
        * env;

    return ret_val;
}

double mask_func(const double lat, const double lon) {
    // 1 indicates water, 0 indicates land
    double ret_val = 1.;
    return ret_val;
}

int main(int argc, char *argv[]) {

    fprintf(stdout, "Beginning filtering tests.\n");

    int thread_safety_provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &thread_safety_provided);

    const int Nlat = 256;
    const int Nlon = Nlat * 2;
    const int Npts = Nlat * Nlon;

    const double scale = 2000e3;

    const double lon_min = -M_PI;
    const double lon_max =  M_PI;

    const double lat_min = -M_PI / 2;
    const double lat_max =  M_PI / 2;

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

    std::vector<double> full_u_r(  Npts);
    std::vector<double> full_u_lon(Npts);
    std::vector<double> full_u_lat(Npts);

    std::vector<double> filt_u_r(  Npts);
    std::vector<double> filt_u_lon(Npts);
    std::vector<double> filt_u_lat(Npts);

    std::vector<double> full_u_x(Npts);
    std::vector<double> full_u_y(Npts);
    std::vector<double> full_u_z(Npts);

    std::vector<double> filt_u_x(Npts);
    std::vector<double> filt_u_y(Npts);
    std::vector<double> filt_u_z(Npts);

    double filt_u_tmp, filt_v_tmp, filt_w_tmp, u_x, u_y, u_z;

    // Initialize mask and velocities
    for (int Ilat = 0; Ilat < Nlat; Ilat++) {
        for (int Ilon = 0; Ilon < Nlon; Ilon++) {
            index = Index(0, 0, Ilat, Ilon, 1, 1, Nlat, Nlon);
            mask.at( index) = mask_func(latitude.at(Ilat), longitude.at(Ilon));

            full_u_r.at(  index) = u_r_field(  latitude.at(Ilat), longitude.at(Ilon));
            full_u_lon.at(index) = u_lon_field(latitude.at(Ilat), longitude.at(Ilon));
            full_u_lat.at(index) = u_lat_field(latitude.at(Ilat), longitude.at(Ilon));

            vel_Spher_to_Cart(
                    u_x, u_y, u_z,
                    full_u_r.at(index), full_u_lon.at(index), full_u_lat.at(index),
                    longitude.at(Ilon), latitude.at(Ilat)
                    );

            full_u_x.at(index) = u_x;
            full_u_y.at(index) = u_y;
            full_u_z.at(index) = u_z;
        }
    }

    // Prepare the output
    std::vector<std::string> vars_to_write;

    vars_to_write.push_back("kernel_sample");
    vars_to_write.push_back("kernel_int_region");

    vars_to_write.push_back("u_r_full");
    vars_to_write.push_back("u_r_filtered");

    vars_to_write.push_back("u_lon_full");
    vars_to_write.push_back("u_lon_filtered");

    vars_to_write.push_back("u_lat_full");
    vars_to_write.push_back("u_lat_filtered");

    vars_to_write.push_back("u_x_full");
    vars_to_write.push_back("u_x_filtered");

    vars_to_write.push_back("u_y_full");
    vars_to_write.push_back("u_y_filtered");

    vars_to_write.push_back("u_z_full");
    vars_to_write.push_back("u_z_filtered");

    char fname [50];
    snprintf(fname, 50, "velocity_filtering_test_%.6gkm.nc", scale/1e3);
    initialize_output_file(time, depth, longitude, latitude, 
            mask, vars_to_write, fname, scale);

    const int ndims = 4;
    size_t starts[ndims] = { 0, 0, 0,            0 };
    size_t counts[ndims] = { 1, 1, size_t(Nlat), size_t(Nlon) };

    //
    //// Kernel sample
    //
    fprintf(stdout, "   kernel sample\n");

    //Ilat = 7*Nlat / 8;
    double arc_lat = constants::R_earth * dlat;
    int len_arc = ceil(scale / arc_lat);
    Ilat = Nlat - 1 * len_arc / 3;
    Ilon = 0;
    get_lat_bounds(LAT_lb, LAT_ub, latitude, Ilat, scale); 
    fprintf(stdout, "    Ilat = %d : LAT_lb, LAT_ub = %d, %d\n", Ilat, LAT_lb, LAT_ub);

    std::vector<double> local_kernel(Npts);
    std::vector<double> kernel_int_region(Npts);

    compute_local_kernel(local_kernel, scale, longitude, latitude,
            Ilat, Ilon, 1, 1, Nlat, Nlon);

    // Create a variable that indicates the region of integration for filtering
    int LON_lb, LON_ub, LON, curr_lon;
    const double centre_lat = latitude.at(Ilat);
    for (Ilat = LAT_lb; Ilat < LAT_ub; Ilat++) {

        get_lon_bounds(LON_lb, LON_ub, longitude, Ilon, 
                centre_lat, latitude.at(Ilat), scale);

        for (LON = LON_lb; LON < LON_ub; LON++) {

            // Handle periodicity
            if (constants::PERIODIC_X) {
                if      (LON <  0   ) { curr_lon = LON + Nlon; }
                else if (LON >= Nlon) { curr_lon = LON - Nlon; }
                else                  { curr_lon = LON; }
            } else {
                curr_lon = LON;
            }

            index = Index(0, 0, Ilat, curr_lon, 1, 1, Nlat, Nlon);
            kernel_int_region.at(index) = 1;
        }
    }

    write_field_to_output(local_kernel, "kernel_sample", starts, counts, fname, &mask);
    write_field_to_output(kernel_int_region, "kernel_int_region", starts, counts, fname, &mask);

    //
    //// Constant Test
    //
    fprintf(stdout, "   beginning filtering\n");

    std::vector<double*> filtered_vals = {
        &filt_u_tmp, &filt_v_tmp, &filt_w_tmp
    };
    std::vector<const std::vector<double>*> filter_fields = {
        &full_u_x, &full_u_y, &full_u_z
    };

    int perc_base = 5;
    int perc = 0; 
    int perc_count=0;
    fprintf(stdout, "      ");

    for (int Ilat = 0; Ilat < Nlat; Ilat++) {
        get_lat_bounds(LAT_lb, LAT_ub, latitude, Ilat, scale); 

        // Every perc_base percent, print a dot, but only the first thread
        if ( ((double)(Ilat) / Nlat) * 100 >= perc ) {
            perc_count++;
            if (perc_count % 5 == 0) { fprintf(stdout, "|"); }
            else                     { fprintf(stdout, "."); }
            fflush(stdout);
            perc += perc_base;
        }
        for (int Ilon = 0; Ilon < Nlon; Ilon++) {

            compute_local_kernel(local_kernel, scale, longitude, latitude,
                    Ilat, Ilon, 1, 1, Nlat, Nlon);

            apply_filter_at_point( 
                    filtered_vals, filter_fields,     
                    1, 1, Nlat, Nlon, 0, 0, Ilat, Ilon,
                    longitude, latitude, LAT_lb, LAT_ub,
                    dAreas, scale, mask, true, &local_kernel );

            index = Index(0, 0, Ilat, Ilon, 1, 1, Nlat, Nlon);

            filt_u_x.at(index) = filt_u_tmp;
            filt_u_y.at(index) = filt_v_tmp;
            filt_u_z.at(index) = filt_w_tmp;

            vel_Cart_to_Spher(
                    u_z, u_x, u_y,
                    filt_u_tmp, filt_v_tmp, filt_w_tmp,
                    longitude.at(Ilon), latitude.at(Ilat)
                    );

            filt_u_r.at(  index) = u_z;
            filt_u_lon.at(index) = u_x;
            filt_u_lat.at(index) = u_y;
        }
    }

    write_field_to_output(full_u_r,   "u_r_full",   starts, counts, fname, &mask);
    write_field_to_output(full_u_lon, "u_lon_full", starts, counts, fname, &mask);
    write_field_to_output(full_u_lat, "u_lat_full", starts, counts, fname, &mask);

    write_field_to_output(filt_u_r,   "u_r_filtered",   starts, counts, fname, &mask);
    write_field_to_output(filt_u_lon, "u_lon_filtered", starts, counts, fname, &mask);
    write_field_to_output(filt_u_lat, "u_lat_filtered", starts, counts, fname, &mask);

    write_field_to_output(full_u_x, "u_x_full", starts, counts, fname, &mask);
    write_field_to_output(full_u_y, "u_y_full", starts, counts, fname, &mask);
    write_field_to_output(full_u_z, "u_z_full", starts, counts, fname, &mask);

    write_field_to_output(filt_u_x, "u_x_filtered", starts, counts, fname, &mask);
    write_field_to_output(filt_u_y, "u_y_filtered", starts, counts, fname, &mask);
    write_field_to_output(filt_u_z, "u_z_filtered", starts, counts, fname, &mask);

    return 0;
}
