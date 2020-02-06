/*
 *
 * This test provides a 2D divergence-free velocity field.
 * 
 * Two scalar fields are then computed:
 *   - the divergence of the filtered fields
 *   - the filtered divergence of the fields
 *
 * The error between these two fields should be small, indicating
 *   that the divergence and filtering operators do indeed commute.
 *
 * In spherical coordinates, if we assume that there is no radial vel.,
 *   then the divergence reduces to:
 *   div(v) = ( ddlon(v_lon) + ddlat( v_lat * cos(lat) )  ) / (r * cos(lat) )
 *
 * Of course, in Cartesian coordinates, the divergence of horizontal vel.
 *   is simply given by:
 *   div(v) = ddx(v_x) + ddy(v_y)
 *
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <algorithm>
#include <math.h>
#include <vector>
#include <omp.h>
#include "../functions.hpp"
#include "../constants.hpp"
#include "../netcdf_io.hpp"
#include "../differentiation_tools.hpp"

const double D2R = M_PI / 180;

const int N_fine1  = 9;
const int N_fine2  = 12;
const int N_coarse = 5;
const double L_env = M_PI / 9;
const double coef1 = 0.5;

double full_lon_vel(const double lat, const double lon) {
    double ret_val = 0;
    double env = exp( - pow(lat / L_env, 2) );

    if (constants::CARTESIAN) {
        ret_val = 
              -env
             * ( 
                 - (2 * lat / pow(L_env, 2) ) 
                   * (
                        -cos(N_coarse * lon) / N_coarse
                        + (coef1 / N_fine1) * sin(N_fine1 * lon) * sin(N_fine2 * lat)
                     )
                 + (coef1 * N_fine2 / N_fine1) * sin(N_fine1 * lon) * cos(N_fine2 * lat)
               );
    } else {
        ret_val = 
              - ( cos(N_coarse * lon) / N_coarse ) 
            * ( sin(lat) + 2 * lat * cos(lat) / pow(L_env, 2)  )
            * env;
        ret_val += 
              coef1 
            * ( sin(N_fine1 * lon) / N_fine1 )
            * (   sin(lat) * sin(N_fine2 * lat)
                + cos(lat) * (
                    - sin(N_fine2 * lat) * 2 * lat / pow(L_env, 2)
                    + cos(N_fine2 * lat) * N_fine2
                    )
                )
            * env;
    }

    return ret_val;
}

double full_lat_vel(const double lat, const double lon) {
    double ret_val = 0;
    double env = exp( - pow(lat / L_env, 2) );

    if (constants::CARTESIAN) {
        ret_val  =         sin(N_coarse * lon) * env;
        ret_val += coef1 * cos(N_fine1 * lon) * sin(N_fine2 * lat) * env;
    } else {
        ret_val  =         sin(N_coarse * lon) * env;
        ret_val += coef1 * cos(N_fine1 * lon) * sin(N_fine2 * lat) * env;
    }

    return ret_val;
}

double mask_func(const double lat, const double lon) {
    // 1 indicates water, 0 indicates land
    double ret_val = 1.;

    // For now, no land, just water
    return ret_val;
}

int main(int argc, char *argv[]) {

    fprintf(stdout, "Beginning DIV commutativity test.\n");

    int thread_safety_provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &thread_safety_provided);

    const int Ntime  = 1;
    const int Ndepth = 1;
    const int Nlat   = 256;
    const int Nlon   = 512;
    const int Npts   = Ntime * Ndepth * Nlat * Nlon;

    const int Itime  = 0;
    const int Idepth = 0;

    const double lon_min = -M_PI;
    const double lon_max =  M_PI;

    const double lat_min = -M_PI / 2;
    const double lat_max =  M_PI / 2;

    double scale;
    if (constants::CARTESIAN) {
        scale = M_PI / 14.;
    } else {
        scale = 2500e3;
    }
    const double dlat = (lat_max - lat_min) / Nlat;
    const double dlon = (lon_max - lon_min) / Nlon;

    // Create the grid
    std::vector<double> times = { 0. };
    std::vector<double> depth = { 0. };
    std::vector<double> longitude(Nlon);
    std::vector<double> latitude( Nlat);
    std::vector<double> dArea( Npts);
    std::vector<double> mask(  Npts);

    std::vector<double> u_lon(Npts);
    std::vector<double> u_lat(Npts);

    double u_tmp, v_tmp, w_tmp;
    std::vector<double> u_x(Npts);
    std::vector<double> u_y(Npts);
    std::vector<double> u_z(Npts);

    for (int II = 0; II < Nlat; II++) { latitude.at( II) = lat_min + (II+0.5) * dlat; }
    for (int II = 0; II < Nlon; II++) { longitude.at(II) = lon_min + (II+0.5) * dlon; }

    compute_areas(dArea, longitude, latitude);
    int Ilat, Ilon, index, perc_base, perc, perc_count;

    // Initialize field
    for (Ilat = 0; Ilat < Nlat; Ilat++) {
        for (Ilon = 0; Ilon < Nlon; Ilon++) {
            index = Index(0, 0, Ilat, Ilon, 1, 1, Nlat, Nlon);

            mask.at( index) = mask_func( latitude.at(Ilat), longitude.at(Ilon));

            if (mask.at(index) == 1) {
                u_lon.at(index) = full_lon_vel(latitude.at(Ilat), longitude.at(Ilon));
                u_lat.at(index) = full_lat_vel(latitude.at(Ilat), longitude.at(Ilon));

                vel_Spher_to_Cart( 
                        u_tmp, v_tmp, w_tmp, 
                        0, u_lon.at(index), u_lat.at(index),
                        longitude.at(Ilon), latitude.at(Ilat)  
                        );

                u_x.at(index) = u_tmp;
                u_y.at(index) = v_tmp;
                u_z.at(index) = w_tmp;
            } else {
                u_lon.at(index) = constants::fill_value;
                u_lat.at(index) = constants::fill_value;

                u_x.at(index) = constants::fill_value;
                u_y.at(index) = constants::fill_value;
                u_z.at(index) = constants::fill_value;
            }
        }
    }

    // Differentiate the field
    fprintf(stdout, "   computing divergence of full velocity field\n");
    fprintf(stdout, "       (both Cartesian and spherical)\n");
    double ddlon_ulon, ddlat_ulat, ddx_ux, ddy_uy, ddz_uz;
    std::vector<double> div_cart_vel(Npts);
    std::vector<double> div_sphe_vel(Npts);

    std::vector<const std::vector<double>*> sphe_deriv_fields {&u_lon, &u_lat};
    std::vector<const std::vector<double>*> cart_deriv_fields {&u_x, &u_y, &u_z};

    std::vector<double*> lon_deriv_vals {&ddlon_ulon, NULL};
    std::vector<double*> lat_deriv_vals {NULL, &ddlat_ulat};

    std::vector<double*> x_deriv_vals {&ddx_ux, NULL, NULL};
    std::vector<double*> y_deriv_vals {NULL, &ddy_uy, NULL};
    std::vector<double*> z_deriv_vals {NULL, NULL, &ddz_uz};

    perc_base = 5;
    perc = 0; 
    perc_count=0;
    fprintf(stdout, "      ");

    double loc_lat;

    for (Ilat = 0; Ilat < Nlat; Ilat++) {

        // Every perc_base percent, print a dot, but only the first thread
        if ( ((double)(Ilat) / Nlat) * 100 >= perc ) {
            perc_count++;
            if (perc_count % 5 == 0) { fprintf(stdout, "|"); }
            else                     { fprintf(stdout, "."); }
            fflush(stdout);
            perc += perc_base;
        }

        for (Ilon = 0; Ilon < Nlon; Ilon++) {
            index = Index(0, 0, Ilat, Ilon, 1, 1, Nlat, Nlon);

            spher_derivative_at_point(
                    lat_deriv_vals, sphe_deriv_fields, latitude, "lat",
                    Itime, Idepth, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon,
                    mask);

            spher_derivative_at_point(
                    lon_deriv_vals, sphe_deriv_fields, longitude, "lon",
                    Itime, Idepth, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon,
                    mask);

            Cart_derivatives_at_point(
                    x_deriv_vals, y_deriv_vals, z_deriv_vals, cart_deriv_fields,
                    latitude, longitude,
                    Itime, Idepth, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon,
                    mask);

            div_cart_vel.at(index) = ddx_ux + ddy_uy + ddz_uz; 

            if (constants::CARTESIAN) {
                div_sphe_vel.at(index) = 0.;
            } else {
                loc_lat = latitude.at(Ilat);
                div_sphe_vel.at(index) = 
                    (ddlon_ulon + ddlat_ulat * cos(loc_lat) - u_lat.at(index) * sin(loc_lat))
                    / (constants::R_earth * cos(loc_lat));
            }

        }
    }
    fprintf(stdout, "\n");

    std::vector<double> local_kernel(Npts);

    std::vector<double> coarse_ux(Npts);
    std::vector<double> coarse_uy(Npts);
    std::vector<double> coarse_uz(Npts);

    std::vector<double> coarse_urad(Npts);
    std::vector<double> coarse_ulon(Npts);
    std::vector<double> coarse_ulat(Npts);

    std::vector<double> coarse_div_cart(Npts);
    std::vector<double> coarse_div_sphe(Npts);

    int LAT_lb, LAT_ub;

    fprintf(stdout, "   filtering\n");

    perc_base = 5;
    perc = 0; 
    perc_count=0;
    fprintf(stdout, "      ");

    // Set up filtering vectors
    double filt_ux, filt_uy, filt_uz, filt_div_cart, filt_div_sphe;
    std::vector<double*> filtered_vals;
    std::vector<const std::vector<double>*> filter_fields = {
        &u_x, &u_y, &u_z, &div_cart_vel, &div_sphe_vel
    };

    // Filter full and derivates of full
    for (Ilat = 0; Ilat < Nlat; Ilat++) {

        get_lat_bounds(LAT_lb, LAT_ub, latitude, Ilat, scale); 

        // Every perc_base percent, print a dot, but only the first thread
        if ( ((double)(Ilat) / Nlat) * 100 >= perc ) {
            perc_count++;
            if (perc_count % 5 == 0) { fprintf(stdout, "|"); }
            else                     { fprintf(stdout, "."); }
            fflush(stdout);
            perc += perc_base;
        }

        #pragma omp parallel \
        default(none) \
        shared(Ilat, latitude, longitude, filter_fields,\
                coarse_ux, coarse_uy, coarse_uz,\
                coarse_urad, coarse_ulon, coarse_ulat,\
                coarse_div_cart, coarse_div_sphe,\
                LAT_lb, LAT_ub, dArea, mask, scale)\
        private(index, Ilon, filtered_vals,\
                u_tmp, v_tmp, w_tmp,\
                filt_ux, filt_uy, filt_uz,\
                filt_div_cart, filt_div_sphe)\
        firstprivate(local_kernel)
        {

            filtered_vals.resize(filter_fields.size());
            filtered_vals.at(0) = &filt_ux;
            filtered_vals.at(1) = &filt_uy;
            filtered_vals.at(2) = &filt_uz;
            filtered_vals.at(3) = &filt_div_cart;
            filtered_vals.at(4) = &filt_div_sphe;

            #pragma omp for collapse(1) schedule(static)
            for (Ilon = 0; Ilon < Nlon; Ilon++) {

                index = Index(0, 0, Ilat, Ilon, 1, 1, Nlat, Nlon);

                compute_local_kernel(
                        local_kernel, scale, longitude, latitude,
                        Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon);
                
                // Filter coarse field
                apply_filter_at_point(
                        filtered_vals, filter_fields,
                        Ntime, Ndepth, Nlat, Nlon,
                        Itime, Idepth, Ilat, Ilon,
                        longitude, latitude,
                        LAT_lb, LAT_ub,
                        dArea, scale, mask, true,
                        &local_kernel);

                coarse_ux.at(index) = filt_ux;
                coarse_uy.at(index) = filt_uy;
                coarse_uz.at(index) = filt_uz;

                vel_Cart_to_Spher(
                        w_tmp, u_tmp, v_tmp,
                        filt_ux, filt_uy, filt_uz,
                        longitude.at(Ilon), latitude.at(Ilat)                    
                        );

                coarse_urad.at(index) = w_tmp;
                coarse_ulon.at(index) = u_tmp;
                coarse_ulat.at(index) = v_tmp;

                coarse_div_cart.at(index) = filt_div_cart;
                coarse_div_sphe.at(index) = filt_div_sphe;
            }
        }
    }
    fprintf(stdout, "\n");

    // Now differentiate the coarse field
    std::vector<double> div_coarse_cart_vel(Npts);
    std::vector<double> div_coarse_sphe_vel(Npts);

    sphe_deriv_fields.at(0) = &coarse_ulon;
    sphe_deriv_fields.at(1) = &coarse_ulat;

    cart_deriv_fields.at(0) = &coarse_ux;
    cart_deriv_fields.at(1) = &coarse_uy;
    cart_deriv_fields.at(2) = &coarse_uz;

    fprintf(stdout, "   differentiating filtered field\n");

    perc_base = 5;
    perc = 0; 
    perc_count=0;
    fprintf(stdout, "      ");

    for (Ilat = 0; Ilat < Nlat; Ilat++) {

        // Every perc_base percent, print a dot, but only the first thread
        if ( ((double)(Ilat) / Nlat) * 100 >= perc ) {
            perc_count++;
            if (perc_count % 5 == 0) { fprintf(stdout, "|"); }
            else                     { fprintf(stdout, "."); }
            fflush(stdout);
            perc += perc_base;
        }

        for (Ilon = 0; Ilon < Nlon; Ilon++) {
            index = Index(0, 0, Ilat, Ilon, 1, 1, Nlat, Nlon);

            spher_derivative_at_point(
                    lat_deriv_vals, sphe_deriv_fields, latitude, "lat",
                    Itime, Idepth, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon,
                    mask);

            spher_derivative_at_point(
                    lon_deriv_vals, sphe_deriv_fields, longitude, "lon",
                    Itime, Idepth, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon,
                    mask);

            Cart_derivatives_at_point(
                    x_deriv_vals, y_deriv_vals, z_deriv_vals, cart_deriv_fields,
                    latitude, longitude,
                    Itime, Idepth, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon,
                    mask);

            div_coarse_cart_vel.at(index) = ddx_ux + ddy_uy + ddz_uz; 

            if (constants::CARTESIAN) {
                div_coarse_sphe_vel.at(index) = 0.;
            } else {
                loc_lat = latitude.at(Ilat);
                div_coarse_sphe_vel.at(index) = 
                    (ddlon_ulon + ddlat_ulat * cos(loc_lat) - u_lat.at(index) * sin(loc_lat))
                    / (constants::R_earth * cos(loc_lat));
            }

        }
    }
    fprintf(stdout, "\n");

    fprintf(stdout, "\nPreparing to output\n");

    // Write the fields to file for interest
    char fname[50];
    snprintf(fname, 50, "commutativity_test_divergence_%d_%d.nc", Nlat, Nlon);
    size_t starts[4] = {0, 0, 0, 0};
    size_t counts[4] = {1, 1, size_t(Nlat), size_t(Nlon)};

    std::vector<std::string> vars_to_write;
    vars_to_write.push_back("u_lon");
    vars_to_write.push_back("u_lat");
    vars_to_write.push_back("u_x");
    vars_to_write.push_back("u_y");
    vars_to_write.push_back("u_z");

    vars_to_write.push_back("div_cart_vel");
    vars_to_write.push_back("div_sphe_vel");

    vars_to_write.push_back("coarse_ux");
    vars_to_write.push_back("coarse_uy");
    vars_to_write.push_back("coarse_uz");

    vars_to_write.push_back("coarse_urad");
    vars_to_write.push_back("coarse_ulon");
    vars_to_write.push_back("coarse_ulat");

    vars_to_write.push_back("coarse_div_cart");
    vars_to_write.push_back("coarse_div_sphe");

    vars_to_write.push_back("div_coarse_cart_vel");
    vars_to_write.push_back("div_coarse_sphe_vel");

    fprintf(stdout, "  Initializing the file\n");
    initialize_output_file(times, depth, longitude, latitude, mask,
            vars_to_write, fname, scale);

    write_field_to_output(u_lon, "u_lon", starts, counts, fname);
    write_field_to_output(u_lat, "u_lat", starts, counts, fname);
    write_field_to_output(u_x,   "u_x",   starts, counts, fname);
    write_field_to_output(u_y,   "u_y",   starts, counts, fname);
    write_field_to_output(u_z,   "u_z",   starts, counts, fname);

    write_field_to_output(coarse_urad, "coarse_urad", starts, counts, fname);
    write_field_to_output(coarse_ulon, "coarse_ulon", starts, counts, fname);
    write_field_to_output(coarse_ulat, "coarse_ulat", starts, counts, fname);
    write_field_to_output(coarse_ux,   "coarse_ux",   starts, counts, fname);
    write_field_to_output(coarse_uy,   "coarse_uy",   starts, counts, fname);
    write_field_to_output(coarse_uz,   "coarse_uz",   starts, counts, fname);

    write_field_to_output(div_cart_vel, "div_cart_vel", starts, counts, fname);
    write_field_to_output(div_sphe_vel, "div_sphe_vel", starts, counts, fname);

    write_field_to_output(coarse_div_cart, "coarse_div_cart", starts, counts, fname);
    write_field_to_output(coarse_div_sphe, "coarse_div_sphe", starts, counts, fname);

    write_field_to_output(div_coarse_cart_vel, "div_coarse_cart_vel", starts, counts, fname);
    write_field_to_output(div_coarse_sphe_vel, "div_coarse_sphe_vel", starts, counts, fname);

    return 0;
}
