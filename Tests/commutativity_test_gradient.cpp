
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <algorithm>
#include <math.h>
#include <vector>
#include <omp.h>
#include <cassert>
#include "../functions.hpp"
#include "../constants.hpp"
#include "../netcdf_io.hpp"
#include "../differentiation_tools.hpp"

const double D2R = M_PI / 180;

double full_field(const double lat, const double lon) {
    double ret_val = 0;
    double coef1, coef2;

    coef1 = 1.5;
    coef2 = 2.0;

    const double lon0 =  0 * D2R;
    const double lon1 = 60 * D2R;
    const double lat0 = 30 * D2R;

    // Gaussian bump centred on (lon0, lat0)
    const double width_rad = 20 * D2R;
    ret_val += coef1 * exp( - pow( (lon - lon0) / width_rad, 2)
                            - pow( (lat - lat0) / width_rad, 2) );

    // Gaussian ring centred on (lon0, lat0) passing through (lon1, lat0)
    const double d0 = distance(lon1, lat0, lon0, lat0);
    const double d  = distance(lon, lat, lon0, lat0);
    const double width_m = distance(lon1, lat0, lon1 + 6*D2R, lat0);
    ret_val += coef2 * exp( -pow( (d-d0) / width_m, 2) );

    // Add some sin/cos waviness around
    ret_val += sin( 16 * lon + 12 * lat) * cos( 10 * lon - 8 * lat );

    return ret_val;
}

double mask_func(const double lat, const double lon) {
    // 1 indicates water, 0 indicates land
    double ret_val = 1.;

    // For now, no land, just water
    
    return ret_val;
}

int main(int argc, char *argv[]) {

    fprintf(stdout, "Beginning commutativity test.\n");

    int thread_safety_provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &thread_safety_provided);

    static_assert( constants::PERIODIC_X);
    static_assert(!constants::PERIODIC_Y);

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
        scale = M_PI / 10.;
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

    std::vector<double>  field(     Npts);

    for (int II = 0; II < Nlat; II++) { latitude.at( II) = lat_min + (II+0.5) * dlat; }
    for (int II = 0; II < Nlon; II++) { longitude.at(II) = lon_min + (II+0.5) * dlon; }

    compute_areas(dArea, longitude, latitude);
    int Ilat, Ilon, index, perc_base, perc, perc_count;

    // Initialize field
    for (Ilat = 0; Ilat < Nlat; Ilat++) {
        for (Ilon = 0; Ilon < Nlon; Ilon++) {
            index = Index(0, 0, Ilat, Ilon, 1, 1, Nlat, Nlon);

            mask.at( index) = mask_func( latitude.at(Ilat), longitude.at(Ilon));
            field.at(index) = full_field(latitude.at(Ilat), longitude.at(Ilon));
        }
    }

    // Differentiate the field
    fprintf(stdout, "   differentiation the field\n");
    double ddlon_field, ddlat_field;
    std::vector<double> gradf_lon(Npts);
    std::vector<double> gradf_lat(Npts);

    double ddx_field, ddy_field, ddz_field;
    std::vector<double> gradf_x(Npts);
    std::vector<double> gradf_y(Npts);
    std::vector<double> gradf_z(Npts);

    std::vector<const std::vector<double>*> deriv_fields {&field};

    std::vector<double*> lon_deriv_vals {&ddlon_field};
    std::vector<double*> lat_deriv_vals {&ddlat_field};

    std::vector<double*> x_deriv_vals {&ddx_field};
    std::vector<double*> y_deriv_vals {&ddy_field};
    std::vector<double*> z_deriv_vals {&ddz_field};

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
                    lat_deriv_vals, deriv_fields, latitude, "lat",
                    Itime, Idepth, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon,
                    mask);

            spher_derivative_at_point(
                    lon_deriv_vals, deriv_fields, longitude, "lon",
                    Itime, Idepth, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon,
                    mask);

            Cart_derivatives_at_point(
                    x_deriv_vals, y_deriv_vals, z_deriv_vals, deriv_fields,
                    latitude, longitude,
                    Itime, Idepth, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon,
                    mask);

            gradf_lon.at(index) = ddlon_field / cos(latitude.at(Ilat));
            gradf_lat.at(index) = ddlat_field;

            gradf_x.at(index) = ddx_field;
            gradf_y.at(index) = ddy_field;
            gradf_z.at(index) = ddz_field;
        }
    }
    fprintf(stdout, "\n");

    std::vector<double> local_kernel(Npts);

    std::vector<double> coarse(      Npts);
    std::vector<double> coarse_ddlon(Npts);
    std::vector<double> coarse_ddlat(Npts);

    std::vector<double> coarse_ddx(Npts);
    std::vector<double> coarse_ddy(Npts);
    std::vector<double> coarse_ddz(Npts);

    int LAT_lb, LAT_ub;

    fprintf(stdout, "   filtering\n");

    perc_base = 5;
    perc = 0; 
    perc_count=0;
    fprintf(stdout, "      ");

    // Set up filtering vectors
    double filt_field, filt_dfdlon, filt_dfdlat, filt_dfdx, filt_dfdy, filt_dfdz;
    std::vector<double*> filtered_vals;
    std::vector<const std::vector<double>*> filter_fields;

    filter_fields.push_back(&field);
    filter_fields.push_back(&gradf_lon);
    filter_fields.push_back(&gradf_lat);
    filter_fields.push_back(&gradf_x);
    filter_fields.push_back(&gradf_y);
    filter_fields.push_back(&gradf_z);

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
        shared(Ilat, field, gradf_lon, gradf_lat,\
                gradf_x, gradf_y, gradf_z,\
                coarse, coarse_ddlon, coarse_ddlat,\
                coarse_ddx, coarse_ddy, coarse_ddz,\
                LAT_lb, LAT_ub, dArea, mask, scale,\
                latitude, longitude, filter_fields) \
        private(index, Ilon, filtered_vals,\
                filt_field, filt_dfdlon, filt_dfdlat,\
                filt_dfdx, filt_dfdy, filt_dfdz)\
        firstprivate(local_kernel)
        {

            filtered_vals.resize(filter_fields.size());
            filtered_vals.at(0) = &filt_field;
            filtered_vals.at(1) = &filt_dfdlon;
            filtered_vals.at(2) = &filt_dfdlat;
            filtered_vals.at(3) = &filt_dfdx;
            filtered_vals.at(4) = &filt_dfdy;
            filtered_vals.at(5) = &filt_dfdz;

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

                coarse.at(index) = filt_field;

                coarse_ddlon.at(index) = filt_dfdlon;
                coarse_ddlat.at(index) = filt_dfdlat;

                coarse_ddx.at(index) = filt_dfdx;
                coarse_ddy.at(index) = filt_dfdy;
                coarse_ddz.at(index) = filt_dfdz;
            }
        }
    }
    fprintf(stdout, "\n");

    // Now differentiate the coarse field
    deriv_fields.resize(1);
    deriv_fields.at(0) = &coarse;

    std::vector<double> ddlon_coarse(Npts);
    std::vector<double> ddlat_coarse(Npts);

    std::vector<double> ddx_coarse(Npts);
    std::vector<double> ddy_coarse(Npts);
    std::vector<double> ddz_coarse(Npts);

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
                    lat_deriv_vals, deriv_fields, latitude, "lat",
                    Itime, Idepth, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon,
                    mask);

            spher_derivative_at_point(
                    lon_deriv_vals, deriv_fields, longitude, "lon",
                    Itime, Idepth, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon,
                    mask);

            Cart_derivatives_at_point(
                    x_deriv_vals, y_deriv_vals, z_deriv_vals, deriv_fields,
                    latitude, longitude,
                    Itime, Idepth, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon,
                    mask);

            ddlon_coarse.at(index) = ddlon_field / cos(latitude.at(Ilat));
            ddlat_coarse.at(index) = ddlat_field;

            ddx_coarse.at(index) = ddx_field;
            ddy_coarse.at(index) = ddy_field;
            ddz_coarse.at(index) = ddz_field;

        }
    }
    fprintf(stdout, "\n");

    fprintf(stdout, "\nPreparing to output\n");

    // Write the fields to file for interest
    char fname[50];
    snprintf(fname, 50, "commutativity_test_gradient_%d_%d.nc", Nlat, Nlon);
    size_t starts[4] = {0, 0, 0, 0};
    size_t counts[4] = {1, 1, size_t(Nlat), size_t(Nlon)};

    std::vector<std::string> vars_to_write;
    vars_to_write.push_back("field");
    vars_to_write.push_back("coarse");

    vars_to_write.push_back("coarse_ddlon");
    vars_to_write.push_back("coarse_ddlat");
    vars_to_write.push_back("coarse_ddx"  );
    vars_to_write.push_back("coarse_ddy"  );
    vars_to_write.push_back("coarse_ddz"  );

    vars_to_write.push_back("ddlon_coarse");
    vars_to_write.push_back("ddlat_coarse");
    vars_to_write.push_back("ddx_coarse"  );
    vars_to_write.push_back("ddy_coarse"  );
    vars_to_write.push_back("ddz_coarse"  );

    fprintf(stdout, "  Initializing the file\n");
    initialize_output_file(times, depth, longitude, latitude, mask,
            vars_to_write, fname, scale);

    fprintf(stdout, "  Writing reference fields\n");
    write_field_to_output(field,  "field",  starts, counts, fname);
    write_field_to_output(coarse, "coarse", starts, counts, fname);

    fprintf(stdout, "  Writing filtered derivatives\n");
    write_field_to_output(coarse_ddlon, "coarse_ddlon", starts, counts, fname);
    write_field_to_output(coarse_ddlat, "coarse_ddlat", starts, counts, fname);
    write_field_to_output(coarse_ddx,   "coarse_ddx",   starts, counts, fname);
    write_field_to_output(coarse_ddy,   "coarse_ddy",   starts, counts, fname);
    write_field_to_output(coarse_ddz,   "coarse_ddz",   starts, counts, fname);

    write_field_to_output(ddlon_coarse, "ddlon_coarse", starts, counts, fname);
    write_field_to_output(ddlat_coarse, "ddlat_coarse", starts, counts, fname);
    write_field_to_output(ddx_coarse,   "ddx_coarse",   starts, counts, fname);
    write_field_to_output(ddy_coarse,   "ddy_coarse",   starts, counts, fname);
    write_field_to_output(ddz_coarse,   "ddz_coarse",   starts, counts, fname);

    return 0;
}
