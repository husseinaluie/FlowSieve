
#include <math.h>
#include <algorithm>
#include "../functions.hpp"
#include "../netcdf_io.hpp"
#include "../constants.hpp"

// If the DEBUG flag hasn't been set,
//   then use default value of 0 
#ifndef DEBUG
    #define DEBUG 0
#endif

void filtering(const double * u_r, const double * u_lon, const double * u_lat,
               const double * scales, const int Nscales,
               const double dlon, const double dlat,
               const int Ntime, const int Ndepth, const int Nlon, const int Nlat,
               const double * dAreas, 
               const double * time, const double * depth, 
               const double * longitude, const double * latitude,
               const double * mask) {

    // Now convert the Spherical velocities to Cartesian
    //   (although we will still be on a spherical
    //     coordinate system)
    #if DEBUG >= 1
    fprintf(stdout, "Converting to Cartesian velocities.\n");
    #endif
    double *u_x, *u_y, *u_z;
    u_x = new double[Ntime * Ndepth * Nlon * Nlat];
    u_y = new double[Ntime * Ndepth * Nlon * Nlat];
    u_z = new double[Ntime * Ndepth * Nlon * Nlat];

    int index, mask_index;
    //for (int Itime = 0; Itime < Ntime; Itime++) {
    for (int Itime = 0; Itime < 1; Itime++) {
        for (int Idepth = 0; Idepth < Ndepth; Idepth++) {
            for (int Ilat = 0; Ilat < Nlat; Ilat++) {
                for (int Ilon = 0; Ilon < Nlon; Ilon++) {

                    // Convert our four-index to a one-index
                    index = Index(Itime, Idepth, Ilat, Ilon,
                                  Ntime, Ndepth, Nlat, Nlon);

                    mask_index = Index(0,     0,      Ilat, Ilon,
                                       Ntime, Ndepth, Nlat, Nlon);

                    if (mask[mask_index] == 1) { // Skip land areas

                        // Note that 0 is in place of u_r, since 
                        // the current datasets don't have u_r
                        vel_Spher_to_Cart(u_x[index], u_y[index],   u_z[index],
                                          u_r[index], u_lon[index], u_lat[index],
                                          longitude[Ilon], latitude[Ilat]);
                    }
                }
            }
        }
    }

    // Create the output file
    #if DEBUG >= 1
    fprintf(stdout, "Creating output file.\n");
    #endif
    initialize_output_file(Ntime, Ndepth, Nlon, Nlat, Nscales,
            time, depth, longitude, latitude, scales, mask);

    // Now prepare to filter
    double scale,
           u_x_tmp, u_y_tmp, u_z_tmp,  // The local coarse-graining result
           u_r_tmp, u_lon_tmp, u_lat_tmp,
           dlat_m, dlon_m; // The grid spacings in metres (for lon, at a specific height)
    int dlat_N, dlon_N;    // The number of grid points required to span the filter radius

    double *coarse_u_r, *coarse_u_lon, *coarse_u_lat;
    coarse_u_r   = new double[Ntime * Ndepth * Nlat * Nlon];
    coarse_u_lon = new double[Ntime * Ndepth * Nlat * Nlon];
    coarse_u_lat = new double[Ntime * Ndepth * Nlat * Nlon];

    // The spacing (in metres) betwee latitude gridpoints
    dlat_m = dlat * constants::R_earth;

    //
    //// Begin the main filtering loop
    //
    #if DEBUG>=2
    fprintf(stdout, "Beginning main filtering loop.\n\n");
    #endif
    for (int Iscale = 0; Iscale < Nscales; Iscale++) {

        #if DEBUG >= 0
        fprintf(stdout, "Scale %d of %d\n", Iscale+1, Nscales);
        #endif

        scale  = scales[Iscale];

        // How many latitude cells are needed to span the filter radius
        dlat_N = ceil( (1.2*scale / dlat_m) / 2 );

        //for (int Itime = 0; Itime < Ntime; Itime++) {
        for (int Itime = 0; Itime < 1; Itime++) {

            #if DEBUG >= 0
            fprintf(stdout, "  Time %d of %d\n", Itime+1, Ntime);
            #endif

            for (int Idepth = 0; Idepth < Ndepth; Idepth++) {

                #if DEBUG >= 0
                fprintf(stdout, "    Depth %d of %d\n", Idepth+1, Ndepth);
                #endif

                for (int Ilat = 0; Ilat < Nlat; Ilat++) {

                    #if DEBUG >= 4
                    fprintf(stdout, "      Ilat %d of %d\n", Ilat+1, Nlat);
                    #endif

                    // Get the (maximum) number of longitude cells that 
                    //   are needed to span the filter radius
                    dlon_m = dlon 
                        * constants::R_earth 
                        * std::min( cos(latitude[Ilat] + dlat_N * dlat), 
                                    cos(latitude[Ilat] - dlat_N * dlat) ) ;
                    dlon_N = ceil( ( 1.2*scale / dlon_m) / 2 );

                    for (int Ilon = 0; Ilon < Nlon; Ilon++) {

                        #if DEBUG >= 4
                        fprintf(stdout, "        Ilon %d of %d\n", Ilon+1, Nlon);
                        #endif

                        // Convert our four-index to a one-index
                        index = Index(Itime, Idepth, Ilat, Ilon,
                                      Ntime, Ndepth, Nlat, Nlon);

                        mask_index = Index(0,     0,      Ilat, Ilon,
                                           Ntime, Ndepth, Nlat, Nlon);

                        if (mask[mask_index] == 1) { // Skip land areas

                            // Apply the filter at the point
                            apply_filter_at_point(
                                    u_x_tmp, u_y_tmp, u_z_tmp,
                                    u_x,     u_y,     u_z,
                                    dlon_N, dlat_N, 
                                    Ntime,  Ndepth, Nlat, Nlon,
                                    Itime,  Idepth, Ilat, Ilon,
                                    longitude, latitude,
                                    dAreas, scale,
                                    mask);

                            // Convert the filtered fields back to spherical
                            vel_Cart_to_Spher(u_r_tmp, u_lon_tmp, u_lat_tmp,
                                              u_x_tmp, u_y_tmp,   u_z_tmp,
                                              longitude[Ilon], latitude[Ilat]);

                            coarse_u_r[  index] = u_r_tmp;
                            coarse_u_lon[index] = u_lon_tmp;
                            coarse_u_lat[index] = u_lat_tmp;
                        }
                    }
                }
            }
        }

        // Write to file
        write_to_output(coarse_u_r, coarse_u_lon, coarse_u_lat, 
                Iscale, Ntime, Ndepth, Nlat, Nlon);

        // Now that we've filtered at the previous scale, 
        //   subtract the high-pass (large scale) from the full so
        //   that we filter the low-pass (small scale) in the next iteration
        for (int Itime = 0; Itime < Ntime; Itime++) {
            for (int Idepth = 0; Idepth < Ndepth; Idepth++) {
                for (int Ilat = 0; Ilat < Nlat; Ilat++) {
                    for (int Ilon = 0; Ilon < Nlon; Ilon++) {

                        // Convert our four-index to a one-index
                        index = Index(Itime, Idepth, Ilat, Ilon,
                                      Ntime, Ndepth, Nlat, Nlon);
                        mask_index = Index(0,     0,      Ilat, Ilon,
                                           Ntime, Ndepth, Nlat, Nlon);

                        if (mask[mask_index] == 1) { // Skip land areas

                            vel_Spher_to_Cart(u_x_tmp,           u_y_tmp,             u_z_tmp,
                                              coarse_u_r[index], coarse_u_lon[index], coarse_u_lat[index],
                                              longitude[Ilon],   latitude[Ilat]);

                            u_x[index] = u_x[index] - u_x_tmp;
                            u_y[index] = u_y[index] - u_y_tmp;
                            u_z[index] = u_z[index] - u_z_tmp;
                        }
                    }
                }
            }
        }

    } // end of scale loop

    // Now write the remaining small scales
    for (int Itime = 0; Itime < Ntime; Itime++) {
        for (int Idepth = 0; Idepth < Ndepth; Idepth++) {
            for (int Ilat = 0; Ilat < Nlat; Ilat++) {
                for (int Ilon = 0; Ilon < Nlon; Ilon++) {

                    // Convert our four-index to a one-index
                    index = Index(Itime, Idepth, Ilat, Ilon,
                                  Ntime, Ndepth, Nlat, Nlon);
                    mask_index = Index(0,     0,      Ilat, Ilon,
                                       Ntime, Ndepth, Nlat, Nlon);

                    if (mask[mask_index] == 1) { // Skip land areas

                        vel_Cart_to_Spher(
                                u_r_tmp,    u_lon_tmp,  u_lat_tmp,
                                u_x[index], u_y[index], u_z[index],
                                longitude[Ilon],   latitude[Ilat]);

                        coarse_u_r[  index] = u_r_tmp;
                        coarse_u_lon[index] = u_lon_tmp;
                        coarse_u_lat[index] = u_lat_tmp;
                    }
                }
            }
        }
    }

    write_to_output(coarse_u_r, coarse_u_lon, coarse_u_lat, 
            Nscales, Ntime, Ndepth, Nlat, Nlon);

    // Free up the arrays
    delete[] coarse_u_r;
    delete[] coarse_u_lon;
    delete[] coarse_u_lat;
    
    delete[] u_x;
    delete[] u_y;
    delete[] u_z;

}
