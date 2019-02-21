
#include <math.h>
#include <algorithm>
#include "../functions.hpp"
#include "../constants.hpp"

#ifndef DEBUG
    #define DEBUG false
#endif

void filtering(const double * u_r, const double * u_lon, const double * u_lat,
               const double * scales, const int Nscales,
               const double dlon, const double dlat,
               const int Ntime, const int Ndepth, const int Nlon, const int Nlat,
               const double * dAreas, const double * longitude, const double * latitude) {

    const bool debug = DEBUG;

    // Now convert the Spherical velocities to Cartesian
    //   (although we will still be on a spherical
    //     coordinate system)
    if (debug) {
        fprintf(stdout, "Converting to Cartesian velocities.\n");
    }
    double *u_x, *u_y, *u_z;
    u_x = new double[Ntime * Ndepth * Nlon * Nlat];
    u_y = new double[Ntime * Ndepth * Nlon * Nlat];
    u_z = new double[Ntime * Ndepth * Nlon * Nlat];

    int index;
    for (int Itime = 0; Itime < Ntime; Itime++) {
        for (int Idepth = 0; Idepth < Ndepth; Idepth++) {
            for (int Ilat = 0; Ilat < Nlat; Ilat++) {
                for (int Ilon = 0; Ilon < Nlon; Ilon++) {

                    // Convert our four-index to a one-index
                    index = Index(Itime, Idepth, Ilat, Ilon,
                                  Ntime, Ndepth, Nlat, Nlon);

                    // Note that 0 is in place of u_r, since 
                    // the current datasets don't have u_r
                    vel_Spher_to_Cart(u_x[index], u_y[index],   u_z[index],
                                      0.,         u_lon[index], u_lat[index],
                                      longitude[Ilon], latitude[Ilat]);
                }
            }
        }
    }

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
    for (int Iscale = 0; Iscale < Nscales; Iscale++) {

        fprintf(stdout, "Scale %d of %d\n", Iscale+1, Nscales);

        scale  = scales[Iscale];

        // How many latitude cells are needed to span the filter radius
        dlat_N = ceil( (scale / dlat_m) / 2 );

        for (int Itime = 0; Itime < Ntime; Itime++) {

            fprintf(stdout, "Time %d of %d\n", Itime+1, Ntime);

            for (int Idepth = 0; Idepth < Ndepth; Idepth++) {

                fprintf(stdout, "Depth %d of %d\n", Idepth+1, Ndepth);

                for (int Ilat = 0; Ilat < Nlat; Ilat++) {

                    for (int Ilon = 0; Ilon < Nlon; Ilon++) {

                        // Convert our four-index to a one-index
                        index = Index(Itime, Idepth, Ilat, Ilon,
                                      Ntime, Ndepth, Nlat, Nlon);

                        // Get the (maximum) number of longitude cells that 
                        //   are needed to span the filter radius
                        dlon_m =   dlon 
                                 * constants::R_earth 
                                 * std::min( cos(latitude[Ilat] + dlat_N * dlat), 
                                             cos(latitude[Ilat] - dlat_N * dlat) ) ;
                        dlon_N = ceil( (scale / dlon_m) / 2 );

                        // Apply the filter at the point
                        apply_filter_at_point(
                                u_x_tmp, u_y_tmp, u_z_tmp,
                                u_x,     u_y,     u_z,
                                dlon_N, dlat_N, 
                                Ntime,  Ndepth, Nlat, Nlon,
                                Itime,  Idepth, Ilat, Ilon,
                                longitude, latitude,
                                dAreas, scale);

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

        // Write to file
        //
        // < code to write this filter to output file >

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

                        vel_Spher_to_Cart(u_x_tmp, u_y_tmp,   u_z_tmp,
                                          u_r_tmp, u_lon_tmp, u_lat_tmp,
                                          longitude[Ilon], latitude[Ilat]);

                        u_x[index] = u_x[index] - u_x_tmp;
                        u_y[index] = u_y[index] - u_y_tmp;
                        u_z[index] = u_z[index] - u_z_tmp;
                    }
                }
            }
        }

    } // end of scale loop

    // Free up the arrays
    delete[] coarse_u_r;
    delete[] coarse_u_lon;
    delete[] coarse_u_lat;
    
    delete[] u_x;
    delete[] u_y;
    delete[] u_z;

}
