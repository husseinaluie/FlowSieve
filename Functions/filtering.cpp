#include <math.h>
#include <algorithm>
#include <vector>
#include "../functions.hpp"
#include "../netcdf_io.hpp"
#include "../constants.hpp"

#ifndef DEBUG
    #define DEBUG 0
#endif

#ifndef COMP_VORT
    #define COMP_VORT true
#endif

#ifndef COMP_TRANSFERS
    #define COMP_TRANSFERS true
#endif

void filtering(
        const std::vector<double> & full_u_r,       /**< [in] Full u_r velocity array */
        const std::vector<double> & full_u_lon,     /**< [in] Full u_lon velocity array */
        const std::vector<double> & full_u_lat,     /**< [in] Full u_lat velocity array */
        const std::vector<double> & scales,    /**< [in] Array of filtering scales */
        const std::vector<double> & dAreas,    /**< [in] Array of cell areas (2D) (compute_areas()) */
        const std::vector<double> & time,      /**< [in] Time dimension (1D) */
        const std::vector<double> & depth,     /**< [in] Depth dimension (1D) */
        const std::vector<double> & longitude, /**< [in] Longitude dimension (1D) */
        const std::vector<double> & latitude,  /**< [in] Latitude dimension (1D) */
        const std::vector<double> & mask       /**< [in] Array to distinguish between land and water cells (2D) */
        ) {


    // Compute grid spacing
    //   for the moment, assume a uniform grid
    double dlat = latitude.at( 1) - latitude.at( 0);
    double dlon = longitude.at(1) - longitude.at(0);

    // Get dimension sizes
    const int Nscales = scales.size();
    const int Ntime   = time.size();
    const int Ndepth  = depth.size();
    const int Nlat    = latitude.size();
    const int Nlon    = longitude.size();

    const int num_pts = Ntime * Ndepth * Nlat * Nlon;

    #if DEBUG >= 1
    fprintf(stdout, "Converting to Cartesian velocities.\n");
    #endif

    std::vector<double> u_x(num_pts);
    std::vector<double> u_y(num_pts);
    std::vector<double> u_z(num_pts);

    std::vector<double> coarse_u_r(  num_pts);
    std::vector<double> coarse_u_lon(num_pts);
    std::vector<double> coarse_u_lat(num_pts);

    // Now convert the Spherical velocities to Cartesian
    //   (although we will still be on a spherical
    //     coordinate system)
    int index, mask_index;
    for (int Itime = 0; Itime < Ntime; Itime++) {
        for (int Idepth = 0; Idepth < Ndepth; Idepth++) {
            for (int Ilat = 0; Ilat < Nlat; Ilat++) {
                for (int Ilon = 0; Ilon < Nlon; Ilon++) {

                    // Convert our four-index to a one-index
                    index = Index(Itime, Idepth, Ilat, Ilon,
                                  Ntime, Ndepth, Nlat, Nlon);

                    mask_index = Index(0,     0,      Ilat, Ilon,
                                       Ntime, Ndepth, Nlat, Nlon);

                    if (mask.at(mask_index) == 1) { // Skip land areas

                        // Note that 0 is in place of u_r, since 
                        // the current datasets don't have u_r
                        vel_Spher_to_Cart(     u_x.at(index),      u_y.at(  index),      u_z.at(  index),
                                          full_u_r.at(index), full_u_lon.at(index), full_u_lat.at(index),
                                          longitude.at(Ilon), latitude.at(Ilat));
                    }

                    // Copy the current full u_r into the initial coarse
                    coarse_u_r.at(  index) = full_u_r.at(  index);
                    coarse_u_lon.at(index) = full_u_lon.at(index);
                    coarse_u_lat.at(index) = full_u_lat.at(index);
                }
            }
        }
    }

    // Create the output file
    #if DEBUG >= 1
    fprintf(stdout, "Creating output file.\n");
    #endif
    initialize_output_file(time, depth, longitude, latitude, scales, mask);

    // Now prepare to filter
    double scale,
           u_x_tmp, u_y_tmp,   u_z_tmp,  // The local coarse-graining result
           u_r_tmp, u_lon_tmp, u_lat_tmp,
           dlat_m, dlon_m; // The grid spacings in metres (for lon, at a specific height)
    int dlat_N, dlon_N;    // The number of grid points required to span the filter radius

    std::vector<double> fine_u_r(  num_pts);
    std::vector<double> fine_u_lon(num_pts);
    std::vector<double> fine_u_lat(num_pts);

    #if COMP_VORT
    std::vector<double> vort_r(  num_pts);
    std::vector<double> vort_lat(num_pts);
    std::vector<double> vort_lon(num_pts);
    #endif

    #if COMP_TRANSFERS
    // If computing energy transfers, we'll need some more arrays
    std::vector<double> coarse_uxux(num_pts);
    std::vector<double> coarse_uxuy(num_pts);
    std::vector<double> coarse_uxuz(num_pts);
    std::vector<double> coarse_uyuy(num_pts);
    std::vector<double> coarse_uyuz(num_pts);
    std::vector<double> coarse_uzuz(num_pts);

    double uxux_tmp, uxuy_tmp, uxuz_tmp, uyuy_tmp, uyuz_tmp, uzuz_tmp;

    // We'll also need to keep the coarse velocities
    std::vector<double> coarse_u_x(num_pts);
    std::vector<double> coarse_u_y(num_pts);
    std::vector<double> coarse_u_z(num_pts);

    // Also an array for the transfer itself
    std::vector<double> energy_transfer(num_pts);

    // If we're computing transfers, then we already have what
    //   we need to computed band-filtered KE, so might as well do it
    std::vector<double> fine_KE(num_pts);
    #endif

    // The spacing (in metres) betwee latitude gridpoints
    dlat_m = dlat * constants::R_earth;

    #if DEBUG >= 2
    int perc;
    #endif

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

        scale  = scales.at(Iscale);

        // How many latitude cells are needed to span the filter radius
        dlat_N = ceil( (1.1*scale / dlat_m) / 2 );

        //for (int Itime = 0; Itime < Ntime; Itime++) {
        for (int Itime = 0; Itime < 1; Itime++) {

            #if DEBUG >= 0
            fprintf(stdout, "  Time %d of %d\n", Itime+1, Ntime);
            #endif

            for (int Idepth = 0; Idepth < Ndepth; Idepth++) {

                #if DEBUG >= 0
                fprintf(stdout, "    Depth %d of %d\n", Idepth+1, Ndepth);
                #endif

                #if DEBUG >= 2
                perc = 0;
                #endif

                for (int Ilat = 0; Ilat < Nlat; Ilat++) {

                    #if DEBUG >= 2
                    // Every 10 percent, print a dot
                    if ( ( ((double) Ilat) / Nlat) * 100 > perc ) {
                        fprintf(stdout, ".");
                        perc += 10;
                    }
                    #endif

                    #if DEBUG >= 3
                    fprintf(stdout, "      Ilat %d of %d\n", Ilat+1, Nlat);
                    #endif

                    // Get the (maximum) number of longitude cells that 
                    //   are needed to span the filter radius
                    dlon_m = dlon 
                        * constants::R_earth 
                        * std::min( cos(latitude.at(Ilat) + dlat_N * dlat), 
                                    cos(latitude.at(Ilat) - dlat_N * dlat) ) ;
                    dlon_N = ceil( ( 1.1*scale / dlon_m) / 2 );

                    for (int Ilon = 0; Ilon < Nlon; Ilon++) {

                        #if DEBUG >= 3
                        fprintf(stdout, "        Ilon %d of %d\n", Ilon+1, Nlon);
                        #endif

                        // Convert our four-index to a one-index
                        index = Index(Itime, Idepth, Ilat, Ilon,
                                      Ntime, Ndepth, Nlat, Nlon);

                        mask_index = Index(0,     0,      Ilat, Ilon,
                                           Ntime, Ndepth, Nlat, Nlon);

                        if (mask.at(mask_index) == 1) { // Skip land areas

                            // Apply the filter at the point
                            #if DEBUG >= 4
                            fprintf(stdout, "          Line %d of %s\n", __LINE__, __FILE__);
                            #endif
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
                            #if DEBUG >= 4
                            fprintf(stdout, "          Line %d of %s\n", __LINE__, __FILE__);
                            #endif
                            vel_Cart_to_Spher(u_r_tmp, u_lon_tmp, u_lat_tmp,
                                              u_x_tmp, u_y_tmp,   u_z_tmp,
                                              longitude.at(Ilon), latitude.at(Ilat));

                            // Subtract current coarse from preceeding coarse to
                            //    get current fine
                            #if DEBUG >= 4
                            fprintf(stdout, "          Line %d of %s\n", __LINE__, __FILE__);
                            #endif
                            fine_u_r.at(  index) = coarse_u_r.at(  index) - u_r_tmp;
                            fine_u_lon.at(index) = coarse_u_lon.at(index) - u_lon_tmp;
                            fine_u_lat.at(index) = coarse_u_lat.at(index) - u_lat_tmp;

                            // Now pass the new coarse along as the preceeding coarse.
                            #if DEBUG >= 4
                            fprintf(stdout, "          Line %d of %s\n", __LINE__, __FILE__);
                            #endif
                            coarse_u_r.at(  index) = u_r_tmp;
                            coarse_u_lon.at(index) = u_lon_tmp;
                            coarse_u_lat.at(index) = u_lat_tmp;

                            #if COMP_TRANSFERS
                            #if DEBUG >= 4
                            fprintf(stdout, "          Line %d of %s\n", __LINE__, __FILE__);
                            #endif
                            apply_filter_at_point_for_quadratics(
                                    uxux_tmp, uxuy_tmp, uxuz_tmp,
                                    uyuy_tmp, uyuz_tmp, uzuz_tmp,
                                    u_x,      u_y,      u_z,
                                    dlon_N, dlat_N, 
                                    Ntime,  Ndepth, Nlat, Nlon,
                                    Itime,  Idepth, Ilat, Ilon,
                                    longitude, latitude,
                                    dAreas, scale,
                                    mask);

                            coarse_uxux.at(index) = uxux_tmp;
                            coarse_uxuy.at(index) = uxuy_tmp;
                            coarse_uxuz.at(index) = uxuz_tmp;
                            coarse_uyuy.at(index) = uyuy_tmp;
                            coarse_uyuz.at(index) = uyuz_tmp;
                            coarse_uzuz.at(index) = uzuz_tmp;

                            coarse_u_x.at(index) = u_x_tmp;
                            coarse_u_y.at(index) = u_y_tmp;
                            coarse_u_z.at(index) = u_z_tmp;
                            
                            fine_KE.at(index) = 0.5 * (   
                                      ( uxux_tmp - u_x_tmp * u_x_tmp )
                                    + ( uyuy_tmp - u_y_tmp * u_y_tmp )
                                    + ( uzuz_tmp - u_z_tmp * u_z_tmp )
                                    );
                            #endif
                        }
                    }
                }
                #if DEBUG >= 2
                fprintf(stdout, "\n");
                #endif
            }
        }

        // Write to file
        write_to_output(fine_u_r, fine_u_lon, fine_u_lat, 
                Iscale, Ntime, Ndepth, Nlat, Nlon);

        #if COMP_VORT
        // Compute and write vorticity
        compute_vorticity(vort_r, vort_lon, vort_lat,
                fine_u_r, fine_u_lon, fine_u_lat,
                Ntime, Ndepth, Nlat, Nlon,
                longitude, latitude, mask);
        write_vorticity(vort_r, vort_lon, vort_lat,
                Iscale, Ntime, Ndepth, Nlat, Nlon);
        #endif

        #if COMP_TRANSFERS
        // Compute the energy transfer through the filter scale
        compute_energy_transfer_through_scale(
                energy_transfer, 
                coarse_u_x,  coarse_u_y,  coarse_u_z,
                coarse_uxux, coarse_uxuy, coarse_uxuz,
                coarse_uyuy, coarse_uyuz, coarse_uzuz,
                Ntime, Ndepth, Nlat, Nlon,
                longitude, latitude, mask);
        write_energy_transfer(energy_transfer, Iscale, Ntime, Ndepth, Nlat, Nlon);
        write_KE(fine_KE, Iscale, Ntime, Ndepth, Nlat, Nlon);
        #endif

        // Now that we've filtered at the previous scale,
        //   just keep the coarse part for the next iteration
        for (int Itime = 0; Itime < Ntime; Itime++) {
            for (int Idepth = 0; Idepth < Ndepth; Idepth++) {
                for (int Ilat = 0; Ilat < Nlat; Ilat++) {
                    for (int Ilon = 0; Ilon < Nlon; Ilon++) {

                        // Convert our four-index to a one-index
                        index = Index(Itime, Idepth, Ilat, Ilon,
                                      Ntime, Ndepth, Nlat, Nlon);
                        mask_index = Index(0,     0,      Ilat, Ilon,
                                           Ntime, Ndepth, Nlat, Nlon);

                        if (mask.at(mask_index) == 1) { // Skip land areas

                            #if DEBUG >= 4
                            fprintf(stdout, "          Line %d of %s\n", __LINE__, __FILE__);
                            #endif
                            vel_Spher_to_Cart(u_x_tmp,           u_y_tmp,             u_z_tmp,
                                              coarse_u_r.at(index), coarse_u_lon.at(index), coarse_u_lat.at(index),
                                              longitude.at(Ilon),   latitude.at(Ilat));

                            u_x.at(index) = u_x_tmp;
                            u_y.at(index) = u_y_tmp;
                            u_z.at(index) = u_z_tmp;
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

                    if (mask.at(mask_index) == 1) { // Skip land areas

                        vel_Cart_to_Spher(
                                u_r_tmp,    u_lon_tmp,  u_lat_tmp,
                                u_x.at(index), u_y.at(index), u_z.at(index),
                                longitude.at(Ilon),   latitude.at(Ilat));

                        coarse_u_r.at(  index) = u_r_tmp;
                        coarse_u_lon.at(index) = u_lon_tmp;
                        coarse_u_lat.at(index) = u_lat_tmp;
                    }
                }
            }
        }
    }

    write_to_output(coarse_u_r, coarse_u_lon, coarse_u_lat, 
            Nscales, Ntime, Ndepth, Nlat, Nlon);

    #if COMP_VORT
    // Compute and write vorticity
    compute_vorticity(vort_r, vort_lon, vort_lat,
            fine_u_r, fine_u_lon, fine_u_lat,
            Ntime, Ndepth, Nlat, Nlon,
            longitude, latitude, mask);
    write_vorticity(vort_r, vort_lon, vort_lat,
            Nscales, Ntime, Ndepth, Nlat, Nlon);
    #endif

}
