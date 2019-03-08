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

#ifndef COMP_BC_TRANSFERS
    #define COMP_BC_TRANSFERS true
#endif

void filtering(
        const std::vector<double> & full_u_r,       /**< [in] Full u_r velocity array */
        const std::vector<double> & full_u_lon,     /**< [in] Full u_lon velocity array */
        const std::vector<double> & full_u_lat,     /**< [in] Full u_lat velocity array */
        const std::vector<double> & full_rho,       /**< [in] Full density array */
        const std::vector<double> & full_p,         /**< [in] Full pressure array */
        const std::vector<double> & scales,         /**< [in] Array of filtering scales */
        const std::vector<double> & dAreas,         /**< [in] Array of cell areas (2D) (compute_areas()) */
        const std::vector<double> & time,           /**< [in] Time dimension (1D) */
        const std::vector<double> & depth,          /**< [in] Depth dimension (1D) */
        const std::vector<double> & longitude,      /**< [in] Longitude dimension (1D) */
        const std::vector<double> & latitude,       /**< [in] Latitude dimension (1D) */
        const std::vector<double> & mask            /**< [in] Array to distinguish between land and water cells (2D) */
        ) {


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

    // Copy the current full u_r into the initial coarse
    std::vector<double> coarse_u_r(  full_u_r);
    std::vector<double> coarse_u_lon(full_u_lon);
    std::vector<double> coarse_u_lat(full_u_lat);

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

                        vel_Spher_to_Cart(     u_x.at(index),      u_y.at(  index),      u_z.at(  index),
                                          full_u_r.at(index), full_u_lon.at(index), full_u_lat.at(index),
                                          longitude.at(Ilon), latitude.at(Ilat));
                    }
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
           u_r_tmp, u_lon_tmp, u_lat_tmp;

    std::vector<double> fine_u_r(  num_pts);
    std::vector<double> fine_u_lon(num_pts);
    std::vector<double> fine_u_lat(num_pts);

    #if COMP_VORT
    std::vector<double> fine_vort_r(  num_pts);
    std::vector<double> fine_vort_lat(num_pts);
    std::vector<double> fine_vort_lon(num_pts);
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

    #if COMP_BC_TRANSFERS
    std::vector<double> coarse_vort_r(num_pts);
    std::vector<double> coarse_vort_lon(num_pts);
    std::vector<double> coarse_vort_lat(num_pts);

    std::vector<double> coarse_rho(full_rho);
    std::vector<double> coarse_p(full_p);

    std::vector<double> baroclinic_transfer(num_pts);

    compute_vorticity(coarse_vort_r, coarse_vort_lon, coarse_vort_lat,
            full_u_r, full_u_lon, full_u_lat,
            Ntime, Ndepth, Nlat, Nlon,
            longitude, latitude, mask);

    double rho_tmp, p_tmp;
    #endif

    #if DEBUG >= 1
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

        //for (int Itime = 0; Itime < Ntime; Itime++) {
        for (int Itime = 0; Itime < 1; Itime++) {

            #if DEBUG >= 0
            fprintf(stdout, "  Time %d of %d\n", Itime+1, Ntime);
            #endif

            for (int Idepth = 0; Idepth < Ndepth; Idepth++) {

                #if DEBUG >= 0
                fprintf(stdout, "    Depth %d of %d\n", Idepth+1, Ndepth);
                #endif

                #if DEBUG >= 1
                perc = 10;
                fprintf(stdout, "      ");
                fflush(stdout);
                #endif

                for (int Ilat = 0; Ilat < Nlat; Ilat++) {

                    #if DEBUG >= 1
                    // Every 10 percent, print a dot
                    if ( ( ((double) Ilat) / Nlat) * 100 >= perc ) {
                        fprintf(stdout, ".");
                        fflush(stdout);
                        perc += 10;
                    }
                    #endif

                    #if DEBUG >= 3
                    fprintf(stdout, "      Ilat %d of %d\n", Ilat+1, Nlat);
                    #endif

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
                                    u_x_tmp, u_x,     
                                    Ntime,  Ndepth, Nlat, Nlon,
                                    Itime,  Idepth, Ilat, Ilon,
                                    longitude, latitude,
                                    dAreas, scale, mask, true);
                            apply_filter_at_point(
                                    u_y_tmp, u_y,     
                                    Ntime,  Ndepth, Nlat, Nlon,
                                    Itime,  Idepth, Ilat, Ilon,
                                    longitude, latitude,
                                    dAreas, scale, mask, true);
                            apply_filter_at_point(
                                    u_z_tmp, u_z,     
                                    Ntime,  Ndepth, Nlat, Nlon,
                                    Itime,  Idepth, Ilat, Ilon,
                                    longitude, latitude,
                                    dAreas, scale, mask, true);

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

                            #if COMP_BC_TRANSFERS
                            apply_filter_at_point(
                                    rho_tmp, full_rho,     
                                    Ntime,  Ndepth, Nlat, Nlon,
                                    Itime,  Idepth, Ilat, Ilon,
                                    longitude, latitude,
                                    dAreas, scale, mask, false);
                            apply_filter_at_point(
                                    p_tmp, full_p,     
                                    Ntime,  Ndepth, Nlat, Nlon,
                                    Itime,  Idepth, Ilat, Ilon,
                                    longitude, latitude,
                                    dAreas, scale, mask, false);
                            coarse_rho.at(index) = rho_tmp;
                            coarse_p.at(  index) = p_tmp;
                            #endif
                        }
                    }
                }
                #if DEBUG >= 1
                fprintf(stdout, "\n");
                #endif
            }
        }

        // Write to file
        write_field_to_output(fine_u_r,   "u_r",   Iscale, Ntime, Ndepth, Nlat, Nlon);
        write_field_to_output(fine_u_lon, "u_lon", Iscale, Ntime, Ndepth, Nlat, Nlon);
        write_field_to_output(fine_u_lat, "u_lat", Iscale, Ntime, Ndepth, Nlat, Nlon);

        #if COMP_VORT
        // Compute and write vorticity
        compute_vorticity(fine_vort_r, fine_vort_lon, fine_vort_lat,
                fine_u_r, fine_u_lon, fine_u_lat,
                Ntime, Ndepth, Nlat, Nlon,
                longitude, latitude, mask);
        write_field_to_output(fine_vort_r,   "vort_r",   Iscale, Ntime, Ndepth, Nlat, Nlon);
        write_field_to_output(fine_vort_lon, "vort_lon", Iscale, Ntime, Ndepth, Nlat, Nlon);
        write_field_to_output(fine_vort_lat, "vort_lat", Iscale, Ntime, Ndepth, Nlat, Nlon);
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
        write_field_to_output(energy_transfer, "energy_transfer", Iscale, Ntime, Ndepth, Nlat, Nlon);
        write_field_to_output(fine_KE, "KE", Iscale, Ntime, Ndepth, Nlat, Nlon);
        #endif

        #if COMP_BC_TRANSFERS
        compute_vorticity(coarse_vort_r, coarse_vort_lon, coarse_vort_lat,
                coarse_u_r, coarse_u_lon, coarse_u_lat,
                Ntime, Ndepth, Nlat, Nlon,
                longitude, latitude, mask);
        compute_baroclinic_transfer(baroclinic_transfer,
                coarse_vort_r, coarse_vort_lon, coarse_vort_lat,
                coarse_rho, coarse_p,
                Ntime, Ndepth, Nlat, Nlon,
                longitude, latitude, mask);
        write_field_to_output(baroclinic_transfer, "baroclinic_transfer", Iscale, Ntime, Ndepth, Nlat, Nlon);
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

        #if DEBUG >= 0
        // Flushing stdout is necessary for SLURM outputs.
        fflush(stdout);
        #endif

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

    write_field_to_output(coarse_u_r,   "u_r",   Nscales, Ntime, Ndepth, Nlat, Nlon);
    write_field_to_output(coarse_u_lon, "u_lon", Nscales, Ntime, Ndepth, Nlat, Nlon);
    write_field_to_output(coarse_u_lat, "u_lat", Nscales, Ntime, Ndepth, Nlat, Nlon);

    #if COMP_VORT
    // Compute and write vorticity
    compute_vorticity(coarse_vort_r, coarse_vort_lon, coarse_vort_lat,
            coarse_u_r, coarse_u_lon, coarse_u_lat,
            Ntime, Ndepth, Nlat, Nlon,
            longitude, latitude, mask);
    write_field_to_output(coarse_vort_r,   "vort_r",   Nscales, Ntime, Ndepth, Nlat, Nlon);
    write_field_to_output(coarse_vort_lon, "vort_lon", Nscales, Ntime, Ndepth, Nlat, Nlon);
    write_field_to_output(coarse_vort_lat, "vort_lat", Nscales, Ntime, Ndepth, Nlat, Nlon);
    #endif

}
