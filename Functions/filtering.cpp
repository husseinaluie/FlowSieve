#include <math.h>
#include <algorithm>
#include <vector>
#include <omp.h>
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

    size_t starts[] = {0, 0, 0, 0, 0};
    size_t counts[] = {1, (size_t)Ntime, (size_t)Ndepth, (size_t)Nlat, (size_t)Nlon};

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

    double KE_tmp;
    std::vector<double> coarse_KE(num_pts);
    std::vector<double> fine_KE_filt(num_pts);

    // Now convert the Spherical velocities to Cartesian
    //   (although we will still be on a spherical
    //     coordinate system)
    int index, mask_index, Ilat, Ilon, tid;
    for (int Itime = 0; Itime < Ntime; Itime++) {
        for (int Idepth = 0; Idepth < Ndepth; Idepth++) {
            #pragma omp parallel private(Ilat, Ilon, index, mask_index)
            {
                #pragma omp for collapse(2) schedule(dynamic)
                for (Ilat = 0; Ilat < Nlat; Ilat++) {
                    for (Ilon = 0; Ilon < Nlon; Ilon++) {

                        // Convert our four-index to a one-index
                        index = Index(Itime, Idepth, Ilat, Ilon,
                                      Ntime, Ndepth, Nlat, Nlon);

                        mask_index = Index(0,     0,      Ilat, Ilon,
                                           Ntime, Ndepth, Nlat, Nlon);

                        if (mask.at(mask_index) == 1) { // Skip land areas

                            vel_Spher_to_Cart(     u_x.at(index),      u_y.at(  index),      u_z.at(  index),
                                              full_u_r.at(index), full_u_lon.at(index), full_u_lat.at(index),
                                              longitude.at(Ilon), latitude.at(Ilat));

                            coarse_KE.at(index) = 0.5 * ( 
                                    pow(u_x.at(index), 2) 
                                    + pow(u_y.at(index), 2) 
                                    + pow(u_z.at(index), 2) );
                        }
                    }
                }
            }
        }
    }

    // Now prepare to filter
    double scale,
           u_x_tmp, u_y_tmp,   u_z_tmp,  // The local coarse-graining result
           u_r_tmp, u_lon_tmp, u_lat_tmp;
    
    // Create the output file
    #if DEBUG >= 1
    fprintf(stdout, "Creating output file.\n");
    #endif
    initialize_output_file(time, depth, longitude, latitude, scales, mask);

    // Add additional files to the output file as necessary
    const char* dim_names[] = {"scale","time", "depth", "latitude", "longitude"};

    std::vector<double> fine_u_r(  num_pts);
    std::vector<double> fine_u_lon(num_pts);
    std::vector<double> fine_u_lat(num_pts);

    #if COMP_VORT
    std::vector<double> fine_vort_r(  num_pts);
    std::vector<double> fine_vort_lat(num_pts);
    std::vector<double> fine_vort_lon(num_pts);

    add_var_to_file("vort_r", dim_names, 5);
    add_var_to_file("vort_lon", dim_names, 5);
    add_var_to_file("vort_lat", dim_names, 5);
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
    add_var_to_file("energy_transfer", dim_names, 5);

    // If we're computing transfers, then we already have what
    //   we need to computed band-filtered KE, so might as well do it
    std::vector<double> fine_KE(num_pts);
    add_var_to_file("KE", dim_names, 5);
    #endif

    #if COMP_BC_TRANSFERS
    std::vector<double> coarse_vort_r(num_pts);
    std::vector<double> coarse_vort_lon(num_pts);
    std::vector<double> coarse_vort_lat(num_pts);

    std::vector<double> coarse_rho(full_rho);
    std::vector<double> coarse_p(full_p);
    std::vector<double> fine_rho(full_rho);
    std::vector<double> fine_p(full_p);
    add_var_to_file("rho", dim_names, 5);
    add_var_to_file("p", dim_names, 5);

    std::vector<double> baroclinic_transfer(num_pts);
    add_var_to_file("baroclinic_transfer", dim_names, 5);

    compute_vorticity(coarse_vort_r, coarse_vort_lon, coarse_vort_lat,
            full_u_r, full_u_lon, full_u_lat,
            Ntime, Ndepth, Nlat, Nlon,
            longitude, latitude, mask);

    double rho_tmp, p_tmp;

    // We also have what we need to compute the PEtoKE term
    //    rho_bar * g * w_bar
    std::vector<double> PEtoKE(num_pts);
    add_var_to_file("PEtoKE", dim_names, 5);
    #endif

    #if DEBUG >= 1
    int perc_base = 10;
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

        for (int Itime = 0; Itime < Ntime; Itime++) {

            #if DEBUG >= 0
            fprintf(stdout, "  Time %d of %d\n", Itime+1, Ntime);
            #endif
            
            for (int Idepth = 0; Idepth < Ndepth; Idepth++) {

                #if DEBUG >= 0
                fprintf(stdout, "    Depth %d of %d\n", Idepth+1, Ndepth);
                #endif
                perc = perc_base;

                #pragma omp parallel \
                default(none) \
                shared(Itime, Idepth, mask, u_x, u_y, u_z,\
                        stdout,\
                        longitude, latitude, dAreas, scale,\
                        coarse_KE, fine_KE_filt,\
                        coarse_u_r, coarse_u_lon, coarse_u_lat,\
                        fine_u_r, fine_u_lon, fine_u_lat, perc_base)\
                private(Ilat, Ilon, index, mask_index,\
                        u_x_tmp, u_y_tmp, u_z_tmp,\
                        u_r_tmp, u_lat_tmp, u_lon_tmp,\
                        KE_tmp, tid) \
                firstprivate(perc)
                {
                    #pragma omp for collapse(2) schedule(dynamic)
                    for (Ilat = 0; Ilat < Nlat; Ilat++) {
                        for (Ilon = 0; Ilon < Nlon; Ilon++) {

                            #if DEBUG >= 3
                            // Need to group these together because of pragma collapse syntax
                            if (Ilon == 0) { fprintf(stdout, "      Ilat %d of %d\n", Ilat+1, Nlat); }
                            fprintf(stdout, "        Ilon %d of %d\n", Ilon+1, Nlon);
                            #endif

                            // Convert our four-index to a one-index
                            index = Index(Itime, Idepth, Ilat, Ilon,
                                          Ntime, Ndepth, Nlat, Nlon);
    
                            mask_index = Index(0,     0,      Ilat, Ilon,
                                               Ntime, Ndepth, Nlat, Nlon);

                            #if DEBUG >= 1
                            tid = omp_get_thread_num();
                            if (tid == 0) {
                                // Every 10 percent, print a dot, but only the first thread
                                if ( ((double)(mask_index+1) / (Nlat*Nlon)) * 100 >= perc ) {
                                    if (perc == perc_base) { fprintf(stdout, "      "); }
                                    fprintf(stdout, ".");
                                    fflush(stdout);
                                    perc += perc_base;
                                }
                            }
                            #endif
    
                            if (mask.at(mask_index) == 1) { // Skip land areas
    
                                // Apply the filter at the point
                                #if DEBUG >= 3
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
                                #if DEBUG >= 3
                                fprintf(stdout, "          Line %d of %s\n", __LINE__, __FILE__);
                                #endif
                                vel_Cart_to_Spher(u_r_tmp, u_lon_tmp, u_lat_tmp,
                                                  u_x_tmp, u_y_tmp,   u_z_tmp,
                                                  longitude.at(Ilon), latitude.at(Ilat));

                                // Also filter KE
                                apply_filter_at_point(
                                        KE_tmp, coarse_KE,     
                                        Ntime,  Ndepth, Nlat, Nlon,
                                        Itime,  Idepth, Ilat, Ilon,
                                        longitude, latitude,
                                        dAreas, scale, mask, true);
                                fine_KE_filt.at(index) = coarse_KE.at(index) - KE_tmp;

                                // Subtract current coarse from preceeding coarse to
                                //    get current fine
                                #if DEBUG >= 3
                                fprintf(stdout, "          Line %d of %s\n", __LINE__, __FILE__);
                                #endif
                                fine_u_r.at(  index) = coarse_u_r.at(  index) - u_r_tmp;
                                fine_u_lon.at(index) = coarse_u_lon.at(index) - u_lon_tmp;
                                fine_u_lat.at(index) = coarse_u_lat.at(index) - u_lat_tmp;

                                // Now pass the new coarse along as the preceeding coarse.
                                #if DEBUG >= 3
                                fprintf(stdout, "          Line %d of %s\n", __LINE__, __FILE__);
                                #endif
                                coarse_u_r.at(  index) = u_r_tmp;
                                coarse_u_lon.at(index) = u_lon_tmp;
                                coarse_u_lat.at(index) = u_lat_tmp;

                            }  // end if(masked) block
                        }  // end for(longitude) block
                    }  // end for(latitude) block
                }  // end pragma parallel block
                #if DEBUG >= 1
                #pragma omp master
                {
                    fprintf(stdout, "\n");
                }
                #endif


                #if COMP_TRANSFERS
                perc = perc_base;
                // If we want energy transfers (Pi), then loop through again to do it.
                //   It needs to be in a separate loop because the list of shared
                //   variables is dependent on the pre-processor flags
                #pragma omp parallel \
                default(none) \
                shared(Itime, Idepth, mask, stdout,\
                        longitude, latitude, dAreas, scale,\
                        coarse_u_r, coarse_u_lat, coarse_u_lon,\
                        coarse_u_x, coarse_u_y, coarse_u_z,\
                        coarse_uxux, coarse_uxuy, coarse_uxuz,\
                        coarse_uyuy, coarse_uyuz, coarse_uzuz,\
                        u_x, u_y, u_z, fine_KE, perc_base) \
                private(Ilat, Ilon, index, mask_index,\
                        u_x_tmp, u_y_tmp, u_z_tmp,\
                        uxux_tmp, uxuy_tmp, uxuz_tmp, uyuy_tmp, uyuz_tmp, uzuz_tmp,\
                        tid)\
                firstprivate(perc)
                {
                    #pragma omp for collapse(2) schedule(dynamic)
                    for (Ilat = 0; Ilat < Nlat; Ilat++) {
                        for (Ilon = 0; Ilon < Nlon; Ilon++) {

                            // Convert our four-index to a one-index
                            index = Index(Itime, Idepth, Ilat, Ilon,
                                          Ntime, Ndepth, Nlat, Nlon);
    
                            mask_index = Index(0,     0,      Ilat, Ilon,
                                               Ntime, Ndepth, Nlat, Nlon);

                            #if DEBUG >= 1
                            tid = omp_get_thread_num();
                            if (tid == 0) {
                                // Every 10 percent, print a dot, but only the first thread
                                if ( ((double)(mask_index+1) / (Nlat*Nlon)) * 100 >= perc ) {
                                    if (perc == perc_base) { fprintf(stdout, "      "); }
                                    fprintf(stdout, ".");
                                    fflush(stdout);
                                    perc += perc_base;
                                }
                            }
                            #endif

                            if (mask.at(mask_index) == 1) { // Skip land areas
                                #if DEBUG >= 3
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

                                vel_Spher_to_Cart(u_x_tmp, u_y_tmp, u_z_tmp,
                                                  coarse_u_r.at(index), 
                                                  coarse_u_lat.at(index),  
                                                  coarse_u_lon.at(index),
                                                  longitude.at(Ilon), latitude.at(Ilat));

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
                            }  // end if(masked) block
                        }  // end for(longitude) block
                    }  // end for(latitude) block
                } // end pragma parallel block
                #if DEBUG >= 1
                #pragma omp master
                {
                    fprintf(stdout, "\n");
                }
                #endif
                #endif
                
                #if COMP_BC_TRANSFERS
                // If we want baroclinic transfers, then loop through again
                // to compute as necessary
                //   It needs to be in a separate loop because the list of shared
                //   variables is dependent on the pre-processor flags
                #pragma omp parallel \
                default(none) \
                private(index, mask_index, Ilat, Ilon,\
                        rho_tmp, p_tmp) \
                shared(Itime, Idepth, scale, mask, dAreas, \
                        longitude, latitude, full_rho, full_p,\
                        coarse_rho, coarse_p, fine_rho, fine_p, \
                        PEtoKE, coarse_u_r,\
                        perc_base, tid, stdout)\
                firstprivate(perc)
                {
                    #pragma omp for collapse(2) schedule(dynamic)
                    for (Ilat = 0; Ilat < Nlat; Ilat++) {
                        for (Ilon = 0; Ilon < Nlon; Ilon++) {

                            // Convert our four-index to a one-index
                            index = Index(Itime, Idepth, Ilat, Ilon,
                                          Ntime, Ndepth, Nlat, Nlon);
    
                            mask_index = Index(0,     0,      Ilat, Ilon,
                                               Ntime, Ndepth, Nlat, Nlon);

                            #if DEBUG >= 1
                            tid = omp_get_thread_num();
                            if (tid == 0) {
                                // Every 10 percent, print a dot, but only the first thread
                                if ( ((double)(mask_index+1) / (Nlat*Nlon)) * 100 >= perc ) {
                                    if (perc == perc_base) { fprintf(stdout, "      "); }
                                    fprintf(stdout, ".");
                                    fflush(stdout);
                                    perc += perc_base;
                                }
                            }
                            #endif

                            if (mask.at(mask_index) == 1) { // Skip land areas
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

                                fine_rho.at(index) = full_rho.at(index) - coarse_rho.at(index);
                                fine_p.at(index)   = full_p.at(  index) - coarse_p.at(  index);

                                PEtoKE.at(index) =   coarse_rho.at(index)
                                                   * (-constants::g)
                                                   * coarse_u_r.at(index);

                            }  // end if(masked) block
                        }  // end for(longitude) block
                    }  // end for(latitude) block
                }  // end pragma parallel block
                #if DEBUG >= 1
                #pragma omp master
                {
                    fprintf(stdout, "\n");
                }
                #endif
                #endif
            }  // end for(depth) block
        }  // end for(time) block

        // Write to file
        starts[0] = Iscale;
        write_field_to_output(fine_u_r,   "u_r",   starts, counts);
        write_field_to_output(fine_u_lon, "u_lon", starts, counts);
        write_field_to_output(fine_u_lat, "u_lat", starts, counts);

        write_field_to_output(fine_KE_filt, "KE_filt", starts, counts);

        #if COMP_VORT
        // Compute and write vorticity
        compute_vorticity(fine_vort_r, fine_vort_lon, fine_vort_lat,
                fine_u_r, fine_u_lon, fine_u_lat,
                Ntime, Ndepth, Nlat, Nlon,
                longitude, latitude, mask);
        write_field_to_output(fine_vort_r,   "vort_r",   starts, counts);
        write_field_to_output(fine_vort_lon, "vort_lon", starts, counts);
        write_field_to_output(fine_vort_lat, "vort_lat", starts, counts);
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
        write_field_to_output(energy_transfer, "energy_transfer", starts, counts);
        write_field_to_output(fine_KE, "KE", starts, counts);
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
        write_field_to_output(baroclinic_transfer, "baroclinic_transfer", starts, counts);
        write_field_to_output(PEtoKE, "PEtoKE", starts, counts);
        write_field_to_output(fine_rho, "rho", starts, counts);
        write_field_to_output(fine_p,   "p", starts, counts);
        #endif

        // Now that we've filtered at the previous scale,
        //   just keep the coarse part for the next iteration
        for (int Itime = 0; Itime < Ntime; Itime++) {
            for (int Idepth = 0; Idepth < Ndepth; Idepth++) {
                #pragma omp parallel \
                default(none) \
                shared(u_x, u_y, u_z, longitude, latitude, coarse_u_r, coarse_u_lon, \
                        coarse_u_lat, coarse_KE, fine_KE_filt, mask, Itime, Idepth,\
                        stdout) \
                private(Ilat, Ilon, index, mask_index, u_x_tmp, u_y_tmp, u_z_tmp)
                {
                    #pragma omp for collapse(2) schedule(dynamic)
                    for (Ilat = 0; Ilat < Nlat; Ilat++) {
                        for (Ilon = 0; Ilon < Nlon; Ilon++) {

                            // Convert our four-index to a one-index
                            index = Index(Itime, Idepth, Ilat, Ilon,
                                          Ntime, Ndepth, Nlat, Nlon);
                            mask_index = Index(0,     0,      Ilat, Ilon,
                                               Ntime, Ndepth, Nlat, Nlon);

                            if (mask.at(mask_index) == 1) { // Skip land areas

                                #if DEBUG >= 3
                                fprintf(stdout, "          Line %d of %s\n", __LINE__, __FILE__);
                                #endif
                                vel_Spher_to_Cart(u_x_tmp,           u_y_tmp,             u_z_tmp,
                                                  coarse_u_r.at(index), coarse_u_lon.at(index), coarse_u_lat.at(index),
                                                  longitude.at(Ilon),   latitude.at(Ilat));

                                u_x.at(index) = u_x_tmp;
                                u_y.at(index) = u_y_tmp;
                                u_z.at(index) = u_z_tmp;

                                coarse_KE.at(index) = coarse_KE.at(index) - fine_KE_filt.at(index);

                            }  // end if(masked) block
                        }  // end for(longitude) block
                    }  // end for(latitude) block
                }  // end pragma parallel block
            }  // end for(depth) block
        }  // end for(time) block

        #if DEBUG >= 0
        // Flushing stdout is necessary for SLURM outputs.
        fflush(stdout);
        #endif

    }  // end for(scale) block

    // Now write the remaining small scales
    for (int Itime = 0; Itime < Ntime; Itime++) {
        for (int Idepth = 0; Idepth < Ndepth; Idepth++) {
            #pragma omp parallel \
            default(none) \
            shared(u_x, u_y, u_z, longitude, latitude, coarse_u_r, coarse_u_lon, \
                    coarse_u_lat, mask, Itime, Idepth) \
            private(Ilat, Ilon, index, mask_index, u_r_tmp, u_lon_tmp, u_lat_tmp)
            {
                #pragma omp for collapse(2) schedule(dynamic)
                for (Ilat = 0; Ilat < Nlat; Ilat++) {
                    for (Ilon = 0; Ilon < Nlon; Ilon++) {

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
                        }  // end if(masked) block
                    }  // end for(longitude) block
                }  // end for(latitude) block
            }  // end pragma parallel block
            #if COMP_TRANSFERS
            #pragma omp parallel \
            default(none) \
            shared(u_x, u_y, u_z, fine_KE, mask, Itime, Idepth) \
            private(Ilat, Ilon, index, mask_index)
            {
                #pragma omp for collapse(2) schedule(dynamic)
                for (Ilat = 0; Ilat < Nlat; Ilat++) {
                    for (Ilon = 0; Ilon < Nlon; Ilon++) {

                        // Convert our four-index to a one-index
                        index = Index(Itime, Idepth, Ilat, Ilon,
                                      Ntime, Ndepth, Nlat, Nlon);
                        mask_index = Index(0,     0,      Ilat, Ilon,
                                           Ntime, Ndepth, Nlat, Nlon);

                        if (mask.at(mask_index) == 1) { // Skip land areas

                            // Technically the coarse KE, but re-use the variable
                            // Following Eyink, G. L., & Aluie, H. (2009). 
                            //     Localness of energy cascade in hydrodynamic turbulence. 
                            //         I. smooth coarse graining. 
                            //     Physics of Fluids, 21(11), 1–9. 
                            // This is the KE above the largest filter scale
                            fine_KE.at(index) = 0.5 * (  pow( u_x.at(index), 2 ) 
                                                       + pow( u_y.at(index), 2 ) 
                                                       + pow( u_z.at(index), 2 ) );
                        }  // end if(masked) block
                    }  // end for(longitude) block
                }  // end for(latitude) block
            }  // end pragma parallel block
            #endif
        }  // end for(depth) block
    }  // end for(time) block

    starts[0] = Nscales;
    write_field_to_output(coarse_u_r,   "u_r",   starts, counts);
    write_field_to_output(coarse_u_lon, "u_lon", starts, counts);
    write_field_to_output(coarse_u_lat, "u_lat", starts, counts);

    write_field_to_output(coarse_KE, "KE_filt",  starts, counts);

    #if COMP_VORT
    // Compute and write vorticity
    //    this is actually coarse vort, not fine, but
    //    COMP_VORT only guarantees that fine_vort variables exist.
    compute_vorticity(fine_vort_r, fine_vort_lon, fine_vort_lat,
            coarse_u_r, coarse_u_lon, coarse_u_lat,
            Ntime, Ndepth, Nlat, Nlon,
            longitude, latitude, mask);
    write_field_to_output(fine_vort_r,   "vort_r",    starts, counts);
    write_field_to_output(fine_vort_lon, "vort_lon",  starts, counts);
    write_field_to_output(fine_vort_lat, "vort_lat",  starts, counts);
    #endif

    #if COMP_TRANSFERS
    write_field_to_output(fine_KE, "KE", starts, counts);
    #endif

    #if COMP_BC_TRANSFERS
    write_field_to_output(coarse_rho, "rho", starts, counts);
    write_field_to_output(coarse_p,   "p", starts, counts);
    #endif
}
