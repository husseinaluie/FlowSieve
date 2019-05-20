#include <math.h>
#include <algorithm>
#include <vector>
#include <omp.h>
#include <mpi.h>
#include "../functions.hpp"
#include "../netcdf_io.hpp"
#include "../constants.hpp"

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
        const std::vector<double> & mask,           /**< [in] Array to distinguish between land and water cells (2D) */
        const std::vector<int>    & myCounts,       /**< [in] Array of dimension sizes */
        const std::vector<int>    & myStarts,       /**< [in] Array of dimension sizes */
        const MPI_Comm comm                         /**< [in] MPI Communicator */
        ) {

    int wRank, wSize;
    MPI_Comm_rank( comm, &wRank );
    MPI_Comm_size( comm, &wSize );

    // Get dimension sizes
    const int Nscales = scales.size();
    const int Ntime   = myCounts.at(0);
    const int Ndepth  = myCounts.at(1);
    const int Nlat    = myCounts.at(2);
    const int Nlon    = myCounts.at(3);

    const unsigned int num_pts = Ntime * Ndepth * Nlat * Nlon;
    char filename [50];
    
    const int ndims = 4;
    size_t starts[ndims] = {
        size_t(myStarts.at(0)), size_t(myStarts.at(1)), 
        size_t(myStarts.at(2)), size_t(myStarts.at(3))};
    size_t counts[ndims] = {
        size_t(Ntime), size_t(Ndepth), 
        size_t(Nlat), size_t(Nlon)};
    std::vector<std::string> vars_to_write;

    #if DEBUG >= 1
    if (wRank == 0) { fprintf(stdout, "Converting to Cartesian velocities.\n"); }
    #endif

    std::vector<double> distances(Nlat * Nlon);

    std::vector<double> u_x(num_pts);
    std::vector<double> u_y(num_pts);
    std::vector<double> u_z(num_pts);

    std::vector<double> coarse_u_r(  num_pts);
    std::vector<double> coarse_u_lon(num_pts);
    std::vector<double> coarse_u_lat(num_pts);
    vars_to_write.push_back("coarse_u_r");
    vars_to_write.push_back("coarse_u_lon");
    vars_to_write.push_back("coarse_u_lat");

    std::vector<double> full_KE(num_pts);

    // Now convert the Spherical velocities to Cartesian
    //   (although we will still be on a spherical
    //     coordinate system)
    int index, mask_index, Itime, Idepth, Ilat, Ilon, tid;
    #pragma omp parallel private(Itime, Idepth, Ilat, Ilon, index, mask_index)
    {
        #pragma omp for collapse(4) schedule(guided)
        for (Itime = 0; Itime < Ntime; Itime++) {
            for (Idepth = 0; Idepth < Ndepth; Idepth++) {
                for (Ilat = 0; Ilat < Nlat; Ilat++) {
                    for (Ilon = 0; Ilon < Nlon; Ilon++) {

                        // Convert our four-index to a one-index
                        index = Index(Itime, Idepth, Ilat, Ilon,
                                      Ntime, Ndepth, Nlat, Nlon);

                        mask_index = Index(0,     0,      Ilat, Ilon,
                                           Ntime, Ndepth, Nlat, Nlon);

                        if (mask.at(mask_index) == 1) { // Skip land areas
                            vel_Spher_to_Cart(     
                                u_x.at(index),      u_y.at(  index),      u_z.at(  index),
                                full_u_r.at(index), full_u_lon.at(index), full_u_lat.at(index),
                                longitude.at(Ilon), latitude.at(Ilat));

                            full_KE.at(index) = 
                                0.5 * ( 
                                          pow(u_x.at(index), 2) 
                                        + pow(u_y.at(index), 2) 
                                        + pow(u_z.at(index), 2) 
                                      );
                        } // done if land
                    } // done lon loop
                } // done lat loop
            } // done depth loop
        } // done time loop
    } // done pragma
    #if DEBUG >= 1
    if (wRank == 0) { fprintf(stdout, "   ... done.\n"); }
    #endif

    // Now prepare to filter
    double scale,
           u_x_tmp, u_y_tmp,   u_z_tmp,
           u_r_tmp, u_lon_tmp, u_lat_tmp;

    std::vector<double> fine_u_r(  num_pts);
    std::vector<double> fine_u_lon(num_pts);
    std::vector<double> fine_u_lat(num_pts);
    vars_to_write.push_back("fine_u_r");
    vars_to_write.push_back("fine_u_lon");
    vars_to_write.push_back("fine_u_lat");

    std::vector<double> div_J(num_pts);
    vars_to_write.push_back("div_Jtransport");

    std::vector<double> fine_vort_r, fine_vort_lat, fine_vort_lon,
        coarse_vort_r, coarse_vort_lon, coarse_vort_lat;
    if (constants::COMP_VORT) {
        #if DEBUG >= 1
        if (wRank == 0) { fprintf(stdout, "Initializing COMP_VORT fields.\n"); }
        #endif
        fine_vort_r.resize(  num_pts);
        fine_vort_lat.resize(num_pts);
        fine_vort_lon.resize(num_pts);

        coarse_vort_r.resize(  num_pts);
        coarse_vort_lon.resize(num_pts);
        coarse_vort_lat.resize(num_pts);

        vars_to_write.push_back("fine_vort_r");
        vars_to_write.push_back("coarse_vort_r");
        //vars_to_write.push_back("vort_lon");
        //vars_to_write.push_back("vort_lat");
        #if DEBUG >= 1
        if (wRank == 0) { fprintf(stdout, "   ... done.\n"); }
        #endif
    }


    double uxux_tmp, uxuy_tmp, uxuz_tmp, uyuy_tmp, uyuz_tmp, uzuz_tmp;
    double KE_tmp;
    std::vector<double> coarse_uxux, coarse_uxuy, coarse_uxuz,
        coarse_uyuy, coarse_uyuz, coarse_uzuz, coarse_u_x,
        coarse_u_y, coarse_u_z, div, energy_transfer,
        fine_KE, coarse_KE;
    if (constants::COMP_TRANSFERS) {
        #if DEBUG >= 1
        if (wRank == 0) { fprintf(stdout, "Initializing COMP_TRANSFERS fields.\n"); }
        #endif
        // If computing energy transfers, we'll need some more arrays
        coarse_uxux.resize(num_pts);
        coarse_uxuy.resize(num_pts);
        coarse_uxuz.resize(num_pts);
        coarse_uyuy.resize(num_pts);
        coarse_uyuz.resize(num_pts);
        coarse_uzuz.resize(num_pts);

        // We'll also need to keep the coarse velocities
        coarse_u_x.resize(num_pts);
        coarse_u_y.resize(num_pts);
        coarse_u_z.resize(num_pts);

        div.resize(num_pts);
        vars_to_write.push_back("full_vel_div");
        vars_to_write.push_back("coarse_vel_div");

        // Also an array for the transfer itself
        energy_transfer.resize(num_pts);
        vars_to_write.push_back("energy_transfer");

        // If we're computing transfers, then we already have what
        //   we need to computed band-filtered KE, so might as well do it
        fine_KE.resize(num_pts);
        coarse_KE.resize(num_pts);
        vars_to_write.push_back("fine_KE");
        vars_to_write.push_back("coarse_KE");
        #if DEBUG >= 1
        if (wRank == 0) { fprintf(stdout, "   ... done.\n"); }
        #endif
    }


    double rho_tmp, p_tmp;
    std::vector<double> coarse_rho, coarse_p, fine_rho, fine_p,
        baroclinic_transfer, PEtoKE;
    if (constants::COMP_BC_TRANSFERS) {
        #if DEBUG >= 1
        if (wRank == 0) { fprintf(stdout, "Initializing COMP_BC_TRANSFERS fields.\n"); }
        #endif
        coarse_rho.resize(num_pts);
        coarse_p.resize(  num_pts);
        fine_rho.resize(  num_pts);
        fine_p.resize(    num_pts);
        vars_to_write.push_back("fine_rho");
        vars_to_write.push_back("fine_p");
        vars_to_write.push_back("coarse_rho");
        vars_to_write.push_back("coarse_p");

        baroclinic_transfer.resize(num_pts);
        vars_to_write.push_back("Lambda_m");

        compute_vorticity(coarse_vort_r, coarse_vort_lon, coarse_vort_lat,
                full_u_r, full_u_lon, full_u_lat,
                Ntime, Ndepth, Nlat, Nlon,
                longitude, latitude, mask);

        // We also have what we need to compute the PEtoKE term
        //    rho_bar * g * w_bar
        PEtoKE.resize(num_pts);
        vars_to_write.push_back("PEtoKE");
        #if DEBUG >= 1
        if (wRank == 0) { fprintf(stdout, "   ... done.\n"); }
        #endif
    }

    #if DEBUG >= 1
    int perc_base = 5;
    int perc;
    #endif

    //
    //// Begin the main filtering loop
    //
    #if DEBUG>=1
    if (wRank == 0) { fprintf(stdout, "Beginning main filtering loop.\n\n"); }
    #endif
    for (int Iscale = 0; Iscale < Nscales; Iscale++) {

        // Create the output file
        snprintf(filename, 50, "filter_%.6gkm.nc", scales.at(Iscale)/1e3);
        initialize_output_file(time, depth, longitude, latitude, 
                mask, vars_to_write, filename, scales.at(Iscale));

        #if DEBUG >= 0
        if (wRank == 0) { 
            fprintf(stdout, "Scale %d of %d (%.5g km)\n", 
                Iscale+1, Nscales, scales.at(Iscale)/1e3); 
        }
        #endif

        scale = scales.at(Iscale);
        perc  = perc_base;

        #if DEBUG >= 0
        if (wRank == 0) { fprintf(stdout, "  filtering: "); }
        fflush(stdout);
        #endif

        #pragma omp parallel \
        default(none) \
        shared(mask, u_x, u_y, u_z, stdout,\
                longitude, latitude, dAreas, scale,\
                full_KE, coarse_KE, fine_KE,\
                full_u_r, full_u_lon, full_u_lat,\
                coarse_u_r, coarse_u_lon, coarse_u_lat,\
                coarse_u_x, coarse_u_y, coarse_u_z,\
                coarse_uxux, coarse_uxuy, coarse_uxuz,\
                coarse_uyuy, coarse_uyuz, coarse_uzuz,\
                full_rho, full_p, coarse_rho, coarse_p,\
                fine_rho, fine_p, PEtoKE,\
                fine_u_r, fine_u_lon, fine_u_lat, perc_base)\
        private(Itime, Idepth, Ilat, Ilon, index, mask_index,\
                u_x_tmp, u_y_tmp, u_z_tmp,\
                u_r_tmp, u_lat_tmp, u_lon_tmp,\
                uxux_tmp, uxuy_tmp, uxuz_tmp,\
                uyuy_tmp, uyuz_tmp, uzuz_tmp,\
                KE_tmp, rho_tmp, p_tmp,\
                tid) \
        firstprivate(perc, wRank, distances)
        {
            #pragma omp for collapse(2) schedule(dynamic)
            for (Ilat = 0; Ilat < Nlat; Ilat++) {
                for (Ilon = 0; Ilon < Nlon; Ilon++) {

                    compute_distances(
                            distances, longitude, latitude,
                            Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon);

                    mask_index = Index(0,     0,      Ilat, Ilon,
                                       Ntime, Ndepth, Nlat, Nlon);

                    #if DEBUG >= 1
                    tid = omp_get_thread_num();
                    if ( (tid == 0) and (wRank == 0) ) {
                        // Every perc_base percent, print a dot, but only the first thread
                        if ( ((double)(mask_index+1) / (Nlon*Nlat)) * 100 >= perc ) {
                            fprintf(stdout, ".");
                            fflush(stdout);
                            perc += perc_base;
                        }
                    }
                    #endif

                    for (Itime = 0; Itime < Ntime; Itime++) {
                        for (Idepth = 0; Idepth < Ndepth; Idepth++) {

                            // Convert our four-index to a one-index
                            index = Index(Itime, Idepth, Ilat, Ilon,
                                          Ntime, Ndepth, Nlat, Nlon);
                            if (mask.at(mask_index) == 1) { // Skip land areas
    
                                // Apply the filter at the point
                                #if DEBUG >= 3
                                if (wRank == 0) { fprintf(stdout, "    Line %d of %s\n", __LINE__, __FILE__); }
                                #endif
                                apply_filter_at_point(
                                        u_x_tmp, u_x,     
                                        Ntime,  Ndepth, Nlat, Nlon,
                                        Itime,  Idepth, Ilat, Ilon,
                                        longitude, latitude,
                                        dAreas, scale, mask, true,
                                        &distances);
                                apply_filter_at_point(
                                        u_y_tmp, u_y,     
                                        Ntime,  Ndepth, Nlat, Nlon,
                                        Itime,  Idepth, Ilat, Ilon,
                                        longitude, latitude,
                                        dAreas, scale, mask, true,
                                        &distances);
                                apply_filter_at_point(
                                        u_z_tmp, u_z,     
                                        Ntime,  Ndepth, Nlat, Nlon,
                                        Itime,  Idepth, Ilat, Ilon,
                                        longitude, latitude,
                                        dAreas, scale, mask, true,
                                        &distances);

                                // Convert the filtered fields back to spherical
                                #if DEBUG >= 3
                                if (wRank == 0) { fprintf(stdout, "          Line %d of %s\n", __LINE__, __FILE__); }
                                #endif
                                vel_Cart_to_Spher(u_r_tmp, u_lon_tmp, u_lat_tmp,
                                                  u_x_tmp, u_y_tmp,   u_z_tmp,
                                                  longitude.at(Ilon), latitude.at(Ilat));

                                coarse_u_r.at(  index) = u_r_tmp;
                                coarse_u_lon.at(index) = u_lon_tmp;
                                coarse_u_lat.at(index) = u_lat_tmp;

                                fine_u_r.at(  index) = full_u_r.at(  index) - coarse_u_r.at(  index);
                                fine_u_lon.at(index) = full_u_lon.at(index) - coarse_u_lon.at(index);
                                fine_u_lat.at(index) = full_u_lat.at(index) - coarse_u_lat.at(index);

                                // Also filter KE
                                apply_filter_at_point(
                                        KE_tmp, full_KE,
                                        Ntime,  Ndepth, Nlat, Nlon,
                                        Itime,  Idepth, Ilat, Ilon,
                                        longitude, latitude,
                                        dAreas, scale, mask, true,
                                        &distances);
                                coarse_KE.at(index) = KE_tmp;
                                fine_KE.at(index) = full_KE.at(index) - coarse_KE.at(index);

                                // If we want energy transfers (Pi), then do those calculations now
                                if (constants::COMP_TRANSFERS) {
                                    apply_filter_at_point_for_quadratics(
                                            uxux_tmp, uxuy_tmp, uxuz_tmp,
                                            uyuy_tmp, uyuz_tmp, uzuz_tmp,
                                            u_x,      u_y,      u_z,
                                            Ntime,  Ndepth, Nlat, Nlon,
                                            Itime,  Idepth, Ilat, Ilon,
                                            longitude, latitude,
                                            dAreas, scale,
                                            mask, &distances);

                                    vel_Spher_to_Cart(u_x_tmp, u_y_tmp, u_z_tmp,
                                            coarse_u_r.at(index), 
                                            coarse_u_lon.at(index),  
                                            coarse_u_lat.at(index),
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
                                }

                                // If we want baroclinic transfers (Lambda^m), 
                                //    then do those calculations now
                                if (constants::COMP_BC_TRANSFERS) {
                                    apply_filter_at_point(
                                            rho_tmp, full_rho,     
                                            Ntime,  Ndepth, Nlat, Nlon,
                                            Itime,  Idepth, Ilat, Ilon,
                                            longitude, latitude,
                                            dAreas, scale, mask, false,
                                            &distances);
                                    apply_filter_at_point(
                                            p_tmp, full_p,     
                                            Ntime,  Ndepth, Nlat, Nlon,
                                            Itime,  Idepth, Ilat, Ilon,
                                            longitude, latitude,
                                            dAreas, scale, mask, false,
                                            &distances);
                                    coarse_rho.at(index) = rho_tmp;
                                    coarse_p.at(  index) = p_tmp;

                                    fine_rho.at(index) = full_rho.at(index) - coarse_rho.at(index);
                                    fine_p.at(index)   = full_p.at(  index) - coarse_p.at(  index);

                                    PEtoKE.at(index) = (coarse_rho.at(index) - constants::rho0)
                                        * (-constants::g)
                                        * coarse_u_r.at(index);
                                }

                            }  // end if(masked) block
                            else { // if not masked
                                fine_KE.at(   index) = constants::fill_value;
                                fine_u_r.at(  index) = constants::fill_value;
                                fine_u_lon.at(index) = constants::fill_value;
                                fine_u_lat.at(index) = constants::fill_value;

                                coarse_KE.at(   index) = constants::fill_value;
                                coarse_u_r.at(  index) = constants::fill_value;
                                coarse_u_lon.at(index) = constants::fill_value;
                                coarse_u_lat.at(index) = constants::fill_value;

                                if (constants::COMP_TRANSFERS) {
                                    coarse_u_x.at(index) = constants::fill_value;
                                    coarse_u_y.at(index) = constants::fill_value;
                                    coarse_u_z.at(index) = constants::fill_value;
                                }

                                if (constants::COMP_BC_TRANSFERS) {
                                    coarse_rho.at(index) = constants::fill_value;
                                    coarse_p.at(  index) = constants::fill_value;
                                    fine_rho.at(index)   = constants::fill_value;
                                    fine_p.at(  index)   = constants::fill_value;
                                    PEtoKE.at(index)     = constants::fill_value;
                                }

                            }  // end not(masked) block
                        }  // end for(depth) block
                    }  // end for(time) block
                }  // end for(longitude) block
            }  // end for(latitude) block
        }  // end pragma parallel block
        #if DEBUG >= 1
        if (wRank == 0) { fprintf(stdout, "\n"); }
        #endif

        #if DEBUG >= 2
        fprintf(stdout, "  = Rank %d finished filtering loop =\n", wRank);
        #endif

        // Write to file
        write_field_to_output(fine_u_r,     "fine_u_r",     starts, counts, filename);
        write_field_to_output(fine_u_lon,   "fine_u_lon",   starts, counts, filename);
        write_field_to_output(fine_u_lat,   "fine_u_lat",   starts, counts, filename);
        write_field_to_output(coarse_u_r,   "coarse_u_r",   starts, counts, filename);
        write_field_to_output(coarse_u_lon, "coarse_u_lon", starts, counts, filename);
        write_field_to_output(coarse_u_lat, "coarse_u_lat", starts, counts, filename);

        write_field_to_output(fine_KE,   "fine_KE",   starts, counts, filename);
        write_field_to_output(coarse_KE, "coarse_KE", starts, counts, filename);

        if (constants::COMP_VORT) {
            // Compute and write vorticity
            compute_vorticity(fine_vort_r, fine_vort_lon, fine_vort_lat,
                    fine_u_r, fine_u_lon, fine_u_lat,
                    Ntime, Ndepth, Nlat, Nlon,
                    longitude, latitude, mask);
            write_field_to_output(fine_vort_r,   "fine_vort_r",   starts, counts, filename);
            //write_field_to_output(fine_vort_lon, "fine_vort_lon", starts, counts, filename);
            //write_field_to_output(fine_vort_lat, "fine_vort_lat", starts, counts, filename);

            compute_vorticity(coarse_vort_r, coarse_vort_lon, coarse_vort_lat,
                    coarse_u_r, coarse_u_lon, coarse_u_lat,
                    Ntime, Ndepth, Nlat, Nlon,
                    longitude, latitude, mask);
            write_field_to_output(coarse_vort_r,   "coarse_vort_r",   starts, counts, filename);
            //write_field_to_output(coarse_vort_lon, "coarse_vort_lon", starts, counts, filename);
            //write_field_to_output(coarse_vort_lat, "coarse_vort_lat", starts, counts, filename);
        }

        if (constants::COMP_TRANSFERS) {
            // Compute the energy transfer through the filter scale
            compute_energy_transfer_through_scale(
                    energy_transfer, 
                    coarse_u_x,  coarse_u_y,  coarse_u_z,
                    coarse_uxux, coarse_uxuy, coarse_uxuz,
                    coarse_uyuy, coarse_uyuz, coarse_uzuz,
                    Ntime, Ndepth, Nlat, Nlon,
                    longitude, latitude, mask);
            write_field_to_output(energy_transfer, "energy_transfer", starts, counts, filename);

            // Compute the divergence of the coarse field and full field (for comparison)
            #if DEBUG >= 1
            if (wRank == 0) { fprintf(stdout, "Starting compute_div_vel (coarse)\n"); }
            #endif
            compute_div_vel(div, coarse_u_x, coarse_u_y, coarse_u_z, longitude, latitude,
                    Ntime, Ndepth, Nlat, Nlon, mask);
            write_field_to_output(div, "coarse_vel_div", starts, counts, filename);

            #if DEBUG >= 1
            if (wRank == 0) { fprintf(stdout, "Starting compute_div_vel (full)\n"); }
            #endif
            compute_div_vel(div, u_x, u_y, u_z, longitude, latitude,
                    Ntime, Ndepth, Nlat, Nlon, mask);
            write_field_to_output(div, "full_vel_div", starts, counts, filename);
        }

        if (constants::COMP_BC_TRANSFERS) {
            #if DEBUG >= 1
            if (wRank == 0) { fprintf(stdout, "Starting compute_baroclinic_transfers\n"); }
            #endif
            compute_baroclinic_transfer(baroclinic_transfer,
                    coarse_vort_r, coarse_vort_lon, coarse_vort_lat,
                    coarse_rho, coarse_p,
                    Ntime, Ndepth, Nlat, Nlon,
                    longitude, latitude, mask);
            write_field_to_output(baroclinic_transfer, "Lambda_m", starts, counts, filename);
            write_field_to_output(PEtoKE,     "PEtoKE",     starts, counts, filename);
            write_field_to_output(coarse_rho, "coarse_rho", starts, counts, filename);
            write_field_to_output(coarse_p,   "coarse_p",   starts, counts, filename);
            write_field_to_output(fine_rho,   "fine_rho",   starts, counts, filename);
            write_field_to_output(fine_p,     "fine_p",     starts, counts, filename);
        }

        #if DEBUG >= 1
        if (wRank == 0) { fprintf(stdout, "Starting compute_div_transport\n"); }
        #endif
        compute_div_transport(
                div_J,
                coarse_u_x,  coarse_u_y,  coarse_u_z,
                coarse_uxux, coarse_uxuy, coarse_uxuz,
                coarse_uyuy, coarse_uyuz, coarse_uzuz,
                coarse_p, longitude, latitude,
                Ntime, Ndepth, Nlat, Nlon,
                mask);
        write_field_to_output(div_J, "div_Jtransport", starts, counts, filename);

        #if DEBUG >= 0
        // Flushing stdout is necessary for SLURM outputs.
        fflush(stdout);
        #endif

    }  // end for(scale) block
} // end filtering
