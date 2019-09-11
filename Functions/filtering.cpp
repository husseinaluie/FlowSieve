#include <math.h>
#include <algorithm>
#include <vector>
#include <omp.h>
#include <mpi.h>
#include "../functions.hpp"
#include "../netcdf_io.hpp"
#include "../constants.hpp"
#include "../postprocess.hpp"

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

    double timing_filt_main, timing_land, timing_kern_precomp,
           timing_filt_comp_transfers, timing_filt_bc_transfers,
           timing_comp_vorticity, timing_comp_pi, timing_writing,
           timing_comp_Lambda, timing_comp_transport,
           clock_on, clock_off;
    //double timing_initialize, timing_writing, timing_comp_vel;
    if (constants::DO_TIMING) {
        // If we've passed the DO_TIMING flag, then create some timing vars
        timing_filt_main = 0.;
        timing_land = 0.;
        timing_kern_precomp = 0.;
        timing_filt_comp_transfers = 0.;
        timing_filt_bc_transfers = 0.;
        timing_comp_vorticity = 0.;
        timing_comp_pi = 0.;
        timing_writing = 0.;
        timing_comp_Lambda = 0.;
        timing_comp_transport = 0.;
    }

    // Get dimension sizes
    const int Nscales = scales.size();
    const int Ntime   = myCounts.at(0);
    const int Ndepth  = myCounts.at(1);
    const int Nlat    = myCounts.at(2);
    const int Nlon    = myCounts.at(3);

    const unsigned int num_pts = Ntime * Ndepth * Nlat * Nlon;
    char fname [50];
    
    const int ndims = 4;
    size_t starts[ndims] = {
        size_t(myStarts.at(0)), size_t(myStarts.at(1)), 
        size_t(myStarts.at(2)), size_t(myStarts.at(3))};
    size_t counts[ndims] = {
        size_t(Ntime), size_t(Ndepth), 
        size_t(Nlat), size_t(Nlon)};
    std::vector<std::string> vars_to_write;

    int LAT_lb, LAT_ub;

    #if DEBUG >= 1
    if (wRank == 0) { fprintf(stdout, "Converting to Cartesian velocities.\n"); }
    #endif

    std::vector<double> local_kernel(Nlat * Nlon);

    std::vector<double> u_x(num_pts);
    std::vector<double> u_y(num_pts);
    std::vector<double> u_z(num_pts);

    std::vector<double> coarse_u_r(  num_pts);
    std::vector<double> coarse_u_lon(num_pts);
    std::vector<double> coarse_u_lat(num_pts);
    if (not(constants::NO_FULL_OUTPUTS)) {
        vars_to_write.push_back("coarse_u_r");
        vars_to_write.push_back("coarse_u_lon");
        vars_to_write.push_back("coarse_u_lat");
    }

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
                                u_x.at(index),      
                                u_y.at(index),      
                                u_z.at(index),
                                full_u_r.at(index), 
                                full_u_lon.at(index), 
                                full_u_lat.at(index),
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

    std::vector<double> fine_u_r, fine_u_lon, fine_u_lat,
        div_J, fine_KE, coarse_KE;

    div_J.resize(num_pts);
    if (not(constants::NO_FULL_OUTPUTS)) {
        vars_to_write.push_back("div_Jtransport");
    }

    if (not(constants::MINIMAL_OUTPUT)) {
        fine_u_r.resize(  num_pts);
        fine_u_lon.resize(num_pts);
        fine_u_lat.resize(num_pts);
        vars_to_write.push_back("fine_u_r");
        vars_to_write.push_back("fine_u_lon");
        vars_to_write.push_back("fine_u_lat");

        // If we're computing transfers, then we already have what
        //   we need to computed band-filtered KE, so might as well do it
        fine_KE.resize(num_pts);
        coarse_KE.resize(num_pts);
        vars_to_write.push_back("fine_KE");
        vars_to_write.push_back("coarse_KE");
    }

    std::vector<double> fine_vort_r, fine_vort_lat, fine_vort_lon,
        coarse_vort_r, coarse_vort_lon, coarse_vort_lat;
    if (constants::COMP_VORT) {
        #if DEBUG >= 1
        if (wRank == 0) { fprintf(stdout, "Initializing COMP_VORT fields.\n"); }
        #endif

        coarse_vort_r.resize(  num_pts);
        coarse_vort_lon.resize(num_pts);
        coarse_vort_lat.resize(num_pts);
        if (not(constants::NO_FULL_OUTPUTS)) {
            vars_to_write.push_back("coarse_vort_r");
        }

        if (not(constants::MINIMAL_OUTPUT)) {
            fine_vort_r.resize(  num_pts);
            fine_vort_lat.resize(num_pts);
            fine_vort_lon.resize(num_pts);
            vars_to_write.push_back("fine_vort_r");
        }

        #if DEBUG >= 1
        if (wRank == 0) { fprintf(stdout, "   ... done.\n"); }
        #endif
    }


    double uxux_tmp, uxuy_tmp, uxuz_tmp, uyuy_tmp, uyuz_tmp, uzuz_tmp;
    double KE_tmp;
    std::vector<double> coarse_uxux, coarse_uxuy, coarse_uxuz,
        coarse_uyuy, coarse_uyuz, coarse_uzuz, coarse_u_x,
        coarse_u_y, coarse_u_z, div, energy_transfer;
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

        if (not(constants::MINIMAL_OUTPUT)) {
            div.resize(num_pts);
            vars_to_write.push_back("full_vel_div");
            vars_to_write.push_back("coarse_vel_div");
        }

        // Also an array for the transfer itself
        energy_transfer.resize(num_pts);
        if (not(constants::NO_FULL_OUTPUTS)) {
            vars_to_write.push_back("energy_transfer");
        }
        #if DEBUG >= 1
        if (wRank == 0) { fprintf(stdout, "   ... done.\n"); }
        #endif
    }


    double rho_tmp, p_tmp;
    std::vector<double> coarse_rho, coarse_p, fine_rho, fine_p,
        lambda_m, PEtoKE;
    if (constants::COMP_BC_TRANSFERS) {
        #if DEBUG >= 1
        if (wRank == 0) { fprintf(stdout, "Initializing COMP_BC_TRANSFERS fields.\n"); }
        #endif
        coarse_rho.resize(num_pts);
        coarse_p.resize(  num_pts);
        if (not(constants::NO_FULL_OUTPUTS)) {
            vars_to_write.push_back("coarse_rho");
            vars_to_write.push_back("coarse_p");
        }

        if (not(constants::MINIMAL_OUTPUT)) {
            fine_rho.resize(  num_pts);
            fine_p.resize(    num_pts);
            vars_to_write.push_back("fine_rho");
            vars_to_write.push_back("fine_p");
        }

        lambda_m.resize(num_pts);
        if (not(constants::NO_FULL_OUTPUTS)) {
            vars_to_write.push_back("Lambda_m");
        }

        // We'll need vorticity, so go ahead and compute it
        compute_vorticity(coarse_vort_r, coarse_vort_lon, coarse_vort_lat,
                full_u_r, full_u_lon, full_u_lat,
                Ntime, Ndepth, Nlat, Nlon,
                longitude, latitude, mask);

        // We also have what we need to compute the PEtoKE term
        //    rho_bar * g * w_bar
        PEtoKE.resize(num_pts);
        if (not(constants::NO_FULL_OUTPUTS)) {
            vars_to_write.push_back("PEtoKE");
        }
        #if DEBUG >= 1
        if (wRank == 0) { fprintf(stdout, "   ... done.\n"); }
        #endif
    }

    int perc_base = 5;
    int perc, perc_count=0;

    // Set up filtering vectors
    std::vector<double*> filtered_vals;
    std::vector<bool> filt_use_mask;
    std::vector<const std::vector<double>*> filter_fields;

    filter_fields.push_back(&u_x);
    filt_use_mask.push_back(true);

    filter_fields.push_back(&u_y);
    filt_use_mask.push_back(true);

    filter_fields.push_back(&u_z);
    filt_use_mask.push_back(true);

    if (not(constants::MINIMAL_OUTPUT)) {
        filter_fields.push_back(&full_KE);
        filt_use_mask.push_back(true);
    }

    if (constants::COMP_BC_TRANSFERS) {
        filter_fields.push_back(&full_rho);
        filt_use_mask.push_back(false);

        filter_fields.push_back(&full_p);
        filt_use_mask.push_back(false);
    }

    //
    //// Begin the main filtering loop
    //
    #if DEBUG>=1
    if (wRank == 0) { fprintf(stdout, "Beginning main filtering loop.\n\n"); }
    #endif
    for (int Iscale = 0; Iscale < Nscales; Iscale++) {

        // Create the output file
        snprintf(fname, 50, "filter_%.6gkm.nc", scales.at(Iscale)/1e3);
        if (not(constants::NO_FULL_OUTPUTS)) {
            initialize_output_file(time, depth, longitude, latitude, 
                    mask, vars_to_write, fname, scales.at(Iscale));
        }

        #if DEBUG >= 0
        if (wRank == 0) { 
            fprintf(stdout, "\nScale %d of %d (%.5g km)\n", 
                Iscale+1, Nscales, scales.at(Iscale)/1e3); 
        }
        #endif

        scale = scales.at(Iscale);
        perc  = perc_base;

        #if DEBUG >= 1
        if (wRank == 0) { fprintf(stdout, "  filtering: "); }
        fflush(stdout);
        #endif

        #pragma omp parallel \
        default(none) \
        shared(mask, u_x, u_y, u_z, stdout, \
                filter_fields, filt_use_mask,\
                timing_filt_main, timing_filt_comp_transfers,\
                timing_filt_bc_transfers, clock_on, clock_off,\
                timing_land, timing_kern_precomp,\
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
                LAT_lb, LAT_ub, tid, filtered_vals) \
        firstprivate(perc, wRank, local_kernel, perc_count)
        {

            filtered_vals.resize(0);

            filtered_vals.push_back(&u_x_tmp);
            filtered_vals.push_back(&u_y_tmp);
            filtered_vals.push_back(&u_z_tmp);

            if (not(constants::MINIMAL_OUTPUT)) {
                filtered_vals.push_back(&KE_tmp);
            }

            if (constants::COMP_BC_TRANSFERS) {
                filtered_vals.push_back(&rho_tmp);
                filtered_vals.push_back(&p_tmp);
            }

            #pragma omp for collapse(2) schedule(dynamic)
            for (Ilat = 0; Ilat < Nlat; Ilat++) {
                for (Ilon = 0; Ilon < Nlon; Ilon++) {

                    get_lat_bounds(LAT_lb, LAT_ub, latitude, Ilat, scale); 
                    mask_index = Index(0,     0,      Ilat, Ilon,
                                       Ntime, Ndepth, Nlat, Nlon);

                    #if DEBUG >= 1
                    tid = omp_get_thread_num();
                    if ( (tid == 0) and (wRank == 0) ) {
                        // Every perc_base percent, print a dot, but only the first thread
                        if ( ((double)(mask_index+1) / (Nlon*Nlat)) * 100 >= perc ) {
                            perc_count++;
                            if (perc_count % 5 == 0) { fprintf(stdout, "|"); }
                            else                     { fprintf(stdout, "."); }
                            fflush(stdout);
                            perc += perc_base;
                        }
                    }
                    #endif

                    if (mask.at(mask_index) == 1) { // Skip land areas

                        if ((constants::DO_TIMING) and (wRank == 0)) { 
                            clock_on = MPI_Wtime(); 
                        }
                        compute_local_kernel(
                                local_kernel, scale, longitude, latitude,
                                Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon);
                        if ((constants::DO_TIMING) and (wRank == 0)) { 
                            clock_off = MPI_Wtime();
                            timing_kern_precomp += clock_off - clock_on;
                        }

                        for (Itime = 0; Itime < Ntime; Itime++) {
                            for (Idepth = 0; Idepth < Ndepth; Idepth++) {

                                // Convert our four-index to a one-index
                                index = Index(Itime, Idepth, Ilat, Ilon,
                                              Ntime, Ndepth, Nlat, Nlon);
    
                                // Apply the filter at the point
                                if ((constants::DO_TIMING) and (wRank == 0)) { 
                                    clock_on = MPI_Wtime(); 
                                }
                                #if DEBUG >= 3
                                if (wRank == 0) { 
                                    fprintf(stdout, "    Line %d of %s\n", 
                                            __LINE__, __FILE__); 
                                }
                                fflush(stdout);
                                #endif

                                apply_filter_at_point(
                                        filtered_vals, filter_fields,
                                        Ntime, Ndepth, Nlat, Nlon,
                                        Itime, Idepth, Ilat, Ilon,
                                        longitude, latitude, LAT_lb, LAT_ub,
                                        dAreas, scale, mask, filt_use_mask,
                                        &local_kernel);

                                // Convert the filtered fields back to spherical
                                #if DEBUG >= 3
                                if (wRank == 0) { 
                                    fprintf(stdout, "          Line %d of %s\n", 
                                            __LINE__, __FILE__); 
                                }
                                fflush(stdout);
                                #endif
                                vel_Cart_to_Spher(u_r_tmp, u_lon_tmp, u_lat_tmp,
                                                  u_x_tmp, u_y_tmp,   u_z_tmp,
                                                  longitude.at(Ilon), latitude.at(Ilat));

                                coarse_u_r.at(  index) = u_r_tmp;
                                coarse_u_lon.at(index) = u_lon_tmp;
                                coarse_u_lat.at(index) = u_lat_tmp;

                                if (not(constants::MINIMAL_OUTPUT)) {
                                    fine_u_r.at(  index) = 
                                        full_u_r.at(  index) - coarse_u_r.at(  index);
                                    fine_u_lon.at(index) = 
                                        full_u_lon.at(index) - coarse_u_lon.at(index);
                                    fine_u_lat.at(index) = 
                                        full_u_lat.at(index) - coarse_u_lat.at(index);
                                }

                                // Also filter KE
                                if (not(constants::MINIMAL_OUTPUT)) {
                                    coarse_KE.at(index) = KE_tmp;
                                    fine_KE.at(index) = 
                                        full_KE.at(index) - coarse_KE.at(index);
                                }

                                if ((constants::DO_TIMING) and (wRank == 0)) { 
                                    clock_off = MPI_Wtime();
                                    timing_filt_main += clock_off - clock_on;
                                }

                                // If we want energy transfers (Pi), 
                                // then do those calculations now
                                if ((constants::DO_TIMING) and (wRank == 0)) { 
                                    clock_on = MPI_Wtime(); 
                                }
                                if (constants::COMP_TRANSFERS) {
                                    apply_filter_at_point_for_quadratics(
                                            uxux_tmp, uxuy_tmp, uxuz_tmp,
                                            uyuy_tmp, uyuz_tmp, uzuz_tmp,
                                            u_x,      u_y,      u_z,
                                            Ntime, Ndepth, Nlat, Nlon,
                                            Itime, Idepth, Ilat, Ilon,
                                            longitude, latitude, LAT_lb, LAT_ub,
                                            dAreas, scale, mask, &local_kernel);

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
                                if ((constants::DO_TIMING) and (wRank == 0)) { 
                                    clock_off = MPI_Wtime();
                                    timing_filt_comp_transfers += clock_off - clock_on;
                                }

                                // If we want baroclinic transfers (Lambda^m), 
                                //    then do those calculations now
                                if ((constants::DO_TIMING) and (wRank == 0)) { 
                                    clock_on = MPI_Wtime(); 
                                }
                                if (constants::COMP_BC_TRANSFERS) {
                                    coarse_rho.at(index) = rho_tmp;
                                    coarse_p.at(  index) = p_tmp;

                                    if (not(constants::MINIMAL_OUTPUT)) {
                                        fine_rho.at(index) = 
                                            full_rho.at(index) - coarse_rho.at(index);
                                        fine_p.at(index)   = 
                                            full_p.at(  index) - coarse_p.at(  index);
                                    }

                                    PEtoKE.at(index) = 
                                        (coarse_rho.at(index) - constants::rho0)
                                        * (-constants::g)
                                        * coarse_u_r.at(index);
                                }
                                if ((constants::DO_TIMING) and (wRank == 0)) { 
                                    clock_off = MPI_Wtime();
                                    timing_filt_bc_transfers += clock_off - clock_on;
                                }

                            }  // end for(depth) block
                        }  // end for(time) block
                    }  // end if(masked) block
                    else { // if not masked
                        if ((constants::DO_TIMING) and (wRank == 0)) { 
                            clock_on = MPI_Wtime(); 
                        }
                        for (Itime = 0; Itime < Ntime; Itime++) {
                            for (Idepth = 0; Idepth < Ndepth; Idepth++) {

                                // Convert our four-index to a one-index
                                index = Index(Itime, Idepth, Ilat, Ilon,
                                              Ntime, Ndepth, Nlat, Nlon);

                                coarse_u_r.at(  index) = constants::fill_value;
                                coarse_u_lon.at(index) = constants::fill_value;
                                coarse_u_lat.at(index) = constants::fill_value;

                                if (not(constants::MINIMAL_OUTPUT)) {
                                    fine_u_r.at(  index) = constants::fill_value;
                                    fine_u_lon.at(index) = constants::fill_value;
                                    fine_u_lat.at(index) = constants::fill_value;

                                    fine_KE.at(   index) = constants::fill_value;
                                    coarse_KE.at( index) = constants::fill_value;
                                }

                                if (constants::COMP_TRANSFERS) {
                                    coarse_u_x.at(index) = constants::fill_value;
                                    coarse_u_y.at(index) = constants::fill_value;
                                    coarse_u_z.at(index) = constants::fill_value;
                                }

                                if (constants::COMP_BC_TRANSFERS) {
                                    coarse_rho.at(index) = constants::fill_value;
                                    coarse_p.at(  index) = constants::fill_value;
                                    PEtoKE.at(index)     = constants::fill_value;
                                    if (not(constants::MINIMAL_OUTPUT)) {
                                        fine_rho.at(index)   = constants::fill_value;
                                        fine_p.at(  index)   = constants::fill_value;
                                    }
                                }
                            }  // end for(depth) block
                        }  // end for(time) block
                        if ((constants::DO_TIMING) and (wRank == 0)) { 
                            clock_off = MPI_Wtime();
                            timing_land += clock_off - clock_on;
                        }
                    }  // end not(masked) block
                }  // end for(longitude) block
            }  // end for(latitude) block
        }  // end pragma parallel block
        #if DEBUG >= 1
        if (wRank == 0) { fprintf(stdout, "\n"); }
        #endif

        #if DEBUG >= 2
        fprintf(stdout, "  = Rank %d finished filtering loop =\n", wRank);
        fflush(stdout);
        #endif

        // Write to file
        if ((constants::DO_TIMING) and (wRank == 0)) { clock_on = MPI_Wtime(); }
        if (not(constants::NO_FULL_OUTPUTS)) {
            write_field_to_output(coarse_u_r,   "coarse_u_r",   
                    starts, counts, fname, &mask);
            write_field_to_output(coarse_u_lon, "coarse_u_lon", 
                    starts, counts, fname, &mask);
            write_field_to_output(coarse_u_lat, "coarse_u_lat", 
                    starts, counts, fname, &mask);
        }

        if (not(constants::MINIMAL_OUTPUT)) {
            write_field_to_output(fine_u_r,   "fine_u_r",   starts, counts, fname, &mask);
            write_field_to_output(fine_u_lon, "fine_u_lon", starts, counts, fname, &mask);
            write_field_to_output(fine_u_lat, "fine_u_lat", starts, counts, fname, &mask);

            write_field_to_output(fine_KE,   "fine_KE",   starts, counts, fname, &mask);
            write_field_to_output(coarse_KE, "coarse_KE", starts, counts, fname, &mask);
        }
        if ((constants::DO_TIMING) and (wRank == 0)) { 
            clock_off = MPI_Wtime();
            timing_writing += clock_off - clock_on;
        }

        if (constants::COMP_VORT) {
            // Compute and write vorticity
            if ((constants::DO_TIMING) and (wRank == 0)) { clock_on = MPI_Wtime(); }

            #if DEBUG >= 1
            if (wRank == 0) { fprintf(stdout, "Starting compute_vorticity\n"); }
            fflush(stdout);
            #endif
            if (not(constants::MINIMAL_OUTPUT)) {
                compute_vorticity(fine_vort_r, fine_vort_lon, fine_vort_lat,
                        fine_u_r, fine_u_lon, fine_u_lat,
                        Ntime, Ndepth, Nlat, Nlon,
                        longitude, latitude, mask);
            }

            compute_vorticity(coarse_vort_r, coarse_vort_lon, coarse_vort_lat,
                    coarse_u_r, coarse_u_lon, coarse_u_lat,
                    Ntime, Ndepth, Nlat, Nlon,
                    longitude, latitude, mask);

            if ((constants::DO_TIMING) and (wRank == 0)) { 
                clock_off = MPI_Wtime();
                timing_comp_vorticity += clock_off - clock_on;
            }

            if ((constants::DO_TIMING) and (wRank == 0)) { clock_on = MPI_Wtime(); }
            if (not(constants::MINIMAL_OUTPUT)) {
                write_field_to_output(fine_vort_r, "fine_vort_r", 
                        starts, counts, fname, &mask);
            }
            if (not(constants::NO_FULL_OUTPUTS)) {
                write_field_to_output(coarse_vort_r, "coarse_vort_r", 
                        starts, counts, fname, &mask);
            }
            if ((constants::DO_TIMING) and (wRank == 0)) { 
                clock_off = MPI_Wtime();
                timing_writing += clock_off - clock_on;
            }
        }

        if (constants::COMP_TRANSFERS) {
            // Compute the energy transfer through the filter scale
            if ((constants::DO_TIMING) and (wRank == 0)) { clock_on = MPI_Wtime(); }
            #if DEBUG >= 1
            if (wRank == 0) { 
                fprintf(stdout, "Starting compute_Pi\n"); 
            }
            fflush(stdout);
            #endif
            compute_Pi(
                    energy_transfer, 
                    coarse_u_x,  coarse_u_y,  coarse_u_z,
                    coarse_uxux, coarse_uxuy, coarse_uxuz,
                    coarse_uyuy, coarse_uyuz, coarse_uzuz,
                    Ntime, Ndepth, Nlat, Nlon,
                    longitude, latitude, mask);
            if ((constants::DO_TIMING) and (wRank == 0)) { 
                clock_off = MPI_Wtime();
                timing_comp_pi += clock_off - clock_on;
            }

            if ((constants::DO_TIMING) and (wRank == 0)) { clock_on = MPI_Wtime(); }
            if (not(constants::NO_FULL_OUTPUTS)) {
                write_field_to_output(energy_transfer, "energy_transfer", 
                        starts, counts, fname, &mask);
            }
            if ((constants::DO_TIMING) and (wRank == 0)) { 
                clock_off = MPI_Wtime();
                timing_writing += clock_off - clock_on;
            }

            // Compute the divergence of the coarse field and full field (for comparison)
            if (not(constants::MINIMAL_OUTPUT)) {
                #if DEBUG >= 1
                if (wRank == 0) { fprintf(stdout, "Starting compute_div_vel (coarse)\n"); }
                #endif
                compute_div_vel(div, coarse_u_x, coarse_u_y, coarse_u_z, 
                        longitude, latitude, Ntime, Ndepth, Nlat, Nlon, mask);
                if ((constants::DO_TIMING) and (wRank == 0)) { clock_on = MPI_Wtime(); }
                write_field_to_output(div, "coarse_vel_div", starts, counts, fname, &mask);
                if ((constants::DO_TIMING) and (wRank == 0)) { 
                    clock_off = MPI_Wtime();
                    timing_writing += clock_off - clock_on;
                }

                #if DEBUG >= 1
                if (wRank == 0) { fprintf(stdout, "Starting compute_div_vel (full)\n"); }
                #endif
                compute_div_vel(div, u_x, u_y, u_z, longitude, latitude,
                        Ntime, Ndepth, Nlat, Nlon, mask);
                if ((constants::DO_TIMING) and (wRank == 0)) { clock_on = MPI_Wtime(); }
                write_field_to_output(div, "full_vel_div", starts, counts, fname, &mask);
                if ((constants::DO_TIMING) and (wRank == 0)) { 
                    clock_off = MPI_Wtime();
                    timing_writing += clock_off - clock_on;
                }
            }
        }

        if (constants::COMP_BC_TRANSFERS) {
            #if DEBUG >= 1
            if (wRank == 0) { fprintf(stdout, "Starting compute_baroclinic_transfers\n"); }
            fflush(stdout);
            #endif
            if ((constants::DO_TIMING) and (wRank == 0)) { clock_on = MPI_Wtime(); }
            compute_baroclinic_transfer(lambda_m,
                    coarse_vort_r, coarse_vort_lon, coarse_vort_lat,
                    coarse_rho, coarse_p,
                    Ntime, Ndepth, Nlat, Nlon,
                    longitude, latitude, mask);
            if ((constants::DO_TIMING) and (wRank == 0)) { 
                clock_off = MPI_Wtime();
                timing_comp_Lambda += clock_off - clock_on;
            }

            if ((constants::DO_TIMING) and (wRank == 0)) { clock_on = MPI_Wtime(); }
            if (not(constants::NO_FULL_OUTPUTS)) {
                write_field_to_output(lambda_m,   "Lambda_m",   
                        starts, counts, fname, &mask);
                write_field_to_output(PEtoKE,     "PEtoKE",     
                        starts, counts, fname, &mask);
                write_field_to_output(coarse_rho, "coarse_rho", 
                        starts, counts, fname, &mask);
                write_field_to_output(coarse_p,   "coarse_p",   
                        starts, counts, fname, &mask);
            }
            if (not(constants::MINIMAL_OUTPUT)) {
                write_field_to_output(fine_rho, "fine_rho", starts, counts, fname, &mask);
                write_field_to_output(fine_p,   "fine_p",   starts, counts, fname, &mask);
            }
            if ((constants::DO_TIMING) and (wRank == 0)) { 
                clock_off = MPI_Wtime();
                timing_writing += clock_off - clock_on;
            }
        }

        #if DEBUG >= 1
        if (wRank == 0) { fprintf(stdout, "Starting compute_div_transport\n"); }
        #endif
        if ((constants::DO_TIMING) and (wRank == 0)) { clock_on = MPI_Wtime(); }
        compute_div_transport(
                div_J,
                coarse_u_x,  coarse_u_y,  coarse_u_z,
                coarse_uxux, coarse_uxuy, coarse_uxuz,
                coarse_uyuy, coarse_uyuz, coarse_uzuz,
                coarse_p, longitude, latitude,
                Ntime, Ndepth, Nlat, Nlon,
                mask);
        if ((constants::DO_TIMING) and (wRank == 0)) { 
            clock_off = MPI_Wtime();
            timing_comp_transport += clock_off - clock_on;
        }

        if ((constants::DO_TIMING) and (wRank == 0)) { clock_on = MPI_Wtime(); }
        if (not(constants::NO_FULL_OUTPUTS)) {
            write_field_to_output(div_J, "div_Jtransport", starts, counts, fname, &mask);
        }
        if ((constants::DO_TIMING) and (wRank == 0)) { 
            clock_off = MPI_Wtime();
            timing_writing += clock_off - clock_on;
        }

        #if DEBUG >= 0
        // Flushing stdout is necessary for SLURM outputs.
        fflush(stdout);
        #endif

        // If we're doing timings, then print out and reset values now
        if ((constants::DO_TIMING) and (wRank == 0)) { 
            fprintf(stdout, "  ::Timings::\n");
            fprintf(stdout, "    filt_kern_precomp = %.13g\n", timing_kern_precomp);
            fprintf(stdout, "\n");
            fprintf(stdout, "    filt_main         = %.13g\n", timing_filt_main);
            fprintf(stdout, "    filt_for_pi       = %.13g\n", timing_filt_comp_transfers);
            fprintf(stdout, "    filt_for_lambda   = %.13g\n", timing_filt_bc_transfers);
            fprintf(stdout, "    filt_land         = %.13g\n", timing_land);
            fprintf(stdout, "\n");
            fprintf(stdout, "    writing           = %.13g\n", timing_writing);
            fprintf(stdout, "\n");
            fprintf(stdout, "    comp_vort         = %.13g\n", timing_comp_vorticity);
            fprintf(stdout, "    comp_pi           = %.13g\n", timing_comp_pi);
            fprintf(stdout, "    comp_lambda       = %.13g\n", timing_comp_Lambda);
            fprintf(stdout, "    comp_transport    = %.13g\n", timing_comp_transport);
            fprintf(stdout, "\n");

            timing_filt_main           = 0.;
            timing_land                = 0.;
            timing_kern_precomp        = 0.;
            timing_filt_comp_transfers = 0.;
            timing_filt_bc_transfers   = 0.;
            timing_comp_vorticity      = 0.;
            timing_comp_pi             = 0.;
            timing_comp_Lambda         = 0.;
            timing_comp_transport      = 0.;
        }

        if (constants::APPLY_POSTPROCESS) {
            Apply_Postprocess_Routines(
                    coarse_u_r, coarse_u_lon, coarse_u_lat, 
                    coarse_vort_r, coarse_vort_lon, coarse_vort_lat,
                    div_J, energy_transfer, lambda_m, PEtoKE,
                    time, depth, latitude, longitude,
                    mask, dAreas,
                    myCounts, myStarts,
                    scales.at(Iscale));
        }

    }  // end for(scale) block
} // end filtering
