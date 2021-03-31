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
        const dataset & source_data,
        const std::vector<double> & scales,
        const MPI_Comm comm
        ) {

    // Create some tidy names for variables
    const std::vector<double>   &time       = source_data.time,
                                &depth      = source_data.depth,
                                &latitude   = source_data.latitude,
                                &longitude  = source_data.longitude,
                                &dAreas     = source_data.areas;

    const std::vector<bool> &mask = source_data.mask;

    const std::vector<int>  &myCounts = source_data.myCounts,
                            &myStarts = source_data.myStarts;

    const std::vector<double>   &full_u_r   = source_data.variables.at("u_r"),
                                &full_u_lon = source_data.variables.at("u_lon"),
                                &full_u_lat = source_data.variables.at("u_lat"),
                                &full_rho   = constants::COMP_BC_TRANSFERS ? source_data.variables.at("rho") : std::vector<double>(),
                                &full_p     = constants::COMP_BC_TRANSFERS ? source_data.variables.at("p")   : std::vector<double>();

    // Get some MPI info
    int wRank, wSize;
    MPI_Comm_rank( comm, &wRank );
    MPI_Comm_size( comm, &wSize );

    // If we've passed the DO_TIMING flag, then create some timing vars
    Timing_Records timing_records;
    double clock_on;

    // Get dimension sizes
    const int   Nscales = scales.size(),
                Ntime   = source_data.Ntime,
                Ndepth  = source_data.Ndepth,
                Nlat    = source_data.Nlat,
                Nlon    = source_data.Nlon;
    const unsigned int num_pts = Ntime * Ndepth * Nlat * Nlon;

    const int OMP_chunksize = get_omp_chunksize(Nlat,Nlon);

    char fname [50];
    
    const int ndims = 4;
    size_t starts[ndims] = { size_t(myStarts.at(0)), size_t(myStarts.at(1)), size_t(myStarts.at(2)), size_t(myStarts.at(3))};
    size_t counts[ndims] = { size_t(Ntime),          size_t(Ndepth),         size_t(Nlat),           size_t(Nlon)};
    std::vector<std::string> vars_to_write;

    // Preset some post-processing variables
    std::vector<const std::vector<double>*> postprocess_fields;
    std::vector<std::string> postprocess_names;

    int LAT_lb, LAT_ub;

    std::vector<double> local_kernel(Nlat * Nlon);

    std::vector<double> u_x(num_pts), u_y(num_pts), u_z(num_pts);
    std::vector<double> coarse_u_r(num_pts), coarse_u_lon(num_pts), coarse_u_lat(num_pts);

    if (not(constants::MINIMAL_OUTPUT)) {
        vars_to_write.push_back("coarse_u_r");
    }
    if (not(constants::NO_FULL_OUTPUTS)) {
        vars_to_write.push_back("coarse_u_lon");
        vars_to_write.push_back("coarse_u_lat");
        vars_to_write.push_back("KE_from_coarse_vels");
    }
    postprocess_names.push_back( "coarse_u_r");
    postprocess_fields.push_back(&coarse_u_r);

    postprocess_names.push_back( "coarse_u_lon");
    postprocess_fields.push_back(&coarse_u_lon);

    postprocess_names.push_back( "coarse_u_lat");
    postprocess_fields.push_back(&coarse_u_lat);

    int index, horiz_index, Itime, Idepth, Ilat, Ilon, tid;
    // Now convert the Spherical velocities to Cartesian
    //   (although we will still be on a spherical
    //     coordinate system)
    vel_Spher_to_Cart(u_x, u_y, u_z,
            full_u_r, full_u_lon, full_u_lat,
            mask, time, depth, latitude, longitude);

    std::vector<double> full_KE(num_pts, 0.), KE_from_coarse_vel(num_pts, 0.);
    postprocess_names.push_back( "KE_from_coarse_vel");
    postprocess_fields.push_back(&KE_from_coarse_vel);
    KE_from_vels(full_KE, &u_x, &u_y, &u_z, mask);

    // Compute the kernal alpha value (for baroclinic transfers)
    const double kern_alpha = kernel_alpha();

    // Now prepare to filter
    double scale,
           u_x_tmp,     u_y_tmp,   u_z_tmp,
           u_r_tmp,     u_lon_tmp, u_lat_tmp,
           u_x_tilde,   u_y_tilde, u_z_tilde;

    std::vector<double> fine_u_r, fine_u_lon, fine_u_lat,
        div_J, fine_KE, filtered_KE;

    div_J.resize(num_pts);
    postprocess_names.push_back( "div_Jtransport");
    postprocess_fields.push_back(&div_J);
    if (not(constants::MINIMAL_OUTPUT)) {
        vars_to_write.push_back("div_Jtransport");
    }

    if (not(constants::MINIMAL_OUTPUT)) {
        fine_u_r.resize(  num_pts);
        fine_u_lon.resize(num_pts);
        fine_u_lat.resize(num_pts);
        vars_to_write.push_back("fine_u_r");
        vars_to_write.push_back("fine_u_lon");
        vars_to_write.push_back("fine_u_lat");
    }

    // If we're computing transfers, then we already have what
    //   we need to computed band-filtered KE, so might as well do it
    filtered_KE.resize(num_pts);
    if (not(constants::MINIMAL_OUTPUT)) {
        vars_to_write.push_back("filtered_KE");
    }

    postprocess_names.push_back( "filtered_KE");
    postprocess_fields.push_back(&filtered_KE);

    std::vector<double> fine_vort_r, fine_vort_lat, fine_vort_lon,
        coarse_vort_r, coarse_vort_lon, coarse_vort_lat,
        div, OkuboWeiss;
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
        postprocess_names.push_back( "coarse_vort_r");
        postprocess_fields.push_back(&coarse_vort_r);

        if (not(constants::MINIMAL_OUTPUT)) {
            fine_vort_r.resize(  num_pts);
            fine_vort_lat.resize(num_pts);
            fine_vort_lon.resize(num_pts);
            vars_to_write.push_back("fine_vort_r");
        }

        div.resize(num_pts);
        OkuboWeiss.resize(num_pts);
        if (not(constants::NO_FULL_OUTPUTS)) {
            vars_to_write.push_back("coarse_vel_div");
            vars_to_write.push_back("OkuboWeiss");
        }

        postprocess_names.push_back( "coarse_vel_div" );
        postprocess_fields.push_back(&div);

        postprocess_names.push_back( "OkuboWeiss" );
        postprocess_fields.push_back(&OkuboWeiss);

        #if DEBUG >= 1
        if (wRank == 0) { fprintf(stdout, "   ... done.\n"); }
        #endif
    }


    double uxux_tmp, uxuy_tmp, uxuz_tmp, uyuy_tmp, uyuz_tmp, uzuz_tmp;
    double KE_tmp;
    std::vector<double> coarse_uxux, coarse_uxuy, coarse_uxuz,
        coarse_uyuy, coarse_uyuz, coarse_uzuz, coarse_u_x,
        coarse_u_y, coarse_u_z, energy_transfer;
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

        // Fine KE (tau(u,u))
        fine_KE.resize(num_pts);

        postprocess_names.push_back( "fine_KE");
        postprocess_fields.push_back(&fine_KE);

        if (not(constants::NO_FULL_OUTPUTS)) {
            vars_to_write.push_back("fine_KE");
        }

        // Also an array for the transfer itself
        energy_transfer.resize(num_pts);
        if (not(constants::NO_FULL_OUTPUTS)) {
            vars_to_write.push_back("energy_transfer");
        }
        postprocess_names.push_back( "energy_transfer");
        postprocess_fields.push_back(&energy_transfer);
        #if DEBUG >= 1
        if (wRank == 0) { fprintf(stdout, "   ... done.\n"); }
        #endif
    }


    double rho_tmp, p_tmp;
    std::vector<double> coarse_rho, coarse_p, fine_rho, fine_p,
        lambda_rot, lambda_nonlin, lambda_full, PEtoKE, 
        tilde_u_r,    tilde_u_lon,    tilde_u_lat,
        tilde_vort_r, tilde_vort_lon, tilde_vort_lat;
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

        // tilde vorticity
        tilde_vort_r.resize(  num_pts);
        tilde_vort_lon.resize(num_pts);
        tilde_vort_lat.resize(num_pts);
        if (not(constants::NO_FULL_OUTPUTS)) {
            vars_to_write.push_back("tilde_vort_r");
        }
        postprocess_names.push_back( "tilde_vort_r");
        postprocess_fields.push_back(&tilde_vort_r);

        lambda_rot.resize(num_pts);
        lambda_nonlin.resize(num_pts);
        lambda_full.resize(num_pts);
        if (not(constants::NO_FULL_OUTPUTS)) {
            vars_to_write.push_back("Lambda_rotational");
            vars_to_write.push_back("Lambda_nonlinear");
            vars_to_write.push_back("Lambda_full");
        }

        tilde_u_r.resize(  num_pts);
        tilde_u_lon.resize(num_pts);
        tilde_u_lat.resize(num_pts);
        if (not(constants::NO_FULL_OUTPUTS)) {
            vars_to_write.push_back("tilde_u_r");
            vars_to_write.push_back("tilde_u_lon");
            vars_to_write.push_back("tilde_u_lat");
        }

        // We'll need vorticity, so go ahead and compute it
        compute_vorticity(
                coarse_vort_r, coarse_vort_lon, coarse_vort_lat,
                div, OkuboWeiss,
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
    std::vector<double*> filtered_vals, tilde_vals;
    std::vector<bool> filt_use_mask;
    std::vector<const std::vector<double>*> filter_fields;

    filter_fields.push_back(&u_x);
    filt_use_mask.push_back(true);

    filter_fields.push_back(&u_y);
    filt_use_mask.push_back(true);

    filter_fields.push_back(&u_z);
    filt_use_mask.push_back(true);

    filter_fields.push_back(&full_KE);
    filt_use_mask.push_back(true);

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
                    dAreas, vars_to_write, fname, scales.at(Iscale));

            // Add some attributes to the file
            add_attr_to_file("kernel_alpha", kern_alpha, fname);
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
                filter_fields, filt_use_mask, \
                timing_records, clock_on, \
                longitude, latitude, dAreas, scale,\
                full_KE, filtered_KE, fine_KE, \
                full_u_r, full_u_lon, full_u_lat,\
                coarse_u_r, coarse_u_lon, coarse_u_lat,\
                tilde_u_r,  tilde_u_lon,  tilde_u_lat,\
                coarse_u_x, coarse_u_y,   coarse_u_z,\
                coarse_uxux, coarse_uxuy, coarse_uxuz,\
                coarse_uyuy, coarse_uyuz, coarse_uzuz,\
                full_rho, full_p, coarse_rho, coarse_p,\
                fine_rho, fine_p, PEtoKE,\
                fine_u_r, fine_u_lon, fine_u_lat, perc_base)\
        private(Itime, Idepth, Ilat, Ilon, index, horiz_index,\
                u_x_tmp, u_y_tmp, u_z_tmp,\
                u_x_tilde, u_y_tilde, u_z_tilde,\
                u_r_tmp, u_lat_tmp, u_lon_tmp,\
                uxux_tmp, uxuy_tmp, uxuz_tmp,\
                uyuy_tmp, uyuz_tmp, uzuz_tmp,\
                KE_tmp, rho_tmp, p_tmp,\
                LAT_lb, LAT_ub, tid, filtered_vals, tilde_vals ) \
        firstprivate(perc, wRank, local_kernel, perc_count)
        {

            filtered_vals.clear();

            filtered_vals.push_back(&u_x_tmp);
            filtered_vals.push_back(&u_y_tmp);
            filtered_vals.push_back(&u_z_tmp);

            filtered_vals.push_back(&KE_tmp);

            if (constants::COMP_BC_TRANSFERS) {
                filtered_vals.push_back(&rho_tmp);
                filtered_vals.push_back(&p_tmp);
            }

            tilde_vals.clear();

            tilde_vals.push_back(&u_x_tilde);
            tilde_vals.push_back(&u_y_tilde);
            tilde_vals.push_back(&u_z_tilde);

            tilde_vals.push_back(NULL);

            if (constants::COMP_BC_TRANSFERS) {
                tilde_vals.push_back(NULL);
                tilde_vals.push_back(NULL);
            }

            #pragma omp for collapse(2) schedule(guided, OMP_chunksize)
            for (Ilat = 0; Ilat < Nlat; Ilat++) {
                for (Ilon = 0; Ilon < Nlon; Ilon++) {

                    get_lat_bounds(LAT_lb, LAT_ub, latitude,  Ilat, scale); 
                    horiz_index = Index(0,     0,      Ilat, Ilon,
                                        Ntime, Ndepth, Nlat, Nlon);

                    #if DEBUG >= 0
                    tid = omp_get_thread_num();
                    if ( (tid == 0) and (wRank == 0) ) {
                        // Every perc_base percent, print a dot, but only the first thread
                        if ( ((double)(horiz_index+1) / (Nlon*Nlat)) * 100 >= perc ) {
                            perc_count++;
                            if (perc_count % 5 == 0) { fprintf(stdout, "|"); }
                            else                     { fprintf(stdout, "."); }
                            fflush(stdout);
                            perc += perc_base;
                        }
                    }
                    #endif


                    if (constants::DO_TIMING) { clock_on = MPI_Wtime(); }
                    std::fill(local_kernel.begin(), local_kernel.end(), 0);
                    compute_local_kernel(
                            local_kernel, scale, longitude, latitude,
                            Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon,
                            LAT_lb, LAT_ub);
                    if (constants::DO_TIMING) { 
                        timing_records.add_to_record(MPI_Wtime() - clock_on,
                                "kernel_precomputation");
                    }

                    for (Itime = 0; Itime < Ntime; Itime++) {
                        for (Idepth = 0; Idepth < Ndepth; Idepth++) {

                            // Convert our four-index to a one-index
                            index = Index(Itime, Idepth, Ilat, Ilon,
                                          Ntime, Ndepth, Nlat, Nlon);

                            if ( mask.at(index) ) { // Skip land areas

                                // Apply the filter at the point
                                if (constants::DO_TIMING) { clock_on = MPI_Wtime(); }
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
                                vel_Cart_to_Spher_at_point(
                                        u_r_tmp, u_lon_tmp, u_lat_tmp,
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
                                filtered_KE.at(index) = KE_tmp;

                                if (constants::DO_TIMING) { 
                                    timing_records.add_to_record(MPI_Wtime() - clock_on,
                                            "filter_main");
                                }

                                // If we want energy transfers (Pi), 
                                // then do those calculations now
                                if (constants::DO_TIMING) { clock_on = MPI_Wtime(); }
                                if (constants::COMP_TRANSFERS) {

                                    apply_filter_at_point_for_quadratics(
                                            uxux_tmp, uxuy_tmp, uxuz_tmp,
                                            uyuy_tmp, uyuz_tmp, uzuz_tmp,
                                            u_x,      u_y,      u_z,
                                            Ntime, Ndepth, Nlat, Nlon,
                                            Itime, Idepth, Ilat, Ilon,
                                            longitude, latitude, LAT_lb, LAT_ub,
                                            dAreas, scale, mask, &local_kernel);

                                    vel_Spher_to_Cart_at_point(
                                            u_x_tmp, u_y_tmp, u_z_tmp,
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

                                    // tau(u,u)
                                    fine_KE.at(index) = 
                                        0.5 * constants::rho0 * (
                                                uxux_tmp - u_x_tmp * u_x_tmp
                                            +   uyuy_tmp - u_y_tmp * u_y_tmp
                                            +   uzuz_tmp - u_z_tmp * u_z_tmp
                                        );
                                }
                                if (constants::DO_TIMING) { 
                                    timing_records.add_to_record(MPI_Wtime() - clock_on,
                                            "filter_for_Pi");
                                }

                                // If we want baroclinic transfers (Lees and Aluie, 2019), 
                                //    then do those calculations now
                                if (constants::DO_TIMING) { clock_on = MPI_Wtime(); }
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

                                    //
                                    // If we have rho, then also compute tilde fields
                                    //
                                    apply_filter_at_point(
                                            tilde_vals, filter_fields,
                                            Ntime, Ndepth, Nlat, Nlon,
                                            Itime, Idepth, Ilat, Ilon,
                                            longitude, latitude, LAT_lb, LAT_ub,
                                            dAreas, scale, mask, filt_use_mask,
                                            &local_kernel, &full_rho);

                                    vel_Cart_to_Spher_at_point(
                                            u_r_tmp,    u_lon_tmp, u_lat_tmp,
                                            u_x_tilde,  u_y_tilde, u_z_tilde,
                                            longitude.at(Ilon), latitude.at(Ilat));

                                    tilde_u_r.at(  index) = u_r_tmp   / rho_tmp;
                                    tilde_u_lon.at(index) = u_lon_tmp / rho_tmp;
                                    tilde_u_lat.at(index) = u_lat_tmp / rho_tmp;
                                }
                                if (constants::DO_TIMING) { 
                                    timing_records.add_to_record(MPI_Wtime() - clock_on,
                                            "filter_for_Lambda");
                                }

                            }  // end if(masked) block
                            else { // if not masked
                                if (constants::DO_TIMING) { clock_on = MPI_Wtime(); }

                                coarse_u_r.at(  index) = constants::fill_value;
                                coarse_u_lon.at(index) = constants::fill_value;
                                coarse_u_lat.at(index) = constants::fill_value;

                                if (not(constants::MINIMAL_OUTPUT)) {
                                    fine_u_r.at(  index) = constants::fill_value;
                                    fine_u_lon.at(index) = constants::fill_value;
                                    fine_u_lat.at(index) = constants::fill_value;

                                }
                                filtered_KE.at(index) = constants::fill_value;

                                if (constants::COMP_TRANSFERS) {
                                    coarse_u_x.at(index) = constants::fill_value;
                                    coarse_u_y.at(index) = constants::fill_value;
                                    coarse_u_z.at(index) = constants::fill_value;

                                    fine_KE.at(index) = constants::fill_value;
                                }

                                if (constants::COMP_BC_TRANSFERS) {
                                    coarse_rho.at(index) = constants::fill_value;
                                    coarse_p.at(  index) = constants::fill_value;
                                    PEtoKE.at(    index) = constants::fill_value;
                                    if (not(constants::MINIMAL_OUTPUT)) {
                                        fine_rho.at(index)   = constants::fill_value;
                                        fine_p.at(  index)   = constants::fill_value;
                                    }
                                }

                                if (constants::DO_TIMING) { 
                                    timing_records.add_to_record(MPI_Wtime() - clock_on, "land");
                                }
                            }  // end not(masked) block
                        }  // end for(depth) block
                    }  // end for(time) block
                }  // end for(longitude) block
            }  // end for(latitude) block
        }  // end pragma parallel block
        #if DEBUG >= 0
        if (wRank == 0) { fprintf(stdout, "\n"); }
        #endif

        #if DEBUG >= 2
        fprintf(stdout, "  = Rank %d finished filtering loop =\n", wRank);
        fflush(stdout);
        #endif

        // Get KE from coarse velocities
        KE_from_vels(KE_from_coarse_vel, &coarse_u_r, &coarse_u_lon, &coarse_u_lat, mask);

        // Write to file
        if (constants::DO_TIMING) { clock_on = MPI_Wtime(); }
        if (not(constants::MINIMAL_OUTPUT)) {
            write_field_to_output(coarse_u_r,         "coarse_u_r",   
                    starts, counts, fname, &mask);
        }
        if (not(constants::NO_FULL_OUTPUTS)) {
            write_field_to_output(coarse_u_lon,       "coarse_u_lon", 
                    starts, counts, fname, &mask);
            write_field_to_output(coarse_u_lat,       "coarse_u_lat", 
                    starts, counts, fname, &mask);
            write_field_to_output(KE_from_coarse_vel, "KE_from_coarse_vels", 
                    starts, counts, fname, &mask);
        }

        if (not(constants::MINIMAL_OUTPUT)) {
            write_field_to_output(fine_u_r,   "fine_u_r",   starts, counts, fname, &mask);
            write_field_to_output(fine_u_lon, "fine_u_lon", starts, counts, fname, &mask);
            write_field_to_output(fine_u_lat, "fine_u_lat", starts, counts, fname, &mask);

        }
        if (not(constants::MINIMAL_OUTPUT)) {
            write_field_to_output(filtered_KE, "filtered_KE", starts, counts, fname, &mask);
        }
        if (constants::DO_TIMING) { 
            timing_records.add_to_record(MPI_Wtime() - clock_on, "writing");
        }

        if (constants::COMP_VORT) {
            // Compute and write vorticity
            if (constants::DO_TIMING) { clock_on = MPI_Wtime(); }

            #if DEBUG >= 1
            if (wRank == 0) { fprintf(stdout, "Starting compute_vorticity\n"); }
            fflush(stdout);
            #endif
            if (not(constants::MINIMAL_OUTPUT)) {
                compute_vorticity(fine_vort_r, fine_vort_lon, fine_vort_lat,
                        div, OkuboWeiss,
                        fine_u_r, fine_u_lon, fine_u_lat,
                        Ntime, Ndepth, Nlat, Nlon,
                        longitude, latitude, mask);
            }

            compute_vorticity(coarse_vort_r, coarse_vort_lon, coarse_vort_lat,
                    div, OkuboWeiss,
                    coarse_u_r, coarse_u_lon, coarse_u_lat,
                    Ntime, Ndepth, Nlat, Nlon,
                    longitude, latitude, mask);

            if (constants::DO_TIMING) { 
                timing_records.add_to_record(MPI_Wtime() - clock_on, "compute_vorticity");
            }

            if (constants::DO_TIMING) { clock_on = MPI_Wtime(); }
            if (not(constants::MINIMAL_OUTPUT)) {
                write_field_to_output(fine_vort_r, "fine_vort_r", starts, counts, fname, &mask);
                write_field_to_output(div, "coarse_vel_div", starts, counts, fname, &mask);
            }
            if (not(constants::NO_FULL_OUTPUTS)) {
                write_field_to_output(coarse_vort_r, "coarse_vort_r", starts, counts, fname, &mask);
                write_field_to_output(OkuboWeiss, "OkuboWeiss", starts, counts, fname, &mask);
            }
            if (constants::DO_TIMING) { 
                timing_records.add_to_record(MPI_Wtime() - clock_on, "writing");
            }
        }

        if (constants::COMP_TRANSFERS) {
            // Compute the energy transfer through the filter scale
            if (constants::DO_TIMING) { clock_on = MPI_Wtime(); }
            #if DEBUG >= 1
            if (wRank == 0) { fprintf(stdout, "Starting compute_Pi\n"); }
            fflush(stdout);
            #endif
            compute_Pi(
                    energy_transfer, 
                    coarse_u_x,  coarse_u_y,  coarse_u_z,
                    coarse_uxux, coarse_uxuy, coarse_uxuz,
                    coarse_uyuy, coarse_uyuz, coarse_uzuz,
                    Ntime, Ndepth, Nlat, Nlon,
                    longitude, latitude, mask);
            if (constants::DO_TIMING) { 
                timing_records.add_to_record(MPI_Wtime() - clock_on, "compute_Pi");
            }

            if (constants::DO_TIMING) { clock_on = MPI_Wtime(); }
            if (not(constants::NO_FULL_OUTPUTS)) {
                write_field_to_output(energy_transfer, "energy_transfer", 
                        starts, counts, fname, &mask);
                if (not(constants::NO_FULL_OUTPUTS)) {
                    write_field_to_output(fine_KE, "fine_KE", 
                            starts, counts, fname, &mask);
                }
            }
            if (constants::DO_TIMING) { 
                timing_records.add_to_record(MPI_Wtime() - clock_on, "writing");
            }
        }

        if (constants::COMP_BC_TRANSFERS) {
            #if DEBUG >= 1
            if (wRank == 0) { fprintf(stdout, "Starting compute_baroclinic_transfers\n"); }
            fflush(stdout);
            #endif
            if (constants::DO_TIMING) { clock_on = MPI_Wtime(); }
            compute_Lambda_rotational(lambda_rot,
                    coarse_vort_r, coarse_vort_lon, coarse_vort_lat,
                    coarse_rho, coarse_p,
                    Ntime, Ndepth, Nlat, Nlon,
                    longitude, latitude, mask,
                    0.5 * kern_alpha * pow(scales.at(Iscale), 2)
                    );

            compute_Lambda_nonlin_model(lambda_nonlin,
                    coarse_u_r, coarse_u_lon, coarse_u_lat,
                    coarse_rho, coarse_p,
                    Ntime, Ndepth, Nlat, Nlon,
                    longitude, latitude, mask,
                    0.5 * kern_alpha * pow(scales.at(Iscale), 2)
                    );

            compute_Lambda_full(lambda_full,
                    coarse_u_r, coarse_u_lon, coarse_u_lat,
                    tilde_u_r, tilde_u_lon, tilde_u_lat,
                    coarse_p,
                    Ntime, Ndepth, Nlat, Nlon,
                    longitude, latitude, mask
                    );

            compute_vorticity(tilde_vort_r, tilde_vort_lon, tilde_vort_lat,
                    div, OkuboWeiss,
                    tilde_u_r, tilde_u_lon, tilde_u_lat,
                    Ntime, Ndepth, Nlat, Nlon,
                    longitude, latitude, mask);

            if (constants::DO_TIMING) { 
                timing_records.add_to_record(MPI_Wtime() - clock_on, "compute_Lambda");
            }

            if (constants::DO_TIMING) { clock_on = MPI_Wtime(); }
            if (not(constants::NO_FULL_OUTPUTS)) {
                write_field_to_output(lambda_rot,  "Lambda_rotational",   
                        starts, counts, fname, &mask);
                write_field_to_output(lambda_nonlin,  "Lambda_nonlinear",   
                        starts, counts, fname, &mask);
                write_field_to_output(lambda_full,"Lambda_full",   
                        starts, counts, fname, &mask);
                write_field_to_output(PEtoKE,     "PEtoKE",     
                        starts, counts, fname, &mask);
                write_field_to_output(coarse_rho, "coarse_rho", 
                        starts, counts, fname, &mask);
                write_field_to_output(coarse_p,   "coarse_p",   
                        starts, counts, fname, &mask);
                write_field_to_output(tilde_vort_r,   "tilde_vort_p",   
                        starts, counts, fname, &mask);
            }
            if (not(constants::MINIMAL_OUTPUT)) {
                write_field_to_output(fine_rho, "fine_rho", starts, counts, fname, &mask);
                write_field_to_output(fine_p,   "fine_p",   starts, counts, fname, &mask);
            }
            if (constants::DO_TIMING) { 
                timing_records.add_to_record(MPI_Wtime() - clock_on, "writing");
            }
        }

        #if DEBUG >= 1
        if (wRank == 0) { fprintf(stdout, "Starting compute_div_transport\n"); }
        fflush(stdout);
        #endif
        if (constants::DO_TIMING) { clock_on = MPI_Wtime(); }
        compute_div_transport(
                div_J,
                coarse_u_x,  coarse_u_y,  coarse_u_z,
                coarse_uxux, coarse_uxuy, coarse_uxuz,
                coarse_uyuy, coarse_uyuz, coarse_uzuz,
                coarse_p, longitude, latitude,
                Ntime, Ndepth, Nlat, Nlon,
                mask);
        if (constants::DO_TIMING) { 
            timing_records.add_to_record(MPI_Wtime() - clock_on, "compute_transport");
        }

        if (constants::DO_TIMING) { clock_on = MPI_Wtime(); }
        if (not(constants::MINIMAL_OUTPUT)) {
            write_field_to_output(div_J, "div_Jtransport", starts, counts, fname, &mask);
        }
        if (constants::DO_TIMING) { 
            timing_records.add_to_record(MPI_Wtime() - clock_on, "writing");
        }

        //
        //// on-line postprocessing, if desired
        //

        if (constants::APPLY_POSTPROCESS) {
            MPI_Barrier(MPI_COMM_WORLD);

            if (wRank == 0) { fprintf(stdout, "Beginning post-process routines\n"); }
            fflush(stdout);

            if (constants::DO_TIMING) { clock_on = MPI_Wtime(); }
            Apply_Postprocess_Routines(
                    source_data, postprocess_fields, postprocess_names, OkuboWeiss,
                    scales.at(Iscale), "postprocess");
            if (constants::DO_TIMING) { 
                timing_records.add_to_record(MPI_Wtime() - clock_on, "postprocess");
            }
        }

        #if DEBUG >= 0
        // Flushing stdout is necessary for SLURM outputs.
        fflush(stdout);
        #endif

        // If we're doing timings, then print out and reset values now
        if (constants::DO_TIMING) { 
            timing_records.print();
            timing_records.reset();
            fflush(stdout);
        }

    }  // end for(scale) block
} // end filtering
