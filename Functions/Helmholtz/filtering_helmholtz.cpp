#include <math.h>
#include <algorithm>
#include <vector>
#include <omp.h>
#include <mpi.h>
#include "../../functions.hpp"
#include "../../netcdf_io.hpp"
#include "../../constants.hpp"
#include "../../postprocess.hpp"
#include "../../preprocess.hpp"

void filtering_helmholtz(
        const std::vector<double> & F_potential,
        const std::vector<double> & F_toroidal,
        const std::vector<double> & scales,
        const std::vector<double> & dAreas,
        const std::vector<double> & time,
        const std::vector<double> & depth,
        const std::vector<double> & longitude,
        const std::vector<double> & latitude,
        const std::vector<bool>   & mask,
        const std::vector<int>    & myCounts,
        const std::vector<int>    & myStarts,
        const MPI_Comm comm
        ) {

    int wRank, wSize;
    MPI_Comm_rank( comm, &wRank );
    MPI_Comm_size( comm, &wSize );

    // If we've passed the DO_TIMING flag, then create some timing vars
    Timing_Records timing_records;
    double clock_on;

    // Get dimension sizes
    const int Nscales = scales.size();
    const int Ntime   = myCounts.at(0);
    const int Ndepth  = myCounts.at(1);
    const int Nlat    = myCounts.at(2);
    const int Nlon    = myCounts.at(3);

    const int OMP_chunksize = get_omp_chunksize(Nlat,Nlon);

    const unsigned int num_pts = Ntime * Ndepth * Nlat * Nlon;
    char fname [50];
    
    const int ndims = 4;
    size_t starts[ndims] = {
        size_t(myStarts.at(0)), size_t(myStarts.at(1)), 
        size_t(myStarts.at(2)), size_t(myStarts.at(3))};
    size_t counts[ndims] = {
        size_t(Ntime), size_t(Ndepth), 
        size_t(Nlat), size_t(Nlon)};
    size_t index;
    std::vector<std::string> vars_to_write;

    int LAT_lb, LAT_ub;

    std::vector<double> local_kernel(Nlat * Nlon);

    std::vector<double> null_vector(0);

    std::vector<double> 
        // Arrays to store filtered Phi and Psi fields (potential, pseudo-potential)
        coarse_F_tor(   num_pts, constants::fill_value ), 
        coarse_F_pot(   num_pts, constants::fill_value ),

        // Original KE
        KE_tor_orig(    num_pts, constants::fill_value ),
        KE_pot_orig(    num_pts, constants::fill_value ),
        KE_tot_orig(    num_pts, constants::fill_value ),

        // Coarse KE (computed from velocities)
        KE_tor_coarse(  num_pts, constants::fill_value ),
        KE_pot_coarse(  num_pts, constants::fill_value ),
        KE_tot_coarse(  num_pts, constants::fill_value ),

        // Fine KE ( tau(uu) = bar(uu) - bar(u)bar(u) )
        KE_tor_fine(    num_pts, constants::fill_value ),
        KE_pot_fine(    num_pts, constants::fill_value ),
        KE_tot_fine(    num_pts, constants::fill_value ),

        // Fine KE modified ( uu - bar(u)bar(u) )
        KE_tor_fine_mod(    num_pts, constants::fill_value ),
        KE_pot_fine_mod(    num_pts, constants::fill_value ),
        KE_tot_fine_mod(    num_pts, constants::fill_value ),

        // Filtered KE (used to compute fine KE)
        KE_tor_filt(    num_pts, constants::fill_value ),
        KE_pot_filt(    num_pts, constants::fill_value ),
        KE_tot_filt(    num_pts, constants::fill_value ),

        // Enstrophy
        Enst_tor(       num_pts, constants::fill_value ),
        Enst_pot(       num_pts, constants::fill_value ),
        Enst_tot(       num_pts, constants::fill_value ),

        // Velocity divergences
        div_tor(        num_pts, constants::fill_value ),
        div_pot(        num_pts, constants::fill_value ),
        div_tot(        num_pts, constants::fill_value ),

        // Cartensian velocities
        u_x_tor( num_pts, constants::fill_value ),
        u_y_tor( num_pts, constants::fill_value ),
        u_z_tor( num_pts, constants::fill_value ),

        u_x_pot( num_pts, constants::fill_value ),
        u_y_pot( num_pts, constants::fill_value ),
        u_z_pot( num_pts, constants::fill_value ),

        u_x_tot( num_pts, constants::fill_value ),
        u_y_tot( num_pts, constants::fill_value ),
        u_z_tot( num_pts, constants::fill_value ),

        u_x_coarse( num_pts, constants::fill_value ),
        u_y_coarse( num_pts, constants::fill_value ),
        u_z_coarse( num_pts, constants::fill_value ),

        //
        //// Diadic (Cartesian) velocity components
        //

        // tor
        ux_ux_tor( num_pts, constants::fill_value ),
        ux_uy_tor( num_pts, constants::fill_value ),
        ux_uz_tor( num_pts, constants::fill_value ),
        uy_uy_tor( num_pts, constants::fill_value ),
        uy_uz_tor( num_pts, constants::fill_value ),
        uz_uz_tor( num_pts, constants::fill_value ),

        // pot
        ux_ux_pot( num_pts, constants::fill_value ),
        ux_uy_pot( num_pts, constants::fill_value ),
        ux_uz_pot( num_pts, constants::fill_value ),
        uy_uy_pot( num_pts, constants::fill_value ),
        uy_uz_pot( num_pts, constants::fill_value ),
        uz_uz_pot( num_pts, constants::fill_value ),

        // tot
        ux_ux_tot( num_pts, constants::fill_value ),
        ux_uy_tot( num_pts, constants::fill_value ),
        ux_uz_tot( num_pts, constants::fill_value ),
        uy_uy_tot( num_pts, constants::fill_value ),
        uy_uz_tot( num_pts, constants::fill_value ),
        uz_uz_tot( num_pts, constants::fill_value ),

        //
        //// Spherical velocity components
        //

        // Spherical - radial velocities (just set to zero)
        u_r_zero(       num_pts, 0.),

        // Spherical - zonal velocities
        u_lon_tor(      num_pts, constants::fill_value ),
        u_lon_pot(      num_pts, constants::fill_value ),
        u_lon_tot(      num_pts, constants::fill_value ),

        // Spherical - meridional velocities
        u_lat_tor(      num_pts, constants::fill_value ),
        u_lat_pot(      num_pts, constants::fill_value ),
        u_lat_tot(      num_pts, constants::fill_value ),

        // Vorticity (only r component)
        vort_tor_r(     num_pts, constants::fill_value ),
        vort_pot_r(     num_pts, constants::fill_value ),
        vort_tot_r(     num_pts, constants::fill_value ),

        // Okubo-Weiss values
        OkuboWeiss_tor( num_pts, constants::fill_value ),
        OkuboWeiss_pot( num_pts, constants::fill_value ),
        OkuboWeiss_tot( num_pts, constants::fill_value ),

        // Pi
        Pi_tor( num_pts, constants::fill_value ),
        Pi_pot( num_pts, constants::fill_value ),
        Pi_tot( num_pts, constants::fill_value );

    //
    //// Compute original (unfiltered) KE
    //
     
    // Get pot and tor velocities
    toroidal_vel_from_F( u_lon_tor, u_lat_tor, F_toroidal,
            longitude, latitude, Ntime, Ndepth, Nlat, Nlon, mask);

    potential_vel_from_F(u_lon_pot, u_lat_pot, F_potential,
            longitude, latitude, Ntime, Ndepth, Nlat, Nlon, mask);

    #pragma omp parallel \
    default( none ) \
    shared( KE_tor_orig, KE_pot_orig, KE_tot_orig, mask, \
            u_lon_tor, u_lat_tor, u_lon_pot, u_lat_pot, u_lon_tot, u_lat_tot) \
    private( index )
    {
        #pragma omp for collapse(1) schedule(guided, OMP_chunksize)
        for (index = 0; index < u_lon_tor.size(); ++index) {
            u_lon_tot.at(index) = u_lon_tor.at(index) + u_lon_pot.at(index);
            u_lat_tot.at(index) = u_lat_tor.at(index) + u_lat_pot.at(index);
            if (mask.at(index)) {

                KE_tor_orig.at(index) = 0.5 * constants::rho0 * ( pow(u_lon_tor.at(index), 2.) + pow(u_lat_tor.at(index), 2.) );
                KE_pot_orig.at(index) = 0.5 * constants::rho0 * ( pow(u_lon_pot.at(index), 2.) + pow(u_lat_pot.at(index), 2.) );
                KE_tot_orig.at(index) = 0.5 * constants::rho0 * ( pow(u_lon_tot.at(index), 2.) + pow(u_lat_tot.at(index), 2.) );
            }
        }
    }

    // Get Cartesian velocities, will need them for Pi
    vel_Spher_to_Cart( u_x_tor, u_y_tor, u_z_tor, u_r_zero, u_lon_tor, u_lat_tor, mask, time, depth, latitude, longitude );
    vel_Spher_to_Cart( u_x_pot, u_y_pot, u_z_pot, u_r_zero, u_lon_pot, u_lat_pot, mask, time, depth, latitude, longitude );
    vel_Spher_to_Cart( u_x_tot, u_y_tot, u_z_tot, u_r_zero, u_lon_tot, u_lat_tot, mask, time, depth, latitude, longitude );

    if (not(constants::NO_FULL_OUTPUTS)) {
        //
        // These variables are output unless full outputs are turned off
        // 

        vars_to_write.push_back("coarse_F_tor");
        vars_to_write.push_back("coarse_F_pot");

        vars_to_write.push_back("u_lon_tor");
        vars_to_write.push_back("u_lat_tor");

        vars_to_write.push_back("u_lon_pot");
        vars_to_write.push_back("u_lat_pot");

        vars_to_write.push_back("KE_tor_fine");
        vars_to_write.push_back("KE_pot_fine");
        vars_to_write.push_back("KE_tot_fine");

        vars_to_write.push_back("KE_tor_fine_mod");
        vars_to_write.push_back("KE_pot_fine_mod");
        vars_to_write.push_back("KE_tot_fine_mod");
    }

    if (not(constants::MINIMAL_OUTPUT)) {
        //
        // There's outputs are only included if not set to minimal outputs
        //

        vars_to_write.push_back("div_tor");
        vars_to_write.push_back("div_pot");
        vars_to_write.push_back("div_tot");

        vars_to_write.push_back("OkuboWeiss_tor");
        vars_to_write.push_back("OkuboWeiss_pot");
        vars_to_write.push_back("OkuboWeiss_tot");

        vars_to_write.push_back("Pi_tor");
        vars_to_write.push_back("Pi_pot");
        vars_to_write.push_back("Pi_tot");

        vars_to_write.push_back("KE_tor_filt");
        vars_to_write.push_back("KE_pot_filt");
        vars_to_write.push_back("KE_tot_filt");
    }

    // Compute the kernal alpha value (for baroclinic transfers)
    const double kern_alpha = kernel_alpha();

    // Now prepare to filter
    double scale;
    int Itime, Idepth, Ilat, Ilon, tid;

    int perc_base = 5;
    int perc, perc_count=0;

    //
    //// Set up filtering vectors
    //
    std::vector<double*> filtered_vals;
    std::vector<bool> filt_use_mask;
    std::vector<const std::vector<double>*> filter_fields;

    double F_pot_tmp;
    filter_fields.push_back(&F_potential);
    filt_use_mask.push_back(false);

    double F_tor_tmp;
    filter_fields.push_back(&F_toroidal);
    filt_use_mask.push_back(false);

    double uxux_tmp, uxuy_tmp, uxuz_tmp, uyuy_tmp, uyuz_tmp, uzuz_tmp;

    //
    //// Set up post-processing variables
    //
    std::vector<const std::vector<double>*> postprocess_fields_tor, postprocess_fields_pot, postprocess_fields_tot;
    std::vector<std::string> postprocess_names;

    postprocess_names.push_back( "F" );
    postprocess_fields_tor.push_back( &coarse_F_tor );
    postprocess_fields_pot.push_back( &coarse_F_pot );
    postprocess_fields_tot.push_back( &u_r_zero     );

    postprocess_names.push_back( "Coarse_KE" );
    postprocess_fields_tor.push_back( &KE_tor_coarse );
    postprocess_fields_pot.push_back( &KE_pot_coarse );
    postprocess_fields_tot.push_back( &KE_tot_coarse );

    postprocess_names.push_back( "Fine_KE" );
    postprocess_fields_tor.push_back( &KE_tor_fine );
    postprocess_fields_pot.push_back( &KE_pot_fine );
    postprocess_fields_tot.push_back( &KE_tot_fine );

    postprocess_names.push_back( "Fine_KE_mod" );
    postprocess_fields_tor.push_back( &KE_tor_fine_mod );
    postprocess_fields_pot.push_back( &KE_pot_fine_mod );
    postprocess_fields_tot.push_back( &KE_tot_fine_mod );

    /*
    postprocess_names.push_back( "Filtered_KE" );
    postprocess_fields_tor.push_back( &KE_tor_filt );
    postprocess_fields_pot.push_back( &KE_pot_filt );
    postprocess_fields_tot.push_back( &KE_tot_filt );
    */

    postprocess_names.push_back( "Enstrophy" );
    postprocess_fields_tor.push_back( &Enst_tor );
    postprocess_fields_pot.push_back( &Enst_pot );
    postprocess_fields_tot.push_back( &Enst_tot );

    postprocess_names.push_back( "u_lon" );
    postprocess_fields_tor.push_back( &u_lon_tor );
    postprocess_fields_pot.push_back( &u_lon_pot );
    postprocess_fields_tot.push_back( &u_lon_tot );

    postprocess_names.push_back( "u_lat" );
    postprocess_fields_tor.push_back( &u_lat_tor );
    postprocess_fields_pot.push_back( &u_lat_pot );
    postprocess_fields_tot.push_back( &u_lat_tot );

    postprocess_names.push_back( "OkuboWeiss" );
    postprocess_fields_tor.push_back( &OkuboWeiss_tor );
    postprocess_fields_pot.push_back( &OkuboWeiss_pot );
    postprocess_fields_tot.push_back( &OkuboWeiss_tot );

    postprocess_names.push_back( "Pi" );
    postprocess_fields_tor.push_back( &Pi_tor );
    postprocess_fields_pot.push_back( &Pi_pot );
    postprocess_fields_tot.push_back( &Pi_tot );

    if (not(constants::MINIMAL_OUTPUT)) {
        postprocess_names.push_back( "velocity_divergence" );
        postprocess_fields_tor.push_back( &div_tor );
        postprocess_fields_pot.push_back( &div_pot );
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
        shared(mask, stdout, perc_base, \
                filter_fields, filt_use_mask, \
                timing_records, clock_on, \
                longitude, latitude, dAreas, scale, \
                F_potential, F_toroidal, coarse_F_tor, coarse_F_pot, \
                u_x_tor, u_y_tor, u_z_tor, u_x_pot, u_y_pot, u_z_pot, u_x_tot, u_y_tot, u_z_tot, \
                ux_ux_tor, ux_uy_tor, ux_uz_tor, uy_uy_tor, uy_uz_tor, uz_uz_tor,\
                ux_ux_pot, ux_uy_pot, ux_uz_pot, uy_uy_pot, uy_uz_pot, uz_uz_pot,\
                ux_ux_tot, ux_uy_tot, ux_uz_tot, uy_uy_tot, uy_uz_tot, uz_uz_tot,\
                KE_tor_filt, KE_pot_filt, KE_tot_filt \
                ) \
        private(Itime, Idepth, Ilat, Ilon, index, \
                F_tor_tmp, F_pot_tmp, uxux_tmp, uxuy_tmp, uxuz_tmp, uyuy_tmp, uyuz_tmp, uzuz_tmp, \
                LAT_lb, LAT_ub, tid, filtered_vals) \
        firstprivate(perc, wRank, local_kernel, perc_count)
        {

            filtered_vals.clear();

            filtered_vals.push_back(&F_pot_tmp);
            filtered_vals.push_back(&F_tor_tmp);

            #pragma omp for collapse(1) schedule(guided)
            for (Ilat = 0; Ilat < Nlat; Ilat++) {

                get_lat_bounds(LAT_lb, LAT_ub, latitude,  Ilat, scale); 

                // If our longitude grid is uniform, and spans the full periodic domain,
                // then we can just compute it once and translate it at each lon index
                if (constants::DO_TIMING) { clock_on = MPI_Wtime(); }
                if ( (constants::UNIFORM_LON_GRID) and (constants::FULL_LON_SPAN) ) {
                    std::fill(local_kernel.begin(), local_kernel.end(), 0);
                    compute_local_kernel(
                            local_kernel, scale, longitude, latitude,
                            Ilat, 0, Ntime, Ndepth, Nlat, Nlon,
                            LAT_lb, LAT_ub);
                }
                if (constants::DO_TIMING) { timing_records.add_to_record(MPI_Wtime() - clock_on, "kernel_precomputation"); }

                for (Ilon = 0; Ilon < Nlon; Ilon++) {


                    #if DEBUG >= 0
                    tid = omp_get_thread_num();
                    if ( (tid == 0) and (wRank == 0) ) {
                        // Every perc_base percent, print a dot, but only the first thread
                        if ( ((double)(Ilat*Nlon + Ilon + 1) / (Nlon*Nlat)) * 100 >= perc ) {
                            perc_count++;
                            if (perc_count % 5 == 0) { fprintf(stdout, "|"); }
                            else                     { fprintf(stdout, "."); }
                            fflush(stdout);
                            perc += perc_base;
                        }
                    }
                    #endif

                    if (constants::DO_TIMING) { clock_on = MPI_Wtime(); }
                    if ( not( (constants::UNIFORM_LON_GRID) and (constants::FULL_LON_SPAN) ) ) {
                        // If we couldn't precompute the kernel earlier, then do it now
                        std::fill(local_kernel.begin(), local_kernel.end(), 0);
                        compute_local_kernel(
                                local_kernel, scale, longitude, latitude,
                                Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon,
                                LAT_lb, LAT_ub);
                    } else {
                        // If we were able to compute the kernel earlier, then just translate it now 
                        // to the new longitude index
                        if (Ilon > 0) { roll_field(local_kernel, "lon", 1, 1, 1, Nlat, Nlon); }
                    }
                    if (constants::DO_TIMING) { timing_records.add_to_record(MPI_Wtime() - clock_on, "kernel_precomputation"); }

                    if (constants::DO_TIMING) { clock_on = MPI_Wtime(); }
                    for (Itime = 0; Itime < Ntime; Itime++) {
                        for (Idepth = 0; Idepth < Ndepth; Idepth++) {

                            // Convert our four-index to a one-index
                            index = Index(Itime, Idepth, Ilat, Ilon,
                                          Ntime, Ndepth, Nlat, Nlon);

                            // The F_tor and F_pot fields exist over land from the projection
                            //     procedure, so do those filtering operations on land as well.
                            // The other stuff (KE, etc), will only be done on water cells

                            // Apply the filter at the point
                            apply_filter_at_point(
                                    filtered_vals, filter_fields,
                                    Ntime, Ndepth, Nlat, Nlon,
                                    Itime, Idepth, Ilat, Ilon,
                                    longitude, latitude, LAT_lb, LAT_ub,
                                    dAreas, scale, mask, filt_use_mask,
                                    &local_kernel);

                            // Store the filtered values in the appropriate arrays
                            coarse_F_pot.at(index) = F_pot_tmp;
                            coarse_F_tor.at(index) = F_tor_tmp;

                            if (mask.at(index)) { // Skip land areas

                                //
                                //// Also get (uiuj)_bar from Cartesian velocities
                                //

                                // tor
                                apply_filter_at_point_for_quadratics(
                                        uxux_tmp, uxuy_tmp, uxuz_tmp,
                                        uyuy_tmp, uyuz_tmp, uzuz_tmp,
                                        u_x_tor,  u_y_tor,  u_z_tor,
                                        Ntime, Ndepth, Nlat, Nlon,
                                        Itime, Idepth, Ilat, Ilon,
                                        longitude, latitude, LAT_lb, LAT_ub,
                                        dAreas, scale, mask, &local_kernel);

                                ux_ux_tor.at(index) = uxux_tmp;
                                ux_uy_tor.at(index) = uxuy_tmp;
                                ux_uz_tor.at(index) = uxuz_tmp;
                                uy_uy_tor.at(index) = uyuy_tmp;
                                uy_uz_tor.at(index) = uyuz_tmp;
                                uz_uz_tor.at(index) = uzuz_tmp;

                                KE_tor_filt.at(index) = 0.5 * constants::rho0 * (uxux_tmp + uyuy_tmp + uzuz_tmp);

                                // pot
                                apply_filter_at_point_for_quadratics(
                                        uxux_tmp, uxuy_tmp, uxuz_tmp,
                                        uyuy_tmp, uyuz_tmp, uzuz_tmp,
                                        u_x_pot,  u_y_pot,  u_z_pot,
                                        Ntime, Ndepth, Nlat, Nlon,
                                        Itime, Idepth, Ilat, Ilon,
                                        longitude, latitude, LAT_lb, LAT_ub,
                                        dAreas, scale, mask, &local_kernel);

                                ux_ux_pot.at(index) = uxux_tmp;
                                ux_uy_pot.at(index) = uxuy_tmp;
                                ux_uz_pot.at(index) = uxuz_tmp;
                                uy_uy_pot.at(index) = uyuy_tmp;
                                uy_uz_pot.at(index) = uyuz_tmp;
                                uz_uz_pot.at(index) = uzuz_tmp;

                                KE_pot_filt.at(index) = 0.5 * constants::rho0 * (uxux_tmp + uyuy_tmp + uzuz_tmp);

                                // tot
                                apply_filter_at_point_for_quadratics(
                                        uxux_tmp, uxuy_tmp, uxuz_tmp,
                                        uyuy_tmp, uyuz_tmp, uzuz_tmp,
                                        u_x_tot,  u_y_tot,  u_z_tot,
                                        Ntime, Ndepth, Nlat, Nlon,
                                        Itime, Idepth, Ilat, Ilon,
                                        longitude, latitude, LAT_lb, LAT_ub,
                                        dAreas, scale, mask, &local_kernel);

                                ux_ux_tot.at(index) = uxux_tmp;
                                ux_uy_tot.at(index) = uxuy_tmp;
                                ux_uz_tot.at(index) = uxuz_tmp;
                                uy_uy_tot.at(index) = uyuy_tmp;
                                uy_uz_tot.at(index) = uyuz_tmp;
                                uz_uz_tot.at(index) = uzuz_tmp;

                                KE_tot_filt.at(index) = 0.5 * constants::rho0 * (uxux_tmp + uyuy_tmp + uzuz_tmp);


                            }  // end if(masked) block
                        }  // end for(depth) block
                    }  // end for(time) block
                    if (constants::DO_TIMING) { 
                        timing_records.add_to_record(MPI_Wtime() - clock_on, "filter_tor_pot");
                    }
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

        // Write to file
        if (not(constants::NO_FULL_OUTPUTS)) {
            // Don't mask these fields, since they are filled over land from the projection
            write_field_to_output(coarse_F_tor, "coarse_F_tor", starts, counts, fname, NULL);
            write_field_to_output(coarse_F_pot, "coarse_F_pot", starts, counts, fname, NULL);
        }

        // Get pot and tor velocities
        toroidal_vel_from_F( u_lon_tor, u_lat_tor, coarse_F_tor,
                longitude, latitude, Ntime, Ndepth, Nlat, Nlon, mask);

        potential_vel_from_F(u_lon_pot, u_lat_pot, coarse_F_pot,
                longitude, latitude, Ntime, Ndepth, Nlat, Nlon, mask);

        #pragma omp parallel \
        default( none ) \
        shared( mask, u_lon_tor, u_lat_tor, u_lon_pot, u_lat_pot, u_lon_tot, u_lat_tot ) \
        private( index )
        {
            #pragma omp for collapse(1) schedule(guided, OMP_chunksize)
            for (index = 0; index < u_lon_tor.size(); ++index) {
                if (mask.at(index)) {
                    u_lon_tot.at(index) = u_lon_tor.at(index) + u_lon_pot.at(index);
                    u_lat_tot.at(index) = u_lat_tor.at(index) + u_lat_pot.at(index);
                }
            }
        }

        if (not(constants::NO_FULL_OUTPUTS)) {
            write_field_to_output(u_lon_tor, "u_lon_tor", starts, counts, fname, &mask);
            write_field_to_output(u_lat_tor, "u_lat_tor", starts, counts, fname, &mask);

            write_field_to_output(u_lon_pot, "u_lon_pot", starts, counts, fname, &mask);
            write_field_to_output(u_lat_pot, "u_lat_pot", starts, counts, fname, &mask);
        }

        // We'll need vorticity, so go ahead and compute it
        compute_vorticity(
                vort_tor_r, null_vector, null_vector,
                div_tor, OkuboWeiss_tor,
                u_r_zero, u_lon_tor, u_lat_tor,
                Ntime, Ndepth, Nlat, Nlon,
                longitude, latitude, mask);

        compute_vorticity(
                vort_pot_r, null_vector, null_vector,
                div_pot, OkuboWeiss_pot,
                u_r_zero, u_lon_pot, u_lat_pot,
                Ntime, Ndepth, Nlat, Nlon,
                longitude, latitude, mask);

        compute_vorticity(
                vort_tot_r, null_vector, null_vector,
                div_tot, OkuboWeiss_tot,
                u_r_zero, u_lon_tot, u_lat_tot,
                Ntime, Ndepth, Nlat, Nlon,
                longitude, latitude, mask);

        if (not(constants::MINIMAL_OUTPUT)) {
            write_field_to_output(div_tor, "div_tor", starts, counts, fname, &mask);
            write_field_to_output(div_pot, "div_pot", starts, counts, fname, &mask);
            write_field_to_output(div_tot, "div_tot", starts, counts, fname, &mask);
        }

        #pragma omp parallel \
        default( none ) \
        shared( KE_tor_coarse, KE_tor_fine, KE_tor_filt, KE_tor_fine_mod, KE_tor_orig, \
                KE_pot_coarse, KE_pot_fine, KE_pot_filt, KE_pot_fine_mod, KE_pot_orig, \
                KE_tot_coarse, KE_tot_fine, KE_tot_filt, KE_tot_fine_mod, KE_tot_orig, \
                Enst_tor, Enst_pot, Enst_tot, mask, \
                u_lon_tor, u_lat_tor, u_lon_pot, u_lat_pot, u_lon_tot, u_lat_tot, \
                vort_tor_r, vort_pot_r, vort_tot_r ) \
        private( index )
        {
            #pragma omp for collapse(1) schedule(guided, OMP_chunksize)
            for (index = 0; index < u_lon_tor.size(); ++index) {
                if (mask.at(index)) {
                    KE_tor_coarse.at(index) = 0.5 * constants::rho0 * ( pow(u_lon_tor.at(index), 2.) + pow(u_lat_tor.at(index), 2.) );
                    KE_pot_coarse.at(index) = 0.5 * constants::rho0 * ( pow(u_lon_pot.at(index), 2.) + pow(u_lat_pot.at(index), 2.) );
                    KE_tot_coarse.at(index) = 0.5 * constants::rho0 * ( pow(u_lon_tot.at(index), 2.) + pow(u_lat_tot.at(index), 2.) );

                    KE_tor_fine.at(index) = KE_tor_filt.at(index) - KE_tor_coarse.at(index);
                    KE_pot_fine.at(index) = KE_pot_filt.at(index) - KE_pot_coarse.at(index);
                    KE_tot_fine.at(index) = KE_tot_filt.at(index) - KE_tot_coarse.at(index);

                    KE_tor_fine_mod.at(index) = KE_tor_orig.at(index) - KE_tor_coarse.at(index);
                    KE_pot_fine_mod.at(index) = KE_pot_orig.at(index) - KE_pot_coarse.at(index);
                    KE_tot_fine_mod.at(index) = KE_tot_orig.at(index) - KE_tot_coarse.at(index);

                    Enst_tor.at(index) = 0.5 * constants::rho0 * ( pow(vort_tor_r.at(index), 2.) );
                    Enst_pot.at(index) = 0.5 * constants::rho0 * ( pow(vort_pot_r.at(index), 2.) );
                    Enst_tot.at(index) = 0.5 * constants::rho0 * ( pow(vort_tot_r.at(index), 2.) );
                }
            }
        }

        if (not(constants::NO_FULL_OUTPUTS)) {
            write_field_to_output( KE_tor_filt, "KE_tor_filt", starts, counts, fname, &mask);
            write_field_to_output( KE_pot_filt, "KE_pot_filt", starts, counts, fname, &mask);
            write_field_to_output( KE_tot_filt, "KE_tot_filt", starts, counts, fname, &mask);

            write_field_to_output( KE_tor_fine, "KE_tor_fine", starts, counts, fname, &mask);
            write_field_to_output( KE_pot_fine, "KE_pot_fine", starts, counts, fname, &mask);
            write_field_to_output( KE_tot_fine, "KE_tot_fine", starts, counts, fname, &mask);
        }

        if (not(constants::MINIMAL_OUTPUT)) {
            write_field_to_output( KE_tor_fine_mod, "KE_tor_fine_mod", starts, counts, fname, &mask);
            write_field_to_output( KE_pot_fine_mod, "KE_pot_fine_mod", starts, counts, fname, &mask);
            write_field_to_output( KE_tot_fine_mod, "KE_tot_fine_mod", starts, counts, fname, &mask);
        }

        // Compute Pi
        vel_Spher_to_Cart( u_x_coarse, u_y_coarse, u_z_coarse, u_r_zero, u_lon_tor, u_lat_tor, mask, time, depth, latitude, longitude );
        compute_Pi( Pi_tor, u_x_coarse,  u_y_coarse,  u_z_coarse,
                    ux_ux_tor, ux_uy_tor, ux_uz_tor, uy_uy_tor, uy_uz_tor, uz_uz_tor,
                    Ntime, Ndepth, Nlat, Nlon, longitude, latitude, mask);

        vel_Spher_to_Cart( u_x_coarse, u_y_coarse, u_z_coarse, u_r_zero, u_lon_pot, u_lat_pot, mask, time, depth, latitude, longitude );
        compute_Pi( Pi_pot, u_x_coarse,  u_y_coarse,  u_z_coarse,
                    ux_ux_pot, ux_uy_pot, ux_uz_pot, uy_uy_pot, uy_uz_pot, uz_uz_pot,
                    Ntime, Ndepth, Nlat, Nlon, longitude, latitude, mask);

        vel_Spher_to_Cart( u_x_coarse, u_y_coarse, u_z_coarse, u_r_zero, u_lon_tot, u_lat_tot, mask, time, depth, latitude, longitude );
        compute_Pi( Pi_tot, u_x_coarse,  u_y_coarse,  u_z_coarse,
                    ux_ux_tot, ux_uy_tot, ux_uz_tot, uy_uy_tot, uy_uz_tot, uz_uz_tot,
                    Ntime, Ndepth, Nlat, Nlon, longitude, latitude, mask);

        if (not(constants::MINIMAL_OUTPUT)) {
            write_field_to_output( Pi_tor, "Pi_tor", starts, counts, fname, &mask);
            write_field_to_output( Pi_pot, "Pi_pot", starts, counts, fname, &mask);
            write_field_to_output( Pi_tot, "Pi_tot", starts, counts, fname, &mask);
        }

        //
        //// on-line postprocessing, if desired
        //

        if (constants::APPLY_POSTPROCESS) {
            MPI_Barrier(MPI_COMM_WORLD);
            if (constants::DO_TIMING) { clock_on = MPI_Wtime(); }

            if (wRank == 0) { fprintf(stdout, "Beginning post-process routines\n"); }
            fflush(stdout);

            Apply_Postprocess_Routines(
                    postprocess_fields_tor, postprocess_names, OkuboWeiss_tor,
                    time, depth, latitude, longitude, mask, dAreas,
                    myCounts, myStarts, scales.at(Iscale), "postprocess_toroidal");

            Apply_Postprocess_Routines(
                    postprocess_fields_pot, postprocess_names, OkuboWeiss_pot,
                    time, depth, latitude, longitude, mask, dAreas,
                    myCounts, myStarts, scales.at(Iscale), "postprocess_potential");

            Apply_Postprocess_Routines(
                    postprocess_fields_tot, postprocess_names, OkuboWeiss_tot,
                    time, depth, latitude, longitude, mask, dAreas,
                    myCounts, myStarts, scales.at(Iscale), "postprocess_full");

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
