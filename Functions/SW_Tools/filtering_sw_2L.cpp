#include <math.h>
#include <algorithm>
#include <vector>
#include <omp.h>
#include <mpi.h>
#include "../../functions.hpp"
#include "../../functions_sw.hpp"
#include "../../netcdf_io.hpp"
#include "../../constants.hpp"
#include "../../differentiation_tools.hpp"

void filtering_sw_2L(
        const std::vector<double> & full_u,     /**< [in] Full u_lon velocity array */
        const std::vector<double> & full_v,     /**< [in] Full u_lat velocity array */
        const std::vector<double> & full_h,     /**< [in] Full pressure array */
        const std::vector<double> & rhos,
        const double nu,
        const std::vector<double> & scales,     /**< [in] Array of filtering scales */
        const std::vector<double> & dAreas,     /**< [in] Array of cell areas (2D) (compute_areas()) */
        const std::vector<double> & time,       /**< [in] Time dimension (1D) */
        const std::vector<double> & depth,      /**< [in] Depth dimension (1D) */
        const std::vector<double> & longitude,  /**< [in] Longitude dimension (1D) */
        const std::vector<double> & latitude,   /**< [in] Latitude dimension (1D) */
        const std::vector<double> & mask,       /**< [in] Array to distinguish between land and water cells (2D) */
        const std::vector<int>    & myCounts,   /**< [in] Array of dimension sizes */
        const std::vector<int>    & myStarts,   /**< [in] Array of dimension sizes */
        const MPI_Comm comm                     /**< [in] MPI Communicator */
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

    const unsigned int num_pts  = Ntime * Ndepth * Nlat * Nlon;
    const int          num_pts2 = Ntime * Ndepth * Nlat * Nlon;
    char fname [50];
    
    // Some IO parameters for writing
    const int ndims = 4;
    size_t starts[ndims] = {
        size_t(myStarts.at(0)), size_t(myStarts.at(1)), 
        size_t(myStarts.at(2)), size_t(myStarts.at(3))};
    size_t counts[ndims] = {
        size_t(Ntime), size_t(Ndepth), 
        size_t(Nlat), size_t(Nlon)};
    std::vector<std::string> vars_to_write;

    // Various indices that we'll use
    int LAT_lb, LAT_ub, Itime, Idepth, Ilat, Ilon, index, ind_UL, ind_LL, tid;

    // Compute the kernal alpha value (for baroclinic transfers)
    const double kern_alpha = kernel_alpha();

    std::vector<double> local_kernel(Nlat * Nlon);

    std::vector<double> u_x(num_pts), u_y(num_pts), 
        v_x(num_pts), v_y(num_pts), full_KE(num_pts);

    // Compute pressure in each layer
    #if DEBUG >= 1
    if (wRank == 0) { fprintf(stdout, "Computing pressure.\n"); }
    #endif
    std::vector<double> full_p(num_pts);
    Compute_pressure_2L(full_p, full_h, rhos,
            Ntime, Ndepth, Nlat, Nlon);

    // We also need the pressure gradients (for one of the transfer terms)
    #if DEBUG >= 1
    if (wRank == 0) { fprintf(stdout, "Computing pressure gradient.\n"); }
    #endif
    std::vector<double> full_p_x(num_pts), full_p_y(num_pts);
    Compute_gradient(full_p_x, full_p_y, full_p, latitude, longitude,
            Ntime, Ndepth, Nlat, Nlon, mask);

    // We'll need the laplacian of our velocity fields to track viscosity
    //   since we want tilde(lap(u)), we will return h * lap(u)
    #if DEBUG >= 1
    if (wRank == 0) { fprintf(stdout, "Computing Laplacians.\n"); }
    #endif
    std::vector<double> h_lap_u(num_pts), h_lap_v(num_pts), lap_h(num_pts);
    Compute_Laplacians(h_lap_u, h_lap_v, lap_h, full_u, full_v, full_h,
                latitude, longitude, Ntime, Ndepth, Nlat, Nlon, mask);

    // Compute product terms (need for tildes and taus)
    #if DEBUG >= 1
    if (wRank == 0) { fprintf(stdout, "Prepping quadratics for taus etc.\n"); }
    #endif
    std::vector<double> 
        uh(num_pts), vh(num_pts), hh(num_pts), 
        hp_x(num_pts), hp_y(num_pts),
        huu(num_pts), huv(num_pts), hvv(num_pts);
    #pragma omp parallel \
    default(none) \
    shared(full_u, full_v, full_h, uh, vh, hh, huu, huv, hvv,\
            hp_x, hp_y, full_p_x, full_p_y) \
    private(index)
    {
        #pragma omp for collapse(1) schedule(static)
        for (index = 0; index < num_pts2; index++) {
            uh.at(index)  = full_u.at(index) * full_h.at(index);
            vh.at(index)  = full_v.at(index) * full_h.at(index);
            hh.at(index)  = full_h.at(index) * full_h.at(index);

            hp_x.at(index)  = full_h.at(index) * full_p_x.at(index);
            hp_y.at(index)  = full_h.at(index) * full_p_y.at(index);

            huu.at(index) = uh.at(index) * full_u.at(index);
            huv.at(index) = uh.at(index) * full_v.at(index);
            hvv.at(index) = vh.at(index) * full_v.at(index);
        }
    }

    // 
    std::vector<double> temp_array(num_pts);

    // Now prepare to filter
    double scale;
    double u_tmp, v_tmp, h_tmp, p_tmp, uh_tmp, vh_tmp, 
           hpx_tmp, hpy_tmp, px_tmp, py_tmp,
           huu_tmp, huv_tmp, hvv_tmp, 
           h_lap_u_tmp, h_lap_v_tmp, lap_h_tmp;

    double u_c, v_c;
    std::vector<double> 
        u_tilde(num_pts), v_tilde(num_pts), u_bar(num_pts), v_bar(num_pts), 
        h_bar(num_pts), p_bar(num_pts),
        uu_tilde(num_pts), uv_tilde(num_pts), vv_tilde(num_pts),
        PE(num_pts), KE(num_pts), omega_bar(num_pts),
        tau_tilde_uiuj(num_pts), tau_bar_uih(num_pts), tau_bar_hp_i(num_pts),
        tau_hpx(num_pts), tau_hpy(num_pts),
        lap_u_tilde(num_pts), lap_v_tilde(num_pts), lap_h_bar(num_pts);

    vars_to_write.push_back("full_u");
    vars_to_write.push_back("full_v");
    vars_to_write.push_back("full_h");

    vars_to_write.push_back("u_bar");
    vars_to_write.push_back("v_bar");
    vars_to_write.push_back("h_bar");
    vars_to_write.push_back("p_bar");

    vars_to_write.push_back("u_tilde");
    vars_to_write.push_back("v_tilde");

    vars_to_write.push_back("KE");
    vars_to_write.push_back("PE");

    vars_to_write.push_back("Lambda_m");
    vars_to_write.push_back("KE2PE");
    vars_to_write.push_back("Pi");
    vars_to_write.push_back("Gamma");

    vars_to_write.push_back("PE_transport_by_u");
    vars_to_write.push_back("PE_transport_by_misc");
    vars_to_write.push_back("KE_transport_by_u");
    vars_to_write.push_back("KE_transport_by_smallscales");

    vars_to_write.push_back("misc1");
    vars_to_write.push_back("misc2");
    vars_to_write.push_back("misc3");
    vars_to_write.push_back("misc_conversion");

    vars_to_write.push_back("viscous_loss_KE");
    vars_to_write.push_back("viscous_loss_PE");

    vars_to_write.push_back("KE_true_bc");
    vars_to_write.push_back("PE_true_bc");
    vars_to_write.push_back("PE_true_bc_parts");

    int perc_base = 5;
    int perc, perc_count=0;

    // Set up filtering vectors
    std::vector<double*> filtered_vals;
    std::vector<bool> filt_use_mask;
    std::vector<const std::vector<double>*> filter_fields;

    filter_fields.push_back(&full_u);   filt_use_mask.push_back(true);
    filter_fields.push_back(&full_v);   filt_use_mask.push_back(true);
    filter_fields.push_back(&full_h);   filt_use_mask.push_back(true);
    filter_fields.push_back(&full_p);   filt_use_mask.push_back(true);

    filter_fields.push_back(&uh);       filt_use_mask.push_back(true);
    filter_fields.push_back(&vh);       filt_use_mask.push_back(true);

    filter_fields.push_back(&full_p_x); filt_use_mask.push_back(true);
    filter_fields.push_back(&full_p_y); filt_use_mask.push_back(true);

    filter_fields.push_back(&hp_x);     filt_use_mask.push_back(true);
    filter_fields.push_back(&hp_y);     filt_use_mask.push_back(true);

    filter_fields.push_back(&huu);      filt_use_mask.push_back(true);
    filter_fields.push_back(&huv);      filt_use_mask.push_back(true);
    filter_fields.push_back(&hvv);      filt_use_mask.push_back(true);

    filter_fields.push_back(&h_lap_u);  filt_use_mask.push_back(true);
    filter_fields.push_back(&h_lap_v);  filt_use_mask.push_back(true);
    filter_fields.push_back(  &lap_h);  filt_use_mask.push_back(true);

    //
    //// Begin the main filtering loop
    //
    #if DEBUG>=1
    if (wRank == 0) { fprintf(stdout, "Beginning main filtering loop.\n\n"); }
    #endif
    for (int Iscale = 0; Iscale < Nscales; Iscale++) {

        // Create the output file
        snprintf(fname, 50, "filter_%.6gkm.nc", scales.at(Iscale)/1e3);
        initialize_output_file(time, depth, longitude, latitude, 
                mask, vars_to_write, fname, scales.at(Iscale));

        // Add some attributes to the file
        add_attr_to_file("kernel_alpha", 
                kern_alpha * pow(scales.at(Iscale), 2), 
                fname);

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
        shared(mask, stdout,\
                longitude, latitude, dAreas, scale, \
                full_h, uh, vh, hh, huu, huv, hvv, \
                h_lap_u, h_lap_v, lap_h, \
                u_bar, v_bar, h_bar, p_bar, \
                u_tilde, v_tilde, uu_tilde, uv_tilde, vv_tilde, \
                lap_u_tilde, lap_v_tilde, lap_h_bar, \
                full_p_x, full_p_y, hp_x, hp_y, tau_hpx, tau_hpy, \
                perc_base, filter_fields, filt_use_mask)\
        private(Itime, Idepth, Ilat, Ilon, index,\
                u_tmp, v_tmp, h_tmp, p_tmp, uh_tmp, vh_tmp,\
                huu_tmp, huv_tmp, hvv_tmp, \
                px_tmp, py_tmp, hpx_tmp, hpy_tmp, \
                h_lap_u_tmp, h_lap_v_tmp, lap_h_tmp, \
                LAT_lb, LAT_ub, tid, filtered_vals) \
        firstprivate(perc, wRank, local_kernel, perc_count)
        {

            filtered_vals.resize(0);

            filtered_vals.push_back(&u_tmp);
            filtered_vals.push_back(&v_tmp);
            filtered_vals.push_back(&h_tmp);
            filtered_vals.push_back(&p_tmp);

            filtered_vals.push_back(&uh_tmp);
            filtered_vals.push_back(&vh_tmp);

            filtered_vals.push_back(&px_tmp);
            filtered_vals.push_back(&py_tmp);

            filtered_vals.push_back(&hpx_tmp);
            filtered_vals.push_back(&hpy_tmp);

            filtered_vals.push_back(&huu_tmp);
            filtered_vals.push_back(&huv_tmp);
            filtered_vals.push_back(&hvv_tmp);

            filtered_vals.push_back(&h_lap_u_tmp);
            filtered_vals.push_back(&h_lap_v_tmp);
            filtered_vals.push_back(  &lap_h_tmp);

            #pragma omp for collapse(2) schedule(static)
            for (Ilat = 0; Ilat < Nlat; Ilat++) {
                for (Ilon = 0; Ilon < Nlon; Ilon++) {

                    get_lat_bounds(LAT_lb, LAT_ub, latitude, Ilat, scale); 

                    #if DEBUG >= 0
                    tid = omp_get_thread_num();
                    if ( (tid == 0) and (wRank == 0) ) {
                        // Every perc_base percent, print a dot, but only the first thread
                        if ( ((double)(Ilat*Nlon + Ilon+1) / (Nlon*Nlat)) * 100 >= perc ) {
                            perc_count++;
                            if (perc_count % 5 == 0) { fprintf(stdout, "|"); }
                            else                     { fprintf(stdout, "."); }
                            fflush(stdout);
                            perc += perc_base;
                        }
                    }
                    #endif

                    std::fill(local_kernel.begin(), local_kernel.end(), 0);
                    compute_local_kernel(
                            local_kernel, scale, longitude, latitude,
                            Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon,
                            LAT_lb, LAT_ub);

                    for (Itime = 0; Itime < Ntime; Itime++) {
                        for (Idepth = 0; Idepth < Ndepth; Idepth++) {

                            // Convert our four-index to a one-index
                            index = Index(Itime, Idepth, Ilat, Ilon,
                                          Ntime, Ndepth, Nlat, Nlon);
    
                            // Apply the filter at the point
                            #if DEBUG >= 3
                            if (wRank == 0) { 
                                fprintf(stdout, "    Line %d of %s\n", 
                                        __LINE__, __FILE__); 
                            }
                            fflush(stdout);
                            #endif

                            // Filter desired fields
                            apply_filter_at_point(
                                    filtered_vals, filter_fields,
                                    Ntime,  Ndepth, Nlat, Nlon,
                                    Itime,  Idepth, Ilat, Ilon,
                                    longitude, latitude, LAT_lb, LAT_ub,
                                    dAreas, scale, mask, filt_use_mask,
                                    &local_kernel);

                            u_bar.at(index) = u_tmp;
                            v_bar.at(index) = v_tmp;
                            h_bar.at(index) = h_tmp;
                            p_bar.at(index) = p_tmp;

                            u_tilde.at(index) = uh_tmp / h_tmp;
                            v_tilde.at(index) = vh_tmp / h_tmp;

                            uu_tilde.at(index) = huu_tmp / h_tmp;
                            uv_tilde.at(index) = huv_tmp / h_tmp;
                            vv_tilde.at(index) = hvv_tmp / h_tmp;

                            lap_u_tilde.at(index) = h_lap_u_tmp / h_tmp;
                            lap_v_tilde.at(index) = h_lap_v_tmp / h_tmp;

                            lap_h_bar.at(index)   =   lap_h_tmp;

                            tau_hpx.at(index) = hpx_tmp - h_tmp*px_tmp;
                            tau_hpy.at(index) = hpy_tmp - h_tmp*py_tmp;

                        }  // end for(depth) block
                    }  // end for(time) block
                }  // end for(longitude) block
            }  // end for(latitude) block
        }  // end pragma parallel block
        #if DEBUG >= 0
        if (wRank == 0) { fprintf(stdout, "\n"); }
        #endif

        // Write to file so far
        write_field_to_output(full_u, "full_u", starts, counts, fname, &mask);
        write_field_to_output(full_v, "full_v", starts, counts, fname, &mask);
        write_field_to_output(full_h, "full_h", starts, counts, fname, &mask);

        write_field_to_output(u_bar, "u_bar", starts, counts, fname, &mask);
        write_field_to_output(v_bar, "v_bar", starts, counts, fname, &mask);
        write_field_to_output(h_bar, "h_bar", starts, counts, fname, &mask);
        write_field_to_output(p_bar, "p_bar", starts, counts, fname, &mask);

        write_field_to_output(u_tilde, "u_tilde", starts, counts, fname, &mask);
        write_field_to_output(v_tilde, "v_tilde", starts, counts, fname, &mask);

        #pragma omp parallel \
        default(none) \
        shared(u_tilde, v_tilde, h_bar, KE, PE, rhos) \
        private(index, Itime, Idepth, Ilat, Ilon, ind_UL, ind_LL, u_c, v_c)
        {
            #pragma omp for collapse(1) schedule(static)
            for (index = 0; index < num_pts2; index++) {

                Index1to4(index, Itime, Idepth, Ilat, Ilon,
                                 Ntime, Ndepth, Nlat, Nlon);

                ind_UL = Index(Itime, 0,      Ilat, Ilon, 
                               Ntime, Ndepth, Nlat, Nlon);
                ind_LL = Index(Itime, 1,      Ilat, Ilon, 
                               Ntime, Ndepth, Nlat, Nlon);

                // KE = 0.5 * bar(h) * (tilde(u)**2 + tilde(v)**2)
                u_c = u_tilde.at(index);
                v_c = v_tilde.at(index);
                KE.at(index) = 0.5 * rhos.at(Idepth) * h_bar.at(index) 
                                    * ( pow(u_c, 2) + pow(v_c, 2) );

                // PE
                if (Idepth == 0) {
                    PE.at(index) = 0.5 * rhos.at(Idepth) * constants::g * h_bar.at(ind_UL) 
                        * ( h_bar.at(ind_UL) + 2 * h_bar.at(ind_LL) );
                } else if (Idepth == 1) {
                    PE.at(index) = 0.5 * rhos.at(Idepth) * constants::g 
                        * pow(h_bar.at(index), 2);
                }
            }
        }

        write_field_to_output(KE, "KE", starts, counts, fname, &mask);
        write_field_to_output(PE, "PE", starts, counts, fname, &mask);

        // Viscous losses (KE)
        Compute_Viscous_Loss_KE(temp_array, h_bar, u_tilde, v_tilde, 
                lap_u_tilde, lap_v_tilde, nu, rhos,
                time, depth, latitude, longitude,
                Ntime, Ndepth, Nlat, Nlon);
        write_field_to_output(temp_array, "viscous_loss_KE", 
                starts, counts, fname, &mask);

        // Viscous losses (PE)
        Compute_Viscous_Loss_PE(temp_array, h_bar, 
                lap_h_bar, nu, rhos, time, depth, latitude, longitude,
                Ntime, Ndepth, Nlat, Nlon);
        write_field_to_output(temp_array, "viscous_loss_PE", 
                starts, counts, fname, &mask);

        // Lambda_m
        Compute_Baroclinic_Transfer_SW(temp_array, u_bar, v_bar, h_bar, p_bar,
                kern_alpha * pow(scales.at(Iscale), 2),
                time, depth, latitude, longitude, mask,
                Ntime, Ndepth, Nlat, Nlon);
        write_field_to_output(temp_array, "Lambda_m", 
                starts, counts, fname, &mask);

        // u_i,i * PE (traditional transfer)
        Compute_Traditional_Transfer_SW(temp_array, u_bar, v_bar, PE,
                time, depth, latitude, longitude, mask,
                Ntime, Ndepth, Nlat, Nlon);
        write_field_to_output(temp_array, "KE2PE", 
                starts, counts, fname, &mask);

        // Pi
        Compute_Pi_SW(temp_array, h_bar, u_tilde, v_tilde, 
                uu_tilde, uv_tilde, vv_tilde, rhos,
                time, depth, latitude, longitude, mask,
                Ntime, Ndepth, Nlat, Nlon);
        write_field_to_output(temp_array, "Pi", 
                starts, counts, fname, &mask);

        // Gamma
        Compute_Gamma_SW(temp_array, h_bar, u_bar, v_bar,
                u_tilde, v_tilde, rhos, 
                time, depth, latitude, longitude, mask,
                Ntime, Ndepth, Nlat, Nlon);
        write_field_to_output(temp_array, "Gamma", 
                starts, counts, fname, &mask);

        // Transport: PE by u
        Compute_Transport_SW(temp_array, u_bar, v_bar, PE,
                time, depth, latitude, longitude, mask,
                Ntime, Ndepth, Nlat, Nlon);
        write_field_to_output(temp_array, "PE_transport_by_u", 
                starts, counts, fname, &mask);

        // Transport: PE - other / mech. unknown
        Compute_Alt_PE_Transport(temp_array, u_tilde, v_tilde, h_bar, rhos,
                time, depth, latitude, longitude, mask,
                Ntime, Ndepth, Nlat, Nlon);
        write_field_to_output(temp_array, "PE_transport_by_misc", 
                starts, counts, fname, &mask);

        // Transport: KE by u
        Compute_Transport_SW(temp_array, u_tilde, v_tilde, KE,
                time, depth, latitude, longitude, mask,
                Ntime, Ndepth, Nlat, Nlon);
        write_field_to_output(temp_array, "KE_transport_by_u", 
                starts, counts, fname, &mask);

        // Transport: KE by small scales
        Compute_KE_Transport_smallscales_SW(temp_array, h_bar,
                u_tilde, v_tilde, uu_tilde, uv_tilde, vv_tilde, rhos,
                time, depth, latitude, longitude, mask,
                Ntime, Ndepth, Nlat, Nlon);
        write_field_to_output(temp_array, "KE_transport_by_smallscales", 
                starts, counts, fname, &mask);

        // Miscellany
        Compute_misc1_SW(temp_array, h_bar, u_bar, v_bar, rhos,
                time, depth, latitude, longitude, mask,
                Ntime, Ndepth, Nlat, Nlon);
        write_field_to_output(temp_array, "misc1", 
                starts, counts, fname, &mask);

        Compute_misc2_SW(temp_array, u_tilde, v_tilde, tau_hpx, tau_hpy,
                time, depth, latitude, longitude, mask,
                Ntime, Ndepth, Nlat, Nlon);
        write_field_to_output(temp_array, "misc2", 
                starts, counts, fname, &mask);

        Compute_misc3_SW(temp_array, h_bar, u_bar, v_bar, rhos,
                kern_alpha * pow(scales.at(Iscale), 2),
                time, depth, latitude, longitude, mask,
                Ntime, Ndepth, Nlat, Nlon);
        write_field_to_output(temp_array, "misc3", 
                starts, counts, fname, &mask);

        Compute_misc_conversion_SW(temp_array, h_bar, u_bar, v_bar, rhos,
                kern_alpha * pow(scales.at(Iscale), 2),
                time, depth, latitude, longitude, mask,
                Ntime, Ndepth, Nlat, Nlon);
        write_field_to_output(temp_array, "misc_conversion", 
                starts, counts, fname, &mask);

        // Full (unapproximated) 'baroclinic' term
        Compute_full_bc_term_KE(temp_array, h_bar, p_bar, u_tilde, u_bar,
                v_tilde, v_bar,
                time, depth, latitude, longitude, mask,
                Ntime, Ndepth, Nlat, Nlon);
        write_field_to_output(temp_array, "KE_true_bc", 
                starts, counts, fname, &mask);

        Compute_full_bc_term_PE(temp_array, h_bar, u_tilde, u_bar,
                v_tilde, v_bar, rhos,
                time, depth, latitude, longitude, mask,
                Ntime, Ndepth, Nlat, Nlon);
        write_field_to_output(temp_array, "PE_true_bc", 
                starts, counts, fname, &mask);

        Compute_full_bc_term_PE_parts(temp_array, h_bar, u_tilde, u_bar,
                v_tilde, v_bar, rhos,
                time, depth, latitude, longitude, mask,
                Ntime, Ndepth, Nlat, Nlon);
        write_field_to_output(temp_array, "PE_true_bc_parts", 
                starts, counts, fname, &mask);

        #if DEBUG >= 0
        // Flushing stdout is necessary for SLURM outputs.
        fflush(stdout);
        #endif

    }  // end for(scale) block
} // end filtering
