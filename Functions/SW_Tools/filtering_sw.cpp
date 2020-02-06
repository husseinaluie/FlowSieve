#include <math.h>
#include <algorithm>
#include <vector>
#include <omp.h>
#include <mpi.h>
#include "../../functions.hpp"
#include "../../netcdf_io.hpp"
#include "../../constants.hpp"
#include "../../differentiation_tools.hpp"

void filtering_sw(
        const std::vector<double> & full_u,     /**< [in] Full u_lon velocity array */
        const std::vector<double> & full_v,     /**< [in] Full u_lat velocity array */
        const std::vector<double> & full_h,     /**< [in] Full pressure array */
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

    const unsigned int num_pts = Ntime * Ndepth * Nlat * Nlon;
    const int num_pts2 = Ntime * Ndepth * Nlat * Nlon;
    char fname [50];
    
    const int ndims = 4;
    size_t starts[ndims] = {
        size_t(myStarts.at(0)), size_t(myStarts.at(1)), 
        size_t(myStarts.at(2)), size_t(myStarts.at(3))};
    size_t counts[ndims] = {
        size_t(Ntime), size_t(Ndepth), 
        size_t(Nlat), size_t(Nlon)};
    std::vector<std::string> vars_to_write;

    int LAT_lb, LAT_ub, Itime, Idepth, Ilat, Ilon, index, tid;
    size_t tmp_ind;

    #if DEBUG >= 1
    if (wRank == 0) { fprintf(stdout, "Computing derivatives.\n"); }
    #endif

    std::vector<double> local_kernel(Nlat * Nlon);

    std::vector<double> u_x(num_pts), u_y(num_pts), 
        v_x(num_pts), v_y(num_pts), full_KE(num_pts);

    // Compute quadratics
    std::vector<double> uh(num_pts), vh(num_pts), hh(num_pts), 
        huu(num_pts), huv(num_pts), hvv(num_pts);
    #pragma omp parallel \
    shared(full_u, full_v, full_h, uh, vh, hh, huu, huv, hvv) \
    private(Itime, Idepth, Ilat, Ilon, index)
    {
        #pragma omp for collapse(1) schedule(static)
        for (index = 0; index < num_pts2; index++) {
            uh.at(index)  = full_u.at(index) * full_h.at(index);
            vh.at(index)  = full_v.at(index) * full_h.at(index);
            hh.at(index)  = full_h.at(index) * full_h.at(index);

            huu.at(index) = uh.at(index) * full_u.at(index);
            huv.at(index) = uh.at(index) * full_v.at(index);
            hvv.at(index) = vh.at(index) * full_v.at(index);
        }
    }

    // Now prepare to filter
    double scale;
    double h_tmp, uh_tmp, vh_tmp, hh_tmp, huu_tmp, huv_tmp, hvv_tmp;

    double u_c, v_c, h_c;
    std::vector<double> 
        coarse_u(num_pts), coarse_v(num_pts), coarse_h(num_pts), 
        coarse_hh(num_pts), KE(num_pts), PE(num_pts),
        uKE(num_pts), vKE(num_pts), uPE(num_pts), vPE(num_pts),
        coarse_uu(num_pts), coarse_uv(num_pts), coarse_vv(num_pts),
        Pi(num_pts), Gamma(num_pts), KE_trans(num_pts), PE_trans(num_pts), KE2PE(num_pts),
        tau_uu(num_pts), tau_uv(num_pts), tau_vv(num_pts), tau_hp(num_pts),
        hui_tau_uiu(num_pts), hui_tau_uiv(num_pts);

    vars_to_write.push_back("full_u");
    vars_to_write.push_back("full_v");
    vars_to_write.push_back("full_h");

    vars_to_write.push_back("coarse_u");
    vars_to_write.push_back("coarse_v");
    vars_to_write.push_back("coarse_h");

    vars_to_write.push_back("Pi");
    vars_to_write.push_back("Gamma");
    vars_to_write.push_back("KE_trans");
    vars_to_write.push_back("PE_trans");
    vars_to_write.push_back("KE2PE");

    vars_to_write.push_back("KE");
    vars_to_write.push_back("PE");

    int perc_base = 5;
    int perc, perc_count=0;

    // Set up filtering vectors
    std::vector<double*> filtered_vals;
    std::vector<bool> filt_use_mask;
    std::vector<const std::vector<double>*> filter_fields;


    std::vector<const std::vector<double>*> deriv_fields;
    std::vector<double*> x_deriv_vals, y_deriv_vals, z_deriv_vals;
    double coarse_u_x, coarse_u_y, coarse_v_x, coarse_v_y, coarse_h_x, coarse_h_y,
           hui_tau_uiu_x, hui_tau_uiv_y, tau_hp_x, tau_hp_y, uKE_x, uPE_x, vKE_y, vPE_y;

    filter_fields.push_back(&full_h); filt_use_mask.push_back(true);

    filter_fields.push_back(&uh);     filt_use_mask.push_back(true);
    filter_fields.push_back(&vh);     filt_use_mask.push_back(true);
    filter_fields.push_back(&hh);     filt_use_mask.push_back(true);

    filter_fields.push_back(&huu);    filt_use_mask.push_back(true);
    filter_fields.push_back(&huv);    filt_use_mask.push_back(true);
    filter_fields.push_back(&hvv);    filt_use_mask.push_back(true);

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
                longitude, latitude, dAreas, scale,\
                full_h, uh, vh, hh, huu, huv, hvv,\
                coarse_u, coarse_v, coarse_h, coarse_hh,\
                coarse_uu, coarse_uv, coarse_vv, \
                perc_base, filter_fields, filt_use_mask, PE, KE)\
        private(Itime, Idepth, Ilat, Ilon, index,\
                h_tmp, uh_tmp, vh_tmp, hh_tmp, u_c, v_c, \
                huu_tmp, huv_tmp, hvv_tmp, \
                LAT_lb, LAT_ub, tid, filtered_vals) \
        firstprivate(perc, wRank, local_kernel, perc_count)
        {

            filtered_vals.resize(0);

            filtered_vals.push_back(&h_tmp);

            filtered_vals.push_back(&uh_tmp);
            filtered_vals.push_back(&vh_tmp);
            filtered_vals.push_back(&hh_tmp);

            filtered_vals.push_back(&huu_tmp);
            filtered_vals.push_back(&huv_tmp);
            filtered_vals.push_back(&hvv_tmp);

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

                            coarse_h.at(index) = h_tmp;

                            coarse_u.at(index) = uh_tmp / h_tmp;
                            coarse_v.at(index) = vh_tmp / h_tmp;

                            coarse_hh.at(index) = hh_tmp;

                            coarse_uu.at(index) = huu_tmp / h_tmp;
                            coarse_uv.at(index) = huv_tmp / h_tmp;
                            coarse_vv.at(index) = hvv_tmp / h_tmp;

                            u_c = coarse_u.at(index);
                            v_c = coarse_v.at(index);

                            // KE = 0.5 * h_c * (u_c**2 + v_c**2)
                            KE.at(index) = 0.5 * h_tmp * ( pow(u_c, 2) + pow(v_c, 2) );
                            
                            // PE = 0.5 * g0 * h_c * h_c
                            PE.at(index) = 0.5 * constants::g * pow(h_tmp, 2);


                        }  // end for(depth) block
                    }  // end for(time) block
                }  // end for(longitude) block
            }  // end for(latitude) block
        }  // end pragma parallel block
        #if DEBUG >= 0
        if (wRank == 0) { fprintf(stdout, "\n"); }
        #endif

        // Compute some taus in advance (since we'll need derivatives)
        #pragma omp parallel \
        shared(tau_uu, tau_uv, tau_vv, tau_hp, \
                hui_tau_uiu, hui_tau_uiv, uKE, vKE, uPE, vPE,\
                coarse_u, coarse_v, coarse_h, coarse_hh,\
                coarse_uu, coarse_uv, coarse_vv)\
        private(index, u_c, v_c, h_c)
        {
            #pragma omp for collapse(1) schedule(static)
            for (index = 0; index < num_pts2; index++) {

                u_c = coarse_u.at(index);
                v_c = coarse_v.at(index);
                h_c = coarse_h.at(index);

                tau_uu.at(index) = coarse_uu.at(index) - u_c * u_c;
                tau_uv.at(index) = coarse_uv.at(index) - u_c * v_c;
                tau_vv.at(index) = coarse_vv.at(index) - v_c * v_c;

                hui_tau_uiu.at(index) =   h_c * u_c * tau_uu.at(index)
                                        + h_c * v_c * tau_uv.at(index);
                hui_tau_uiv.at(index) =   h_c * u_c * tau_uv.at(index)
                                        + h_c * v_c * tau_vv.at(index);

                tau_hp.at(index) = constants::g * (coarse_hh.at(index) - h_c * h_c);

                uKE.at(index) = u_c * KE.at(index);
                vKE.at(index) = v_c * KE.at(index);

                uPE.at(index) = u_c * PE.at(index);
                vPE.at(index) = v_c * PE.at(index);
            }
        }

        deriv_fields.clear();
        deriv_fields.push_back(&coarse_u);
        deriv_fields.push_back(&coarse_v);
        deriv_fields.push_back(&coarse_h);
        deriv_fields.push_back(&hui_tau_uiu);
        deriv_fields.push_back(&hui_tau_uiv);
        deriv_fields.push_back(&tau_hp);
        deriv_fields.push_back(&uKE);
        deriv_fields.push_back(&vKE);
        deriv_fields.push_back(&uPE);
        deriv_fields.push_back(&vPE);

        // Now go ahead and compute the various transfer fields
        #pragma omp parallel \
        shared(coarse_u, coarse_v, coarse_h, hui_tau_uiu, hui_tau_uiv,\
                tau_hp, uKE, uPE, vKE, vPE, deriv_fields,\
                tau_uu, tau_uv, tau_vv,\
                Pi, KE_trans, Gamma, PE_trans, KE2PE) \
        private(index, Itime, Idepth, Ilat, Ilon, tmp_ind, \
                x_deriv_vals, y_deriv_vals, z_deriv_vals, \
                coarse_u_x, coarse_u_y, coarse_v_x, coarse_v_y,\
                coarse_h_x, coarse_h_y, hui_tau_uiu_x, hui_tau_uiv_y,\
                tau_hp_x, tau_hp_y, uKE_x, uPE_x, vKE_y, vPE_y)
        {

            x_deriv_vals.clear();
            x_deriv_vals.push_back(&coarse_u_x);
            x_deriv_vals.push_back(&coarse_v_x);
            x_deriv_vals.push_back(&coarse_h_x);
            x_deriv_vals.push_back(&hui_tau_uiu_x);
            x_deriv_vals.push_back(NULL);
            x_deriv_vals.push_back(&tau_hp_x);
            x_deriv_vals.push_back(&uKE_x);
            x_deriv_vals.push_back(NULL);
            x_deriv_vals.push_back(&uPE_x);
            x_deriv_vals.push_back(NULL);

            y_deriv_vals.clear();
            y_deriv_vals.push_back(&coarse_u_y);
            y_deriv_vals.push_back(&coarse_v_y);
            y_deriv_vals.push_back(&coarse_h_y);
            y_deriv_vals.push_back(NULL);
            y_deriv_vals.push_back(&hui_tau_uiv_y);
            y_deriv_vals.push_back(&tau_hp_y);
            y_deriv_vals.push_back(NULL);
            y_deriv_vals.push_back(&vKE_y);
            y_deriv_vals.push_back(NULL);
            y_deriv_vals.push_back(&vPE_y);

            z_deriv_vals.clear();
            for (tmp_ind = 0; tmp_ind < y_deriv_vals.size(); ++tmp_ind) {
                z_deriv_vals.push_back(NULL);
            }

            #pragma omp for collapse(1) schedule(static)
            for (index = 0; index < num_pts2; index++) {

                // Convert our four-index to a one-index
                Index1to4(index, Itime, Idepth, Ilat, Ilon,
                                 Ntime, Ndepth, Nlat, Nlon);

                // We'll need some derivatives
                Cart_derivatives_at_point(
                        x_deriv_vals, y_deriv_vals,
                        z_deriv_vals, deriv_fields,
                        latitude, longitude,
                        Itime, Idepth, Ilat, Ilon,
                        Ntime, Ndepth, Nlat, Nlon,
                        mask, 1);

                //
                //// KE
                //

                // Pi
                // h u_{i,j} tau(u_i,u_j)  - 0.5 *  u_i * tau(h, p)_{,i}
                Pi.at(index) = coarse_h.at(index) * (
                                 coarse_u_x * tau_uu.at(index)
                               + coarse_u_y * tau_uv.at(index)
                               + coarse_v_x * tau_uv.at(index)
                               + coarse_v_y * tau_vv.at(index)
                        )
                    - 0.5 * coarse_u.at(index) * tau_hp_x
                    - 0.5 * coarse_v.at(index) * tau_hp_y;

                // Transport
                // -( u_j * KE + h u_i tau(u_i,u_j) )_,j
                KE_trans.at(index) = - ( uKE_x + vKE_y + hui_tau_uiu_x + hui_tau_uiv_y );

                //
                //// PE
                //

                // Gamma
                // 0
                Gamma.at(index) = 0.;

                // Transport
                // -( u_j * PE )_,j
                PE_trans.at(index) = - ( uPE_x + vPE_y );

                //
                //// Conversion
                //

                // h * u_i * p_{,i} 
                KE2PE.at(index) = constants::g * coarse_h.at(index) * (
                            coarse_u.at(index) * coarse_h_x
                          + coarse_v.at(index) * coarse_h_y
                        );
            }
        }

        #if DEBUG >= 2
        fprintf(stdout, "  = Rank %d finished filtering loop =\n", wRank);
        fflush(stdout);
        #endif

        // Write to file
        write_field_to_output(full_u, "full_u", starts, counts, fname, &mask);
        write_field_to_output(full_v, "full_v", starts, counts, fname, &mask);
        write_field_to_output(full_h, "full_h", starts, counts, fname, &mask);

        write_field_to_output(coarse_u, "coarse_u", starts, counts, fname, &mask);
        write_field_to_output(coarse_v, "coarse_v", starts, counts, fname, &mask);
        write_field_to_output(coarse_h, "coarse_h", starts, counts, fname, &mask);

        write_field_to_output(Pi,       "Pi",       starts, counts, fname, &mask);
        write_field_to_output(Gamma,    "Gamma",    starts, counts, fname, &mask);
        write_field_to_output(KE_trans, "KE_trans", starts, counts, fname, &mask);
        write_field_to_output(PE_trans, "PE_trans", starts, counts, fname, &mask);
        write_field_to_output(KE2PE,    "KE2PE",    starts, counts, fname, &mask);

        write_field_to_output(KE, "KE", starts, counts, fname, &mask);
        write_field_to_output(PE, "PE", starts, counts, fname, &mask);

        #if DEBUG >= 0
        // Flushing stdout is necessary for SLURM outputs.
        fflush(stdout);
        #endif

    }  // end for(scale) block
} // end filtering
