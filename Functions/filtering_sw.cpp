#include <math.h>
#include <algorithm>
#include <vector>
#include <omp.h>
#include <mpi.h>
#include "../functions.hpp"
#include "../netcdf_io.hpp"
#include "../constants.hpp"
#include "../differentiation_tools.hpp"

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

    int LAT_lb, LAT_ub;

    #if DEBUG >= 1
    if (wRank == 0) { fprintf(stdout, "Computing derivatives.\n"); }
    #endif

    std::vector<double> local_kernel(Nlat * Nlon);

    std::vector<double> u_x(num_pts);
    std::vector<double> u_y(num_pts);

    std::vector<double> v_x(num_pts);
    std::vector<double> v_y(num_pts);

    std::vector<double> full_KE(num_pts);

    double ux, uy, vx, vy;
    double tau_uh_x, tau_vh_y, uKE_x, uPE_x, vKE_y, vPE_y;
    std::vector<double*> x_deriv_vals, y_deriv_vals, z_deriv_vals;
    std::vector<const std::vector<double>*> deriv_fields;

    deriv_fields.push_back(&full_u);
    deriv_fields.push_back(&full_v);

    // Compute derivatives
    int index, Itime, Idepth, Ilat, Ilon, tid;
    #pragma omp parallel \
    shared(full_u, full_v, deriv_fields, latitude, longitude, mask) \
    private(Itime, Idepth, Ilat, Ilon, index, \
            ux, vx, uy, vy,\
            x_deriv_vals, y_deriv_vals, z_deriv_vals)
    {
        x_deriv_vals.push_back(&ux);
        x_deriv_vals.push_back(&vx);

        y_deriv_vals.push_back(&uy);
        y_deriv_vals.push_back(&vy);

        z_deriv_vals.push_back(NULL);
        z_deriv_vals.push_back(NULL);

        #pragma omp for collapse(1) schedule(static)
        for (index = 0; index < num_pts2; index++) {

            // Convert our four-index to a one-index
            Index1to4(index, Itime, Idepth, Ilat, Ilon,
                             Ntime, Ndepth, Nlat, Nlon);

                Cart_derivatives_at_point(
                        x_deriv_vals, y_deriv_vals,
                        z_deriv_vals, deriv_fields,
                        latitude, longitude,
                        Itime, Idepth, Ilat, Ilon,
                        Ntime, Ndepth, Nlat, Nlon,
                        mask);

                u_x.at(index) = ux;
                u_y.at(index) = uy;
                v_x.at(index) = vx;
                v_y.at(index) = vy;

        } // done index loop
    } // done pragma
    #if DEBUG >= 1
    if (wRank == 0) { fprintf(stdout, "   ... done.\n"); }
    #endif

    // Compute quadratics
    std::vector<double> uh, vh, uux, uvx, vuy, vvy;

    uh.resize(num_pts);
    vh.resize(num_pts);

    uux.resize(num_pts);
    uvx.resize(num_pts);
    vuy.resize(num_pts);
    vvy.resize(num_pts);

    #pragma omp parallel \
    shared(full_u, full_v, full_h, uh, vh, uux, uvx, vuy, vvy) \
    private(Itime, Idepth, Ilat, Ilon, index)
    {
        #pragma omp for collapse(1) schedule(static)
        for (index = 0; index < num_pts2; index++) {
            uh.at(index)  = full_u.at(index) * full_h.at(index);
            vh.at(index)  = full_v.at(index) * full_h.at(index);

            uux.at(index) = full_u.at(index) * u_x.at(index);
            uvx.at(index) = full_u.at(index) * v_x.at(index);
            vuy.at(index) = full_v.at(index) * u_y.at(index);
            vvy.at(index) = full_v.at(index) * v_y.at(index);
        }
    }

    // Now prepare to filter
    double scale;
    double u_tmp, v_tmp, h_tmp, uh_tmp, vh_tmp, 
           ux_tmp, uy_tmp, vx_tmp, vy_tmp,
           uux_tmp, uvx_tmp, vuy_tmp, vvy_tmp;

    double u_c, v_c;
    std::vector<double> 
        coarse_u(num_pts), coarse_v(num_pts), coarse_h(num_pts), 
        tau_uh(num_pts), tau_vh(num_pts), 
        uKE(num_pts), vKE(num_pts), uPE(num_pts), vPE(num_pts),
        coarse_uh(num_pts), coarse_vh(num_pts), KE(num_pts), PE(num_pts),
        coarse_ux(num_pts), coarse_uy(num_pts), coarse_vx(num_pts), coarse_vy(num_pts),
        coarse_uux(num_pts), coarse_uvx(num_pts), coarse_vuy(num_pts), coarse_vvy(num_pts),
        Pi(num_pts), Gamma(num_pts), KE_trans(num_pts), PE_trans(num_pts), KE2PE(num_pts);

    vars_to_write.push_back("full_u");
    vars_to_write.push_back("full_v");
    vars_to_write.push_back("full_h");

    vars_to_write.push_back("coarse_u");
    vars_to_write.push_back("coarse_v");
    vars_to_write.push_back("coarse_h");

    vars_to_write.push_back("coarse_uh");
    vars_to_write.push_back("coarse_vh");

    vars_to_write.push_back("coarse_uux");
    vars_to_write.push_back("coarse_uvx");
    vars_to_write.push_back("coarse_vuy");
    vars_to_write.push_back("coarse_vvy");

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

    filter_fields.push_back(&full_u); filt_use_mask.push_back(true);
    filter_fields.push_back(&full_v); filt_use_mask.push_back(true);
    filter_fields.push_back(&full_h); filt_use_mask.push_back(true);

    filter_fields.push_back(&uh);     filt_use_mask.push_back(true);
    filter_fields.push_back(&vh);     filt_use_mask.push_back(true);

    filter_fields.push_back(&u_x);    filt_use_mask.push_back(true);
    filter_fields.push_back(&u_y);    filt_use_mask.push_back(true);
    filter_fields.push_back(&v_x);    filt_use_mask.push_back(true);
    filter_fields.push_back(&v_y);    filt_use_mask.push_back(true);

    filter_fields.push_back(&uux);    filt_use_mask.push_back(true);
    filter_fields.push_back(&uvx);    filt_use_mask.push_back(true);
    filter_fields.push_back(&vuy);    filt_use_mask.push_back(true);
    filter_fields.push_back(&vvy);    filt_use_mask.push_back(true);

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
                full_u, full_v, full_h, uh, vh, uux, uvx, vuy, vvy,\
                coarse_u, coarse_v, coarse_h, coarse_uh, coarse_vh,\
                coarse_ux, coarse_uy, coarse_vx, coarse_vy,\
                coarse_uux, coarse_uvx, coarse_vuy, coarse_vvy,\
                perc_base, filter_fields, filt_use_mask, PE, KE)\
        private(Itime, Idepth, Ilat, Ilon, index,\
                u_tmp, v_tmp, h_tmp, uh_tmp, vh_tmp,\
                ux_tmp, uy_tmp, vx_tmp, vy_tmp,\
                uux_tmp, uvx_tmp, vuy_tmp, vvy_tmp,\
                LAT_lb, LAT_ub, tid, filtered_vals) \
        firstprivate(perc, wRank, local_kernel, perc_count)
        {

            filtered_vals.resize(0);

            filtered_vals.push_back(&u_tmp);
            filtered_vals.push_back(&v_tmp);
            filtered_vals.push_back(&h_tmp);

            filtered_vals.push_back(&uh_tmp);
            filtered_vals.push_back(&vh_tmp);

            filtered_vals.push_back(&ux_tmp);
            filtered_vals.push_back(&uy_tmp);
            filtered_vals.push_back(&vx_tmp);
            filtered_vals.push_back(&vy_tmp);

            filtered_vals.push_back(&uux_tmp);
            filtered_vals.push_back(&uvx_tmp);
            filtered_vals.push_back(&vuy_tmp);
            filtered_vals.push_back(&vvy_tmp);

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

                    compute_local_kernel(
                            local_kernel, scale, longitude, latitude,
                            Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon);

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

                            coarse_u.at(  index) = u_tmp;
                            coarse_v.at(  index) = v_tmp;
                            coarse_h.at(  index) = h_tmp;

                            coarse_uh.at( index) = uh_tmp;
                            coarse_vh.at( index) = vh_tmp;

                            coarse_ux.at( index) = ux_tmp;
                            coarse_uy.at( index) = uy_tmp;
                            coarse_vx.at( index) = vx_tmp;
                            coarse_vy.at( index) = vy_tmp;

                            coarse_uux.at(index) = uux_tmp;
                            coarse_uvx.at(index) = uvx_tmp;
                            coarse_vuy.at(index) = vuy_tmp;
                            coarse_vvy.at(index) = vvy_tmp;

                            // KE = 0.5 * h_c * (u_c**2 + v_c**2)
                            KE.at(index) = 0.5 * h_tmp * ( pow(u_tmp, 2) + pow(v_tmp, 2) );
                            
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
        shared(tau_uh, tau_vh, coarse_uh, coarse_vh, coarse_u, coarse_v, coarse_h)\
        private(index, u_c, v_c)
        {
            #pragma omp for collapse(1) schedule(static)
            for (index = 0; index < num_pts2; index++) {

                u_c = coarse_u.at(index);
                v_c = coarse_v.at(index);

                tau_uh.at(index) = coarse_uh.at(index) - u_c * coarse_h.at(index);
                tau_vh.at(index) = coarse_vh.at(index) - v_c * coarse_h.at(index);

                uKE.at(index) = u_c * KE.at(index);
                vKE.at(index) = v_c * KE.at(index);

                uPE.at(index) = u_c * PE.at(index);
                vPE.at(index) = v_c * PE.at(index);

            }
        }

        deriv_fields.resize(6);
        deriv_fields.at(0) = &tau_uh;
        deriv_fields.at(1) = &tau_vh;
        deriv_fields.at(2) = &uKE;
        deriv_fields.at(3) = &vKE;
        deriv_fields.at(4) = &uPE;
        deriv_fields.at(5) = &vPE;

        // Now go ahead and compute the various transfer fields
        #pragma omp parallel \
        shared(full_u, full_v, full_h, uh, vh, uux, uvx, vuy, vvy, deriv_fields) \
        private(index, Itime, Idepth, Ilat, Ilon, u_c, v_c, \
                x_deriv_vals, y_deriv_vals, z_deriv_vals, \
                tau_uh_x, tau_vh_y, uKE_x, uPE_x, vKE_y, vPE_y)
        {

            x_deriv_vals.resize(6);
            x_deriv_vals.at(0) = &tau_uh_x;
            x_deriv_vals.at(1) = NULL;
            x_deriv_vals.at(2) = &uKE_x;
            x_deriv_vals.at(3) = NULL;
            x_deriv_vals.at(4) = &uPE_x;
            x_deriv_vals.at(5) = NULL;

            y_deriv_vals.resize(6);
            y_deriv_vals.at(0) = NULL;
            y_deriv_vals.at(1) = &tau_vh_y;
            y_deriv_vals.at(2) = NULL;
            y_deriv_vals.at(3) = &vKE_y;
            y_deriv_vals.at(4) = NULL;
            y_deriv_vals.at(5) = &vPE_y;

            z_deriv_vals.resize(6);
            z_deriv_vals.at(0) = NULL;
            z_deriv_vals.at(1) = NULL;
            z_deriv_vals.at(2) = NULL;
            z_deriv_vals.at(3) = NULL;
            z_deriv_vals.at(4) = NULL;
            z_deriv_vals.at(5) = NULL;

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
                        mask);

                u_c = coarse_u.at(index);
                v_c = coarse_v.at(index);

                //
                //// KE
                //

                // Pi
                // - h * u_j * tau(u_i, u_j,i)
                // - 0.5 * u_i * u_i * [ tau(u_j, h) ]_,j
                Pi.at(index) = 
                    - coarse_h.at(index) * ( 
                         u_c * (   coarse_uux.at(index) - u_c * coarse_ux.at(index)
                                 + coarse_vuy.at(index) - v_c * coarse_uy.at(index)
                               )
                       + v_c * (   coarse_uvx.at(index) - u_c * coarse_vx.at(index) 
                                 + coarse_vvy.at(index) - v_c * coarse_vy.at(index) 
                               )
                       );

                Pi.at(index) += - 0.5 * ( pow(u_c, 2) + pow(v_c, 2) ) 
                                      * ( tau_uh_x + tau_vh_y );


                // Transport
                // -( u_j * (KE + PE) )_,j
                KE_trans.at(index) = - ( uKE_x + uPE_x + vKE_y + vPE_y );

                //
                //// PE
                //

                // Gamma
                // - g * h * tau(u_i, h)_,i
                Gamma.at(index) = - constants::g * coarse_h.at(index) 
                                                 * ( tau_uh_x + tau_vh_y);

                // Transport
                // -( u_j * PE )_,j
                PE_trans.at(index) = - ( uPE_x + vPE_y );

                //
                //// Conversion
                //

                // PE * u_j,j
                KE2PE.at(index) = PE.at(index) * (  coarse_ux.at(index) 
                                                  + coarse_vy.at(index) );
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

        write_field_to_output(coarse_uh, "coarse_uh", starts, counts, fname, &mask);
        write_field_to_output(coarse_vh, "coarse_vh", starts, counts, fname, &mask);

        write_field_to_output(coarse_uux, "coarse_uux", starts, counts, fname, &mask);
        write_field_to_output(coarse_uvx, "coarse_uvx", starts, counts, fname, &mask);
        write_field_to_output(coarse_vuy, "coarse_vuy", starts, counts, fname, &mask);
        write_field_to_output(coarse_vvy, "coarse_vvy", starts, counts, fname, &mask);

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
