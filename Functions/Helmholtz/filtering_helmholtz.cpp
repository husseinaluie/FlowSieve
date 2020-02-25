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
        const std::vector<double> & mask,
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

    std::vector<double> local_kernel(Nlat * Nlon);

    std::vector<double> 
        coarse_F_tor(num_pts, 0.), 
        coarse_F_pot(num_pts, 0.),
        KE_tor(      num_pts, 0.),
        KE_pot(      num_pts, 0.),
        KE_all(      num_pts, 0.),
        div_tor(     num_pts, 0.),
        div_pot(     num_pts, 0.),
        u_x_tmp(     num_pts, 0.),
        u_y_tmp(     num_pts, 0.),
        u_z_tmp(     num_pts, 0.),
        u_r_zero(    num_pts, 0.),
        u_lon_tor(   num_pts, 0.),
        u_lat_tor(   num_pts, 0.),
        u_lon_pot(   num_pts, 0.),
        u_lat_pot(   num_pts, 0.);

    vars_to_write.push_back("coarse_F_tor");
    vars_to_write.push_back("coarse_F_pot");
    if (not(constants::MINIMAL_OUTPUT)) {
        vars_to_write.push_back("u_lon_tor");
        vars_to_write.push_back("u_lat_tor");

        vars_to_write.push_back("u_lon_pot");
        vars_to_write.push_back("u_lat_pot");

        vars_to_write.push_back("div_tor");
        vars_to_write.push_back("div_pot");
    }

    // Now prepare to filter
    double scale, F_tor_tmp, F_pot_tmp;
    int Itime, Idepth, Ilat, Ilon, index, mask_index, tid;
    size_t index_st;

    int perc_base = 5;
    int perc, perc_count=0;

    // Set up filtering vectors
    std::vector<double*> filtered_vals;
    std::vector<bool> filt_use_mask;
    std::vector<const std::vector<double>*> filter_fields;

    filter_fields.push_back(&F_potential);
    filt_use_mask.push_back(false);

    filter_fields.push_back(&F_toroidal);
    filt_use_mask.push_back(false);

    // Preset some post-processing variables
    std::vector<const std::vector<double>*> postprocess_fields;
    std::vector<std::string> postprocess_names;

    postprocess_names.push_back( "F_toroidal");
    postprocess_fields.push_back(&coarse_F_tor);

    postprocess_names.push_back( "F_potential");
    postprocess_fields.push_back(&coarse_F_pot);

    postprocess_names.push_back( "KE_tor");
    postprocess_fields.push_back(&KE_tor);

    postprocess_names.push_back( "KE_pot");
    postprocess_fields.push_back(&KE_pot);

    postprocess_names.push_back( "KE_all");
    postprocess_fields.push_back(&KE_all);

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
        shared(mask, stdout, perc_base, \
                filter_fields, filt_use_mask, \
                timing_records, clock_on, \
                longitude, latitude, dAreas, scale, \
                F_potential, F_toroidal, coarse_F_tor, coarse_F_pot) \
        private(Itime, Idepth, Ilat, Ilon, index, mask_index, \
                F_tor_tmp, F_pot_tmp, \
                LAT_lb, LAT_ub, tid, filtered_vals) \
        firstprivate(perc, wRank, local_kernel, perc_count)
        {

            filtered_vals.clear();

            filtered_vals.push_back(&F_pot_tmp);
            filtered_vals.push_back(&F_tor_tmp);

            #pragma omp for collapse(2) schedule(guided)
            for (Ilat = 0; Ilat < Nlat; Ilat++) {
                for (Ilon = 0; Ilon < Nlon; Ilon++) {

                    get_lat_bounds(LAT_lb, LAT_ub, latitude,  Ilat, scale); 
                    mask_index = Index(0,     0,      Ilat, Ilon,
                                       Ntime, Ndepth, Nlat, Nlon);

                    #if DEBUG >= 0
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

                                coarse_F_pot.at(index) = F_pot_tmp;
                                coarse_F_tor.at(index) = F_tor_tmp;

                                if (constants::DO_TIMING) { 
                                    timing_records.add_to_record(MPI_Wtime() - clock_on,
                                            "filter_tor_pot");
                                }

                            }  // end for(depth) block
                        }  // end for(time) block
                    }  // end if(masked) block
                    else { // if not masked
                        if (constants::DO_TIMING) { clock_on = MPI_Wtime(); }
                        for (Itime = 0; Itime < Ntime; Itime++) {
                            for (Idepth = 0; Idepth < Ndepth; Idepth++) {

                                // Convert our four-index to a one-index
                                index = Index(Itime, Idepth, Ilat, Ilon,
                                              Ntime, Ndepth, Nlat, Nlon);

                                coarse_F_pot.at(index) = constants::fill_value;
                                coarse_F_tor.at(index) = constants::fill_value;

                            }  // end for(depth) block
                        }  // end for(time) block
                        if (constants::DO_TIMING) { 
                            timing_records.add_to_record(MPI_Wtime() - clock_on, "land");
                        }
                    }  // end not(masked) block
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
            write_field_to_output(coarse_F_tor, "coarse_F_tor", 
                    starts, counts, fname, &mask);
            write_field_to_output(coarse_F_pot, "coarse_F_pot", 
                    starts, counts, fname, &mask);
        }

        // Get pot and tor velocities
        toroidal_vel_from_F( u_lon_tor, u_lat_tor, coarse_F_tor,
                longitude, latitude, Ntime, Ndepth, Nlat, Nlon, mask);

        potential_vel_from_F(u_lon_pot, u_lat_pot, coarse_F_pot,
                longitude, latitude, Ntime, Ndepth, Nlat, Nlon, mask);

        // Compute the divergence of the tor field
        vel_Spher_to_Cart(u_x_tmp,  u_y_tmp,   u_z_tmp,
                          u_r_zero, u_lon_tor, u_lat_tor,
                          mask, time, depth, latitude, longitude);
        compute_div_vel(div_tor, u_x_tmp, u_y_tmp, u_z_tmp, 
                longitude, latitude, Ntime, Ndepth, Nlat, Nlon, mask);

        // Compute the divergence of the tor field
        vel_Spher_to_Cart(u_x_tmp,  u_y_tmp,   u_z_tmp,
                          u_r_zero, u_lon_pot, u_lat_pot,
                          mask, time, depth, latitude, longitude);
        compute_div_vel(div_pot, u_x_tmp, u_y_tmp, u_z_tmp, 
                longitude, latitude, Ntime, Ndepth, Nlat, Nlon, mask);


        if (not(constants::MINIMAL_OUTPUT)) {
            write_field_to_output(u_lon_tor, "u_lon_tor", 
                    starts, counts, fname, &mask);
            write_field_to_output(u_lat_tor, "u_lat_tor", 
                    starts, counts, fname, &mask);
            
            write_field_to_output(u_lon_pot, "u_lon_pot", 
                    starts, counts, fname, &mask);
            write_field_to_output(u_lat_pot, "u_lat_pot", 
                    starts, counts, fname, &mask);

            write_field_to_output(div_tor, "div_tor", 
                    starts, counts, fname, &mask);
            write_field_to_output(div_pot, "div_pot", 
                    starts, counts, fname, &mask);
        }


        #pragma omp parallel \
        default( none ) \
        shared( KE_tor, KE_pot, KE_all, mask, \
                u_lon_tor, u_lat_tor, u_lon_pot, u_lat_pot) \
        private( index_st )
        {
            #pragma omp for collapse(1) schedule(guided)
            for (index_st = 0; index_st < u_lon_tor.size(); ++index_st) {
                if (mask.at(index_st) == 1) {
                    KE_tor.at(index_st) =    pow(u_lon_tor.at(index_st), 2.) 
                                           + pow(u_lat_tor.at(index_st), 2.);

                    KE_pot.at(index_st) =    pow(u_lon_pot.at(index_st), 2.) 
                                           + pow(u_lat_pot.at(index_st), 2.);

                    KE_all.at(index_st) =    pow(   u_lon_pot.at(index_st) 
                                                  + u_lon_tor.at(index_st) , 2.) 
                                           + pow(   u_lat_pot.at(index_st) 
                                                  + u_lat_tor.at(index_st) , 2.);

                    KE_tor.at(index_st) *= 0.5 * constants::rho0;
                    KE_pot.at(index_st) *= 0.5 * constants::rho0;
                    KE_all.at(index_st) *= 0.5 * constants::rho0;

                } else {
                    KE_tor.at(index_st) = constants::fill_value;
                    KE_pot.at(index_st) = constants::fill_value;
                    KE_all.at(index_st) = constants::fill_value;
                }
            }
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
                    postprocess_fields, postprocess_names,
                    time, depth, latitude, longitude,
                    mask, dAreas,
                    myCounts, myStarts,
                    scales.at(Iscale));

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
