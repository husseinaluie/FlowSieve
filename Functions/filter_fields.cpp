#include <math.h>
#include <algorithm>
#include <vector>
#include <omp.h>
#include <mpi.h>
#include "../functions.hpp"
#include "../netcdf_io.hpp"
#include "../constants.hpp"
#include "../postprocess.hpp"

void filter_fields(
        const std::vector<const std::vector<double>*> & fields,
        const std::vector<std::string> var_names,
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
    
    // Set up some arrays to store the filtered fields
    const int num_fields = fields.size();
    std::vector<bool> filt_use_mask(num_fields, false);
    std::vector< std::vector<double> > filtered_fields(num_fields);
    for (int Ifield = 0; Ifield < num_fields; ++Ifield) {
        filtered_fields.at(Ifield).resize(num_pts);
    }
    std::vector<double*> filtered_vals_ptrs;
    std::vector<double>  filtered_vals;



    //
    int LAT_lb, LAT_ub,
        Itime, Idepth, Ilat, Ilon, mask_index,
        tid, index;

    std::vector<double> local_kernel(Nlat * Nlon);

    // Now prepare to filter
    double scale;

    int perc_base = 5;
    int perc, perc_count=0;

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
                    mask, var_names, fname, scales.at(Iscale));
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
        shared(mask, stdout, fields, filtered_fields, filt_use_mask, \
                longitude, latitude, dAreas, scale, perc_base,\
                clock_on, timing_records)\
        private(Itime, Idepth, Ilat, Ilon, index, mask_index,\
                tid, filtered_vals, filtered_vals_ptrs, LAT_lb, LAT_ub) \
        firstprivate(perc, wRank, local_kernel, perc_count)
        {

            filtered_vals.clear();
            filtered_vals.resize(num_fields);

            for (int Ifield = 0; Ifield < num_fields; ++Ifield) {
                filtered_vals_ptrs.push_back(&filtered_vals.at(Ifield));
            }

            #pragma omp for collapse(2) schedule(dynamic)
            for (Ilat = 0; Ilat < Nlat; Ilat++) {
                for (Ilon = 0; Ilon < Nlon; Ilon++) {

                    // Get the integration bounds (to speed things up)
                    get_lat_bounds(LAT_lb, LAT_ub, latitude,  Ilat, scale); 
                    mask_index = Index(0,     0,      Ilat, Ilon,
                                       Ntime, Ndepth, Nlat, Nlon);

                    // Every perc_base percent, print a dot, but only the first thread
                    #if DEBUG >= 0
                    tid = omp_get_thread_num();
                    if ( (tid == 0) and (wRank == 0) ) {
                        if ( ((double)(mask_index+1) / (Nlon*Nlat)) * 100 >= perc ) {
                            perc_count++;
                            if (perc_count % 5 == 0) { fprintf(stdout, "|"); }
                            else                     { fprintf(stdout, "."); }
                            fflush(stdout);
                            perc += perc_base;
                        }
                    }
                    #endif

                    // If not a land cell, then apply filter
                    if (mask.at(mask_index) == 1) {

                        // Pre-compute the kernel at this lat-lon point
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

                        // Loop through time/depth for this lat-lon points
                        for (Itime = 0; Itime < Ntime; Itime++) {
                            for (Idepth = 0; Idepth < Ndepth; Idepth++) {

                                // Convert our four-index to a one-index
                                index = Index(Itime, Idepth, Ilat, Ilon,
                                              Ntime, Ndepth, Nlat, Nlon);
    
                                // Apply the filter at the point
                                if (constants::DO_TIMING) { clock_on = MPI_Wtime(); }

                                apply_filter_at_point(
                                        filtered_vals_ptrs, fields,
                                        Ntime, Ndepth, Nlat, Nlon,
                                        Itime, Idepth, Ilat, Ilon,
                                        longitude, latitude, LAT_lb, LAT_ub,
                                        dAreas, scale, mask, filt_use_mask,
                                        &local_kernel);

                                for (int Ifield = 0; Ifield < num_fields; ++Ifield) {
                                    filtered_fields.at(Ifield).at(index) =
                                        filtered_vals.at(Ifield);
                                }
                            }
                        }
                    }  // end if(masked) block
                    else { // if not masked
                        if (constants::DO_TIMING) { clock_on = MPI_Wtime(); }
                        for (Itime = 0; Itime < Ntime; Itime++) {
                            for (Idepth = 0; Idepth < Ndepth; Idepth++) {

                                // Convert our four-index to a one-index
                                index = Index(Itime, Idepth, Ilat, Ilon,
                                              Ntime, Ndepth, Nlat, Nlon);

                                for (int Ifield = 0; Ifield < num_fields; ++Ifield) {
                                    filtered_fields.at(Ifield).at(index) =
                                        constants::fill_value;
                                }
                            }
                        }
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
        if (constants::DO_TIMING) { clock_on = MPI_Wtime(); }

        if (not(constants::NO_FULL_OUTPUTS)) {
            for (int Ifield = 0; Ifield < num_fields; ++Ifield) {
                write_field_to_output(filtered_fields.at(Ifield), 
                        var_names.at(Ifield).c_str(), starts, counts, fname, &mask);
            }
        }

        if (constants::DO_TIMING) { 
            timing_records.add_to_record(MPI_Wtime() - clock_on, "writing");
        }

    }  // end for(scale) block
} // end filtering
