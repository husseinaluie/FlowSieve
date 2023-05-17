#include <math.h>
#include <algorithm>
#include <vector>
#include <omp.h>
#include <mpi.h>
#include <bitset>
#include <deque>
#include "../../functions.hpp"
#include "../../netcdf_io.hpp"
#include "../../constants.hpp"
#include "../../postprocess.hpp"
#include "../../preprocess.hpp"

/*!
 * \brief Main filtering driver for Helmholtz decomposed data
 *
 * This function is the main filtering driver. It sets up the appropriate
 * loop sequences, calls the other funcations (velocity conversions), and
 * calls the IO functionality.
 *
 * @param[in]   source_data     dataset class instance containing data (Psi, Phi, etc)
 * @param[in]   scales          scales at which to filter the data
 * @param[in]   comm            MPI communicator (default MPI_COMM_WORLD)
 *
 */
void LLC_filtering_helmholtz(
        const dataset & source_data,
        const std::vector<double> & scales,
        const MPI_Comm comm
        ) {

    // Get dimension sizes
    const int   Nscales = scales.size(),
                Ntime   = source_data.Ntime,    // this is the MPI-local Ntime, not the full Ntime
                Ndepth  = source_data.Ndepth,   // this is the MPI-local Ndepth, not the full Ndepth
                Nlatlon = source_data.latitude.size();
    const size_t num_pts = Ntime * Ndepth * Nlatlon;

    const std::vector<double> zero_vector( num_pts, 0. );

    // Create some tidy names for variables
    const std::vector<double>   &latitude   = source_data.latitude,
                                &longitude  = source_data.longitude;

    const std::vector<double>   &F_potential    = source_data.variables.at("F_potential"),
                                &F_toroidal     = source_data.variables.at("F_toroidal");

    const std::vector<bool> &mask = source_data.mask;

    const std::vector<int>  &myStarts = source_data.myStarts;

    // Get some MPI info
    int wRank, wSize;
    MPI_Comm_rank( comm, &wRank );
    MPI_Comm_size( comm, &wSize );

    #if DEBUG >= 2
    if (wRank == 0) { fprintf(stdout, "\nEntered filtering_helmholtz\n\n"); }
    #endif

    // If we've passed the DO_TIMING flag, then create some timing vars
    Timing_Records timing_records;
    double clock_on;

    #if DEBUG >= 1
    if (wRank == 0) { fprintf( stdout, "\nPreparing to apply %d filters to data with (MPI-local) sizes (%'d - %'d - %'d) \n", Nscales, Ntime, Ndepth, Nlatlon ); }
    #endif

    char fname [50];
    
    const int ndims = 4;
    size_t starts[ndims] = { size_t(myStarts.at(0)), size_t(myStarts.at(1)), size_t(myStarts.at(2)) };
    size_t counts[ndims] = { size_t(Ntime),          size_t(Ndepth),         size_t(Nlatlon),       };
    size_t index;
    std::vector<std::string> vars_to_write;

    int LAT_lb, LAT_ub;

    std::vector<double> null_vector(0);

    #if DEBUG >= 2
    if (wRank == 0) { fprintf(stdout, "\nInitializing storage arrays.\n"); }
    #endif
    std::vector<double> 
        // Arrays to store filtered Phi and Psi fields (potential, pseudo-potential)
        coarse_F_tor(   num_pts, 0. ),
        coarse_F_pot(   num_pts, 0. ),
        u_r_coarse(     num_pts, 0. ),

        // Coarse KE (computed from velocities)
        KE_tor_coarse(  num_pts, 0. ),
        KE_pot_coarse(  num_pts, 0. ),
        KE_tot_coarse(  num_pts, 0. ),

        //
        //// Spherical velocity components
        //
        zero_array( num_pts, 0.),

        // Spherical - radial velocity
        //  by definition this is potential-only since toroidal is incompressible
        u_r( num_pts, 0.),

        // Spherical - zonal velocities
        u_lon_tor( num_pts, 0. ),
        u_lon_pot( num_pts, 0. ),
        u_lon_tot( num_pts, 0. ),

        // Spherical - meridional velocities
        u_lat_tor( num_pts, 0. ),
        u_lat_pot( num_pts, 0. ),
        u_lat_tot( num_pts, 0. ),

        // Vorticity (only r component)
        vort_tor_r( num_pts, 0. ),
        vort_pot_r( num_pts, 0. ),
        vort_tot_r( num_pts, 0. );


    //
    //// Compute original (unfiltered) KE
    //
     
    #if DEBUG >= 2
    if (wRank == 0) { fprintf(stdout, "\nExtracting velocities from Phi and Psi\n"); }
    #endif
    // Get pot and tor velocities
    if (constants::DO_TIMING) { clock_on = MPI_Wtime(); }
    toroidal_vel_from_F(  u_lon_tor, u_lat_tor, F_toroidal,  source_data, mask );
    potential_vel_from_F( u_lon_pot, u_lat_pot, F_potential, source_data, mask );
    if (constants::DO_TIMING) { timing_records.add_to_record(MPI_Wtime() - clock_on, "compute velocities from F"); }

    #if DEBUG >= 2
    if (wRank == 0) { fprintf(stdout, "\nComputing KE of unfiltered velocities\n"); }
    #endif
    if (constants::DO_TIMING) { clock_on = MPI_Wtime(); }
    #pragma omp parallel \
    default( none ) \
    shared( u_r, u_lon_tor, u_lat_tor, u_lon_pot, u_lat_pot, u_lon_tot, u_lat_tot) \
    private( index )
    {
        #pragma omp for collapse(1) schedule(guided)
        for (index = 0; index < u_lon_tor.size(); ++index) {
            u_lon_tot.at(index) = u_lon_tor.at(index) + u_lon_pot.at(index);
            u_lat_tot.at(index) = u_lat_tor.at(index) + u_lat_pot.at(index);
        }
    }
    if (constants::DO_TIMING) { timing_records.add_to_record(MPI_Wtime() - clock_on, "compute KE and Enstrophy"); }

    #if DEBUG >= 2
    if (wRank == 0) { fprintf(stdout, "\nFlagging variables for output\n"); }
    #endif
    vars_to_write.push_back("coarse_F_tor");
    vars_to_write.push_back("coarse_F_pot");

    vars_to_write.push_back("u_lon_tor");
    vars_to_write.push_back("u_lat_tor");

    vars_to_write.push_back("u_lon_pot");
    vars_to_write.push_back("u_lat_pot");

    // Compute the kernal alpha value (for baroclinic transfers)
    const double kern_alpha = kernel_alpha();

    // Now prepare to filter
    double scale;
    int Itime, Idepth, Ilat, Ilon, Ilatlon, thread_id, num_threads, prev_Ilat = -1;

    int perc_base = 5;
    int perc, perc_count=0;

    //
    //// Set up filtering vectors
    //
    std::vector<double*> filtered_vals;
    std::vector<const std::vector<double>*> filter_fields;

    double F_pot_tmp;
    filter_fields.push_back(&F_potential);

    double F_tor_tmp;
    filter_fields.push_back(&F_toroidal);

    const size_t Nfields = filter_fields.size();
    size_t Ifield;


    // Filter-loop variables
    size_t target_index, search_index, neighbour_index;
    const size_t num_neighbours = source_data.num_neighbours;
    double target_lat, target_lon, kern_val, dArea, local_val, kernel_normalization,
           local_lat, local_lon, local_dist;
    
    std::deque<size_t> points_to_test;
    std::vector<bool> was_rejected(num_pts, false), was_accepted(num_pts, false);
    // Want to use bitset directly, but does not accept runtime-specified size
    // so will use std::vector<bool>, which is a bitset under the hood
    //std::bitset< num_pts > was_rejected, was_accepted;

    //
    //// Begin the main filtering loop
    //
    #if DEBUG>=1
    if (wRank == 0) { fprintf(stdout, "\nBeginning main filtering loop.\n\n"); }
    #endif
    for (int Iscale = 0; Iscale < Nscales; Iscale++) {

        // Rest our timing records
        timing_records.reset();

        // Create the output file
        snprintf(fname, 50, "filter_%.6gkm.nc", scales.at(Iscale)/1e3);
        if (not(constants::NO_FULL_OUTPUTS)) {
            initialize_output_file( source_data, vars_to_write, fname, scales.at(Iscale));
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
        shared( source_data, mask, stdout, perc_base, filter_fields, timing_records, clock_on, \
                longitude, latitude, scale, F_potential, F_toroidal, coarse_F_tor, coarse_F_pot, \
                scales, Iscale \
                ) \
        private(Itime, Idepth, index, Ilatlon, F_tor_tmp, F_pot_tmp, filtered_vals, \
                target_index, search_index, neighbour_index, kern_val, dArea, local_val, Ifield, \
                was_rejected, was_accepted, points_to_test, target_lat, target_lon, kernel_normalization, \
                local_lat, local_lon, local_dist \
                )\
        firstprivate(perc, wRank, perc_count, Nlatlon, Ndepth, Ntime, Nfields )
        {

            filtered_vals.clear();

            filtered_vals.push_back(&F_pot_tmp);
            filtered_vals.push_back(&F_tor_tmp);

            was_rejected.resize( num_pts );
            was_accepted.resize( num_pts );

            #pragma omp for collapse(1) schedule(guided)
            for (target_index = 0; target_index < Nlatlon; target_index++ ) {

                std::fill(was_rejected.begin(), was_rejected.end(), false);
                std::fill(was_accepted.begin(), was_accepted.end(), false);

                target_lat = latitude.at(  target_index );
                target_lon = longitude.at( target_index );


                kernel_normalization = 0.0;
                dArea = source_data.areas.at(target_index);
                kernel_normalization += 1. * dArea;
                for (Ifield = 0; Ifield < Nfields; Ifield++) {
                    local_val = filter_fields.at(Ifield)->at(target_index);
                    *(filtered_vals.at(Ifield)) = dArea * local_val;
                    // not +=, here we re-set the values to the local value
                }

                // Next, seed the 'points to test' with the adjacent points
                points_to_test.clear();
                for (neighbour_index = 0; neighbour_index < num_neighbours; neighbour_index++ ) {
                    points_to_test.push_back( source_data.adjacency_indices.at(target_index)[neighbour_index] );
                }

                // So long as we still have points to test, keep testing!
                // This implements a breadth-first search through the adjacency matrix to build
                //      filtering kernel. If a point is within distance, add it to the kernel,
                //      and then test it's neighbours. If those are in, test their neighbours,
                //      and so on. If a point is too far away, we do not test it's neighbours.
                // This then assumes that for any point X within distance L of Y, that there
                //      is an adjacency path from Y to X strictly using points within distance
                //      L of Y.
                // Along the way, test for double inclusing to make sure that we don't test 
                //      points repeatedly etc.
                while ( points_to_test.size() > 0 ) {

                    // Pull out the most-recently-added point, and remove it from the 'to test' list,
                    // since we're testing it now.
                    // Since we're pulling out the most-recently-added, that effectively makes this a
                    // depth-first search to build the kernel.
                    search_index = points_to_test.front();
                    points_to_test.pop_front();

                    local_lat = source_data.latitude.at(  search_index );
                    local_lon = source_data.longitude.at( search_index );
                    local_dist = distance( target_lon, target_lat, local_lon, local_lat );

                    kern_val = kernel( local_dist, scales.at(Iscale) );

                    if ( ( kern_val > 1e-10 ) or ( local_dist <= 0.5 * scales.at(Iscale) )  ) {
                        was_accepted[search_index] = true;

                        // Accumulate the coarse values
                        dArea = source_data.areas.at(search_index);
                        kernel_normalization += kern_val * dArea;
                        for (Ifield = 0; Ifield < Nfields; Ifield++) {
                            local_val = filter_fields.at(Ifield)->at(search_index);
                            *(filtered_vals.at(Ifield)) += kern_val * dArea * local_val;
                        }

                        // Since this point is in the kernel
                        // add its neighbours to the 'to search' list
                        for ( int ii = 0; ii < num_neighbours; ii++ ) {
                            neighbour_index = source_data.adjacency_indices.at(search_index)[ii];

                            // but first check if that neighbour has been rejected already
                            if (was_rejected[neighbour_index]) { continue; }

                            // then check if that neighbour is already accepted
                            if (was_accepted[neighbour_index]) { continue; }

                            // then check if that neighbour is already on the search list
                            if ( std::find( points_to_test.begin(),
                                            points_to_test.end(),
                                            neighbour_index )
                                    != points_to_test.end()
                               ) { continue; }

                            // if it's not on any of those list already
                            // then add it to the 'points to test'
                            points_to_test.push_back( neighbour_index );

                        }
                    } else {
                        // Otherwise, record this point as rejected, and move on
                        was_rejected[search_index] = true;
                    }

                } // loop through to make the kernel

                coarse_F_pot.at(target_index) = (kernel_normalization == 0) ? 0 : F_pot_tmp / kernel_normalization;
                coarse_F_tor.at(target_index) = (kernel_normalization == 0) ? 0 : F_tor_tmp / kernel_normalization;

            }  // end pragma parallel block
        }
        #if DEBUG >= 0
        if (wRank == 0) { fprintf(stdout, "\n"); }
        #endif

        #if DEBUG >= 2
        fprintf(stdout, "  = Rank %d finished filtering loop =\n", wRank);
        fflush(stdout);
        #endif

        // Write to file
        if (constants::DO_TIMING) { clock_on = MPI_Wtime(); }
        write_field_to_output(coarse_F_tor, "coarse_F_tor", starts, counts, fname, NULL);
        write_field_to_output(coarse_F_pot, "coarse_F_pot", starts, counts, fname, NULL);
        if (constants::DO_TIMING) { timing_records.add_to_record(MPI_Wtime() - clock_on, "writing");  }

        // Get pot and tor velocities
        if (constants::DO_TIMING) { clock_on = MPI_Wtime(); }
        toroidal_vel_from_F( u_lon_tor, u_lat_tor, coarse_F_tor, source_data, mask);
        potential_vel_from_F(u_lon_pot, u_lat_pot, coarse_F_pot, source_data, mask);
        if (constants::DO_TIMING) { timing_records.add_to_record(MPI_Wtime() - clock_on, "compute velocities from F"); }

        if (constants::DO_TIMING) { clock_on = MPI_Wtime(); }
        write_field_to_output(u_lon_tor, "u_lon_tor", starts, counts, fname, &mask);
        write_field_to_output(u_lat_tor, "u_lat_tor", starts, counts, fname, &mask);

        write_field_to_output(u_lon_pot, "u_lon_pot", starts, counts, fname, &mask);
        write_field_to_output(u_lat_pot, "u_lat_pot", starts, counts, fname, &mask);
        if (constants::DO_TIMING) { timing_records.add_to_record(MPI_Wtime() - clock_on, "writing");  }

        // If we're doing timings, then print out and reset values now
        if (constants::DO_TIMING) { 
            timing_records.print();
            timing_records.reset();
            fflush(stdout);
        }

    }  // end for(scale) block
} // end filtering
