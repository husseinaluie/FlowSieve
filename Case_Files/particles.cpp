#include <fenv.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <algorithm>
#include <math.h>
#include <vector>
#include <mpi.h>
#include <omp.h>
#include <cassert>

#include "../netcdf_io.hpp"
#include "../functions.hpp"
#include "../constants.hpp"
#include "../particles.hpp"

int main(int argc, char *argv[]) {
    
    // PERIODIC_Y implies UNIFORM_LAT_GRID
    static_assert( (constants::UNIFORM_LAT_GRID) or (not(constants::PERIODIC_Y)),
            "PERIODIC_Y requires UNIFORM_LAT_GRID.\n"
            "Please update constants.hpp accordingly.\n");

    static_assert( ( (constants::PERIODIC_X) and (not(constants::PERIODIC_Y)) ),
            "The particles routine currently requires globe-like periodicity.\n"
            "Please update constants.hpp accordingly.\n");


    // Specify the number of OpenMP threads
    //   and initialize the MPI world
    int thread_safety_provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &thread_safety_provided);
    //MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI::ERRORS_THROW_EXCEPTIONS);

    int wRank=-1, wSize=-1;
    MPI_Comm_rank( MPI_COMM_WORLD, &wRank );
    MPI_Comm_size( MPI_COMM_WORLD, &wSize );

    //
    //// Parse command-line arguments
    //
    InputParser input(argc, argv);
    if(input.cmdOptionExists("--version")){
        if (wRank == 0) { print_compile_info(NULL); } 
        return 0;
    }

    // first argument is the flag, second argument is default value (for when flag is not present)
    const std::string &zonal_vel_name   = input.getCmdOption("--zonal_vel",   "uo");
    const std::string &merid_vel_name   = input.getCmdOption("--merid_vel",   "vo");
    const std::string &input_fname      = input.getCmdOption("--input_file",  "input.nc");
    const std::string &output_fname     = input.getCmdOption("--output_file", "particles.nc");

    // Set OpenMP thread number
    const int max_threads = omp_get_max_threads();
    omp_set_num_threads( max_threads );

    // Print some header info, depending on debug level
    print_header_info();

    std::vector<double> longitude, latitude, time, depth;
    std::vector<double> u_lon, u_lat;
    std::vector<double> mask;
    size_t II;

    // Read in source data / get size information
    #if DEBUG >= 1
    if (wRank == 0) { fprintf(stdout, "Reading in source data.\n\n"); }
    #endif

    // Read in the grid coordinates
    read_var_from_file(longitude, "longitude", input_fname);
    read_var_from_file(latitude,  "latitude",  input_fname);
    read_var_from_file(time,      "time",      input_fname);
    read_var_from_file(depth,     "depth",     input_fname);
     
    convert_coordinates(longitude, latitude);

    // Read in the velocity fields
    read_var_from_file(u_lon, zonal_vel_name, input_fname, &mask, NULL, NULL, false);
    read_var_from_file(u_lat, merid_vel_name, input_fname, &mask, NULL, NULL, false);

    // A sanity check to make sure that the loaded data is nonzero
    double u_max = 0., v_max = 0.;
    size_t land_count = 0, water_count = 0;
    #pragma omp parallel \
    default(none) shared( u_lon, u_lat, mask ) private( II )\
    reduction(max : u_max, v_max) reduction(+ : land_count, water_count)
    {
        #pragma omp for collapse(1) schedule(static)
        for ( II = 0; II < u_lon.size(); ++II ){
            if ( mask.at(II) == 1 ) {
                u_max = std::max( u_max, std::fabs(u_lon.at(II)) );
                v_max = std::max( v_max, std::fabs(u_lat.at(II)) );
                water_count++;
            } else {
                land_count++;
            }
        }
    }
    assert( u_max > 0 );
    assert( v_max > 0 );

    // Set the output times
    const double start_time = time.front(),
                 final_time = time.at(time.size()-2);
    const double out_freq = 10. * 60.;
    const int Nouts = (int) ( (final_time - start_time) / out_freq );

    std::vector<double> target_times(Nouts);
    for ( II = 0; II < target_times.size(); ++II ) {
        target_times.at(II) = start_time 
            + (double)II * ( final_time - start_time ) / Nouts;
    }

    // Get particle positions
    const int Npts = 1000;  // particles per MPI process
    std::vector<double> starting_lat(Npts), starting_lon(Npts);
    particles_initial_positions(starting_lat, starting_lon, Npts,
            latitude, longitude, mask);

    // Trajectories dimension (essentially just a numbering)
    std::vector<double> trajectories(Npts);
    for (II = 0; II < Npts; ++II) { trajectories.at(II) = (double) II + wRank * Npts; };

    // List the fields to track along particle trajectories
    std::vector<const std::vector<double>*> fields_to_track;
    std::vector<std::string> names_of_tracked_fields;

    names_of_tracked_fields.push_back( "vel_lon");
    fields_to_track.push_back(&u_lon);

    names_of_tracked_fields.push_back( "vel_lat");
    fields_to_track.push_back(&u_lat);

    // Storage for tracked fields
    std::vector< std::vector< double > >     field_trajectories(fields_to_track.size()),
                                         rev_field_trajectories(fields_to_track.size());

    for (size_t Ifield = 0; Ifield < fields_to_track.size(); ++Ifield) {
            field_trajectories.at(Ifield).resize(Npts * Nouts, constants::fill_value);
        rev_field_trajectories.at(Ifield).resize(Npts * Nouts, constants::fill_value);
    }

    // Initialize particle output file
    initialize_particle_file(target_times, trajectories, names_of_tracked_fields, output_fname);

    std::vector<double>     part_lon_hist(Npts * Nouts, constants::fill_value), 
                            part_lat_hist(Npts * Nouts, constants::fill_value),
                        rev_part_lon_hist(Npts * Nouts, constants::fill_value),
                        rev_part_lat_hist(Npts * Nouts, constants::fill_value),
                        trajectory_dists( Npts * Nouts, constants::fill_value);

    for (II = 0; II < Npts; ++II) {
        part_lon_hist.at(II) = starting_lon.at(II);
        part_lat_hist.at(II) = starting_lat.at(II);
    }

    size_t starts[2], counts[2];
    starts[0] = 0;
    counts[0] = Nouts;

    starts[1] = wRank * Npts;
    counts[1] = Npts;

    /*
    // Convert to seconds
    for (II = 0; II < time.size();         ++II) { time.at(II)         *= 60 * 60; }
    for (II = 0; II < target_times.size(); ++II) { target_times.at(II) *= 60 * 60; }
    */

    // Now do the particle routine
    particles_evolve_trajectories(
            part_lon_hist,     part_lat_hist,
        rev_part_lon_hist, rev_part_lat_hist,
            field_trajectories,
        rev_field_trajectories,
        starting_lat,  starting_lon,
        target_times,
        u_lon, u_lat,
        fields_to_track, names_of_tracked_fields,
        time, latitude, longitude,        
        mask);

    fprintf(stdout, "\nProcessor %d of %d finished stepping particles.\n", wRank+1, wSize);

    std::vector<double> out_mask(part_lon_hist.size());
    for (II = 0; II < out_mask.size(); ++II) {
        out_mask.at(II) = part_lon_hist.at(II) == constants::fill_value ? 0 : 1;
    }

    particles_fore_back_difference(trajectory_dists,     part_lon_hist,     part_lat_hist, 
                                                     rev_part_lon_hist, rev_part_lat_hist);

    MPI_Barrier(MPI_COMM_WORLD);
    write_field_to_output(part_lon_hist, "longitude", starts, counts, output_fname, &out_mask);
    write_field_to_output(part_lat_hist, "latitude",  starts, counts, output_fname, &out_mask);

    #if DEBUG >= 1
    write_field_to_output(rev_part_lon_hist, "rev_longitude", starts, counts, output_fname, &out_mask);
    write_field_to_output(rev_part_lat_hist, "rev_latitude",  starts, counts, output_fname, &out_mask);
    #endif

    write_field_to_output(trajectory_dists, "fore_back_dists", starts, counts, output_fname, &out_mask);

    for (size_t Ifield = 0; Ifield < fields_to_track.size(); ++Ifield) {
        MPI_Barrier(MPI_COMM_WORLD);
        write_field_to_output(field_trajectories.at(Ifield), 
                names_of_tracked_fields.at(Ifield),  
                starts, counts, output_fname, &out_mask);
    }

    MPI_Finalize();
    return 0;
}
