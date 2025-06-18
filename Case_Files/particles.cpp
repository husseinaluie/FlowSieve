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

    const bool asked_help = input.cmdOptionExists("--help");
    if (asked_help) {
        fprintf( stdout, "\033[1;4mThe command-line input arguments [and default values] are:\033[0m\n" );
    }
    // first argument is the flag, second argument is default value (for when flag is not present)
    const std::string   &zonal_vel_name   = input.getCmdOption("--zonal_vel",   
                                                               "uo",
                                                               asked_help,
                                                               "Name of zonal velocity variable in the input file."),
                        &merid_vel_name   = input.getCmdOption("--merid_vel",   
                                                               "vo",
                                                               asked_help,
                                                               "Name of the meridional velocity variable in the input file."),
                        &input_fname      = input.getCmdOption("--input_file",  
                                                               "input.nc",
                                                               asked_help,
                                                               "Path to netCDF file containing velocity information."),
                        &output_fname     = input.getCmdOption("--output_file", 
                                                               "particles.nc",
                                                               asked_help,
                                                               "Path to netCDF file to be created for output. Will clobber existing file."),
                        &time_units       = input.getCmdOption("--time_unit",   
                                                               "hours",
                                                               asked_help,
                                                               "Time unit used in the input file. Options are: seconds, minutes, hours, days");

    // particles per MPI process
    const std::string &particles_string = input.getCmdOption("--particle_per_mpi",
                                                             "1000",
                                                             asked_help,
                                                             "Number of particles to be generated and evolved by each MPI rank.");
    const size_t Npts = stoi(particles_string);  

    const std::string &output_frequency_string = input.getCmdOption("--output_frequency",
                                                                    "3600",
                                                                    asked_help,
                                                                    "Frequency [in seconds] to output particle positions.");
    const double out_freq = stod(output_frequency_string);  // in seconds

    const std::string &final_time_string = input.getCmdOption("--final_time", "-1", asked_help);
    const double final_time_input = stod(final_time_string);  // in seconds

    const std::string &particle_lifespan_string = input.getCmdOption("--particle_lifespan",
                                                                     "-1",
                                                                     asked_help,
                                                                     "Specifies how often particles are 'recycled' [i.e. randomly re-seeded]");
    const double particle_lifespan = stod(particle_lifespan_string);  // in seconds


    if (asked_help) { return 0; }

    // Set OpenMP thread number
    const int max_threads = omp_get_max_threads();
    omp_set_num_threads( max_threads );

    // Print some header info, depending on debug level
    print_header_info();

    std::vector<double> longitude, latitude, time, depth;
    std::vector<bool> mask;
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

    // Convert time units, if needed
    double time_scale_factor = 1.;
    if (time_units == "minutes") {
        time_scale_factor = 60.;
    } else if (time_units == "hours") {
        time_scale_factor = 60 *60.;
    } else if (time_units == "days") {
        time_scale_factor = 24 * 60 * 60.;
    }
    for ( II = 0; II < time.size(); ++II ) { time[II] = time[II] * time_scale_factor; }

    //  Verify final times make sense
    const int Ntime = time.size();
    if ( (final_time_input * time_scale_factor <= time.front()) and (Ntime == 1) ) {
        assert(false);  // Since only one time point is given, need a final particle time specified
                        // that is larger than the initial time
    } else if ( ( final_time_input > 0 ) and (Ntime > 1) and ( final_time_input * time_scale_factor > time.back() ) ) {
        assert(false); // Final time specified goes beyond provided data time
    }

    // Read in the velocity fields
    //  only load in two time-instances at a time to save on memory
    //  we'll simply roll through time and read in each new time as we need it
    // If Ntime == 1, then we're doing streamlines, so just load in the first time twice
    std::vector<double> u_lon_0, u_lon_1, u_lat_0, u_lat_1;
    read_var_from_file_at_time( u_lon_0, 0, zonal_vel_name, input_fname, &mask );
    read_var_from_file_at_time( u_lat_0, 0, merid_vel_name, input_fname, &mask );
    if (Ntime == 1) {
        fprintf(stdout, "Velocity data has a single time instance. Will compute streamlines.\n");
        u_lon_1 = u_lon_0;
        u_lat_1 = u_lat_0;
    } else {
        read_var_from_file_at_time( u_lon_1, 1, zonal_vel_name, input_fname, &mask );
        read_var_from_file_at_time( u_lat_1, 1, merid_vel_name, input_fname, &mask );
    }

    // Set the output times
    const double start_time = time.front(),
                 final_time = ( Ntime == 1 ) ? final_time_input * time_scale_factor : time.back();
    const size_t Nouts = std::max( (int) ( (final_time - start_time) / out_freq ), 2);

    #if DEBUG >= 1
    fprintf(stdout, " Output every %'g seconds, between %'g and %'g. Total of %'zu outputs.\n",
            out_freq, start_time, final_time, Nouts);
    #endif
    std::vector<double> target_times(Nouts);
    for ( II = 0; II < target_times.size(); ++II ) {
        target_times.at(II) = start_time + (double)II * ( final_time - start_time ) / Nouts;
    }

    // Get particle positions
    std::vector<double> starting_lat(Npts), starting_lon(Npts);
    particles_initial_positions(starting_lat, starting_lon, Npts, latitude, longitude, mask);

    // Trajectories dimension (essentially just a numbering)
    std::vector<double> trajectories(Npts);
    for (II = 0; II < Npts; ++II) { trajectories.at(II) = (double) II + wRank * Npts; };

    // List the fields to track along particle trajectories
    std::vector<const std::vector<double>*> fields_to_track;
    std::vector<std::string> names_of_tracked_fields;

    // Tracked fields is currently broken. Need to figure out how
    //  to make it work with the leap-frog loading system
    // Only vels are tracked by default
    names_of_tracked_fields.push_back( "vel_lon");
    names_of_tracked_fields.push_back( "vel_lat");

    // Storage for tracked fields
    std::vector< std::vector< double > > field_trajectories(names_of_tracked_fields.size());

    for (size_t Ifield = 0; Ifield < names_of_tracked_fields.size(); ++Ifield) {
        field_trajectories.at(Ifield).resize(Npts * Nouts, constants::fill_value);
    }

    std::vector<double> part_lon_hist(Npts * Nouts, constants::fill_value), 
                        part_lat_hist(Npts * Nouts, constants::fill_value);

    #if DEBUG >= 2
    fprintf(stdout, "Setting particle initial positions.\n");
    #endif
    for (II = 0; II < Npts; ++II) {
        part_lon_hist.at(II) = starting_lon.at(II);
        part_lat_hist.at(II) = starting_lat.at(II);
    }

    #if DEBUG >= 2
    fprintf(stdout, "Beginning evolution routine.\n");
    #endif
    // Now do the particle routine
    particles_evolve_trajectories(
        part_lon_hist,     part_lat_hist,
        field_trajectories,
        starting_lat,  starting_lon,
        target_times,
        particle_lifespan,
        u_lon_0, u_lon_1, u_lat_0, u_lat_1,
        zonal_vel_name, merid_vel_name, input_fname,
        fields_to_track, names_of_tracked_fields,
        time, latitude, longitude,        
        mask);

    fprintf(stdout, "\nProcessor %d of %d finished stepping particles.\n", wRank+1, wSize);

    std::vector<bool> out_mask(part_lon_hist.size());
    #pragma omp parallel \
    default (none) shared(out_mask, part_lon_hist) private(II) 
    {
        #pragma omp for collapse(1) schedule(static)
        for (II = 0; II < out_mask.size(); ++II) {
            out_mask.at(II) = (part_lon_hist.at(II) == constants::fill_value) ? false : true;
        }
    }


    // Initialize particle output file
    initialize_particle_file(target_times, trajectories, names_of_tracked_fields, output_fname);

    size_t starts[2], counts[2];
    starts[0] = 0;
    counts[0] = Nouts;

    starts[1] = wRank * Npts;
    counts[1] = Npts;

    MPI_Barrier(MPI_COMM_WORLD);
    write_field_to_output(part_lon_hist, "longitude", starts, counts, output_fname, &out_mask);
    write_field_to_output(part_lat_hist, "latitude",  starts, counts, output_fname, &out_mask);

    for (size_t Ifield = 0; Ifield < field_trajectories.size(); ++Ifield) {
        write_field_to_output(field_trajectories.at(Ifield), 
                names_of_tracked_fields.at(Ifield),  
                starts, counts, output_fname, &out_mask);
    }
    fprintf(stdout, "Processor %d of %d finished.\n", wRank+1, wSize);

    MPI_Finalize();
    return 0;
}
