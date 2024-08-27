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
#include "../differentiation_tools.hpp"

int main(int argc, char *argv[]) {
    
    // PERIODIC_Y implies UNIFORM_LAT_GRID
    static_assert( (constants::UNIFORM_LAT_GRID) or (not(constants::PERIODIC_Y)),
            "PERIODIC_Y requires UNIFORM_LAT_GRID.\n"
            "Please update constants.hpp accordingly.\n");

    static_assert( ( (constants::PERIODIC_X) and (not(constants::PERIODIC_Y)) ),
            "The particles routine currently requires globe-like periodicity.\n"
            "Please update constants.hpp accordingly.\n");

    std::setlocale(LC_ALL, "");

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
        fprintf( stdout, "The command-line input arguments [and default values] are:\n" );
    }

    // first argument is the flag, second argument is default value (for when flag is not present)
    const std::string   &time_dim_name      = input.getCmdOption("--time",        "time",      asked_help),
                        &depth_dim_name     = input.getCmdOption("--depth",       "depth",     asked_help),
                        &latitude_dim_name  = input.getCmdOption("--latitude",    "latitude",  asked_help),
                        &longitude_dim_name = input.getCmdOption("--longitude",   "longitude", asked_help);

    const std::string &latlon_in_degrees  = input.getCmdOption("--is_degrees",   "true", asked_help);

    const std::string   &part_latitude_name  = input.getCmdOption("--particle_latitude",  "latitude", asked_help),
                        &part_longitude_name = input.getCmdOption("--particle_longitude", "longitude", asked_help);

    const std::string   &Nprocs_in_time_string  = input.getCmdOption("--Nprocs_in_time",  "1", asked_help),
                        &Nprocs_in_depth_string = input.getCmdOption("--Nprocs_in_depth", "1", asked_help);
    const int   Nprocs_in_time_input  = stoi(Nprocs_in_time_string),
                Nprocs_in_depth_input = stoi(Nprocs_in_depth_string);

    const std::string   &input_fields     = input.getCmdOption("--input_fields",       "input.nc", asked_help),
                        &input_trajectory = input.getCmdOption("--input_trajectories", "particles.nc", asked_help),
                        &output_name      = input.getCmdOption("--output",             "interpolated_particles.nc", asked_help);

    // Also read in the fields to be filtered from commandline
    //   e.g. --variables "rho salinity p" (names must match with input netcdf file)
    std::vector< std::string > vars_to_filter;
    input.getListofStrings( vars_to_filter, "--variables", asked_help );
    const size_t Nvars = vars_to_filter.size();

    if (asked_help) { return 0; }

    // Print some header info, depending on debug level
    print_header_info();

    // Set OpenMP thread number
    const int max_threads = omp_get_max_threads();
    omp_set_num_threads( max_threads );

    // Initialize dataset class instance
    dataset source_data;

    // Read in the grid coordinates
    source_data.load_time(      time_dim_name,      input_fields );
    source_data.load_depth(     depth_dim_name,     input_fields );
    source_data.load_latitude(  latitude_dim_name,  input_fields );
    source_data.load_longitude( longitude_dim_name, input_fields );

    // Apply some cleaning to the processor allotments if necessary. 
    source_data.check_processor_divisions( Nprocs_in_time_input, Nprocs_in_depth_input );
     
    // Convert to radians, if appropriate
    if ( ( latlon_in_degrees == "true" ) and not( constants::CARTESIAN ) ) { 
        convert_coordinates( source_data.longitude, source_data.latitude ); 
    }

    // Read in the scalar fields to filter
    /*
    #if DEBUG >= 1
    if (wRank == 0) { fprintf( stdout, "Reading in original fields.\n" ); }
    #endif
    for (size_t field_ind = 0; field_ind < vars_to_filter.size(); field_ind++) {
        source_data.load_variable( vars_to_filter.at(field_ind), vars_to_filter.at(field_ind), input_fields, true, true );
    }
    */


    // Read in the particle trajectories
    std::vector<double> particle_time, particle_traj, particle_lon, particle_lat;
    std::vector<bool> particle_mask, mask;

    read_var_from_file( particle_time, "time",       input_trajectory.c_str() );
    read_var_from_file( particle_traj, "trajectory", input_trajectory.c_str() );

    std::vector<int> myCounts_particles, myStarts_particles;
    read_var_from_file( particle_lon, part_longitude_name, input_trajectory.c_str(), 
            &particle_mask, &myCounts_particles, &myStarts_particles, 
            source_data.Nprocs_in_time, 1, true, 1 );
    read_var_from_file( particle_lat, part_latitude_name,  input_trajectory.c_str(),
            &particle_mask, &myCounts_particles, &myStarts_particles, 
            source_data.Nprocs_in_time, 1, true, 1 );

    //const int Ntraj = particle_traj.size();
    //const int Ntime_traj = particle_time.size();
    const size_t Ntime_traj = myCounts_particles[0];
    const size_t Ntraj      = myCounts_particles[1];
    if (wRank == 0) { fprintf(stdout, "Ntraj = %'zu    Ntime_traj = %'zu\n", Ntraj, Ntime_traj); }

    // Storage for tracked fields
    std::vector< std::vector< double > > field_trajectories( Nvars );
    for (size_t Ifield = 0; Ifield < Nvars; ++Ifield) {
        field_trajectories.at(Ifield).resize(Ntraj * Ntime_traj, constants::fill_value);
    }

    size_t Itime_traj, Itraj, field_ind, Itime = 0, prev_Itime = source_data.time.size() + 10;
    size_t Nlat = source_data.latitude.size(),
           Nlon = source_data.longitude.size(),
           Ntime = source_data.time.size();
    if (wRank == 0) { fprintf(stdout, "Beginning main interpolation loop.\n"); }
    std::vector<double> var_time_0( Nlat*Nlon, 0 ),
                        var_time_1( Nlat*Nlon, 0 );
    for ( Itime_traj = 0; Itime_traj < Ntime_traj; Itime_traj++ ) {
        Itime = std::upper_bound( source_data.time.begin(), source_data.time.end(), particle_time[Itime_traj] ) - source_data.time.begin();
        if (Itime > 0) { Itime--; }
        if ( Itime >= Ntime - 1 ) { Itime = Ntime - 2; }

        double t_p =   ( particle_time[Itime_traj] - source_data.time[Itime] ) 
                     / ( source_data.time[Itime+1] - source_data.time[Itime] );

        if (wRank == 0) {
            fprintf( stdout, " Particle -> ( %'g, %'zu ),  Field -> ( %'g - %'g, %'zu - %'zu / %'zu ), t_p = %'g\n",
                  particle_time[Itime_traj], Itime_traj,
                  source_data.time[Itime], source_data.time[Itime+1], 
                  Itime, Itime+1, source_data.time.size(),
                  t_p );
        }

        for ( field_ind = 0; field_ind < vars_to_filter.size(); field_ind++) {
            if ( (Itime_traj == 0) or (Itime != prev_Itime) ) {
                // If we're at a new Eulerian time, load the relevant data
                read_var_from_file_at_time( var_time_0, Itime,   vars_to_filter[field_ind], input_fields, &mask );
                read_var_from_file_at_time( var_time_1, Itime+1, vars_to_filter[field_ind], input_fields, &mask );
            }

            int left, right, bottom, top;
            double f0, f1, lon_p, lat_p;
            size_t particle_index;
            #pragma omp parallel \
            default(none) \
            shared( source_data, var_time_0, var_time_1, mask, field_trajectories, \
                    particle_lon, particle_lat, particle_mask ) \
            private( Itraj, particle_index, lon_p, lat_p, f0, f1,\
                     left, right, bottom, top ) \
            firstprivate( Itime_traj, Ntraj, field_ind, t_p )
            {
                #pragma omp for collapse(1) schedule(static)
                for ( Itraj = 0; Itraj < Ntraj; Itraj++ ) {

                    particle_index = Itime_traj * Ntraj + Itraj;

                    if ( not(particle_mask.at(particle_index)) ) {
                        field_trajectories.at(field_ind).at(particle_index) = constants::fill_value_double;
                        continue;
                    }

                    lon_p = particle_lon.at( particle_index );
                    lat_p = particle_lat.at( particle_index );

                    particles_get_edges(left, right, bottom, top, lat_p, lon_p, 
                            source_data.latitude, source_data.longitude);

                    f0 = particles_interp_from_edges( 
                                lat_p, lon_p,
                                source_data.latitude, source_data.longitude,
                                &var_time_0, mask,
                                left, right, bottom, top
                            );

                    f1 = particles_interp_from_edges( 
                                lat_p, lon_p,
                                source_data.latitude, source_data.longitude,
                                &var_time_1, mask,
                                left, right, bottom, top
                            );

                    field_trajectories.at(field_ind).at(particle_index) = f0 * ( 1 - t_p ) + f1 * t_p;
                }
            }
        }

        prev_Itime = Itime;
    }

    //
    //// Initialize particle output file
    ////    and write what we have so far
    //

    if (wRank == 0) { fprintf(stdout, "Done projection, outputting partial results.\n"); fflush(stdout); }
    size_t starts[2], counts[2];
    starts[0] = size_t(myStarts_particles[0]);
    counts[0] = size_t(myCounts_particles[0]);

    starts[1] = size_t(myStarts_particles[1]);
    counts[1] = size_t(myCounts_particles[1]);

    if (wRank == 0) { fprintf(stdout, "Initializing output file\n"); fflush(stdout); }
    initialize_projected_particle_file( particle_time, particle_traj, vars_to_filter, 
            output_name.c_str(), starts, counts);

    MPI_Barrier(MPI_COMM_WORLD);

    write_field_to_output(particle_lon, "longitude", starts, counts, output_name.c_str(), &particle_mask);
    write_field_to_output(particle_lat, "latitude",  starts, counts, output_name.c_str(), &particle_mask);

    for (size_t Ifield = 0; Ifield < Nvars; ++Ifield) {
        MPI_Barrier(MPI_COMM_WORLD);
        write_field_to_output(  field_trajectories.at(Ifield), 
                                vars_to_filter.at(Ifield),  
                                starts, counts, output_name, &particle_mask
                                );
    }

    MPI_Finalize();
    return 0;
}
