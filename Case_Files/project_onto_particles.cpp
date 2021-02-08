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

    // first argument is the flag, second argument is default value (for when flag is not present)
    const std::string &time_dim_name      = input.getCmdOption("--time",        "time");
    const std::string &depth_dim_name     = input.getCmdOption("--depth",       "depth");
    const std::string &latitude_dim_name  = input.getCmdOption("--latitude",    "latitude");
    const std::string &longitude_dim_name = input.getCmdOption("--longitude",   "longitude");

    const std::string &part_latitude_name  = input.getCmdOption("--particle_latitude",  "latitude");
    const std::string &part_longitude_name = input.getCmdOption("--particle_longitude", "longitude");

    const std::string &Pi_var_name        = input.getCmdOption("--Pi",        "energy_transfer");
    const std::string &vort_var_name      = input.getCmdOption("--vort",      "vort_r");
    const std::string &fine_KE_var_name   = input.getCmdOption("--fine_KE",   "fine_KE");
    const std::string &coarse_KE_var_name = input.getCmdOption("--coarse_KE", "KE_from_coarse_vels");

    const std::string &Lambda_rot_var_name     = input.getCmdOption("--Lambda_rot",     "Lambda_rotational");
    const std::string &Lambda_nonlin_var_name  = input.getCmdOption("--Lambda_nonlin",  "Lambda_nonlinear");
    const std::string &Lambda_full_var_name    = input.getCmdOption("--Lambda_full",    "Lambda_full");

    const std::string &input_fields     = input.getCmdOption("--input_fields",       "input.nc");
    const std::string &input_trajectory = input.getCmdOption("--input_trajectories", "particles.nc");
    const std::string &output_name      = input.getCmdOption("--output",             "interpolated_particles.nc");

    // Print some header info, depending on debug level
    print_header_info();

    // Set OpenMP thread number
    const int max_threads = omp_get_max_threads();
    omp_set_num_threads( max_threads );

    std::vector<double> longitude, latitude, time, depth,
                        Pi, vort, fine_KE, coarse_KE,
                        Lambda_rot, Lambda_nonlin, Lambda_full,
                        particle_time, particle_traj, particle_lon, particle_lat;
    std::vector<bool>   mask, particle_mask;

    // Read in source data / get size information
    #if DEBUG >= 1
    if (wRank == 0) { fprintf(stdout, "Reading in source data.\n\n"); }
    #endif

    // Read in the grid coordinates
    read_var_from_file( longitude, longitude_dim_name, input_fields.c_str() );
    read_var_from_file( latitude,  latitude_dim_name,  input_fields.c_str() );
    read_var_from_file( time,      time_dim_name,      input_fields.c_str() );
    read_var_from_file( depth,     depth_dim_name,     input_fields.c_str() );
     
    convert_coordinates(longitude, latitude);

    // Read in the velocity fields
    read_var_from_file( Pi,         Pi_var_name.c_str(),        input_fields.c_str(), &mask, NULL, NULL, false );
    read_var_from_file( vort,       vort_var_name.c_str(),      input_fields.c_str(), NULL,  NULL, NULL, false );
    read_var_from_file( fine_KE,    fine_KE_var_name.c_str(),   input_fields.c_str(), NULL,  NULL, NULL, false );
    read_var_from_file( coarse_KE,  coarse_KE_var_name.c_str(), input_fields.c_str(), NULL,  NULL, NULL, false );

    read_var_from_file( Lambda_rot,     Lambda_rot_var_name.c_str(),     input_fields.c_str(), NULL,  NULL, NULL, false );
    read_var_from_file( Lambda_nonlin,  Lambda_nonlin_var_name.c_str(),  input_fields.c_str(), NULL,  NULL, NULL, false );
    read_var_from_file( Lambda_full,    Lambda_full_var_name.c_str(),    input_fields.c_str(), NULL,  NULL, NULL, false );

    // Read in the particle trajectories
    read_var_from_file( particle_time, "time",       input_trajectory.c_str() );
    read_var_from_file( particle_traj, "trajectory", input_trajectory.c_str() );

    std::vector<int> myCounts, myStarts;
    read_var_from_file( particle_lon, part_longitude_name, input_trajectory.c_str(), 
            &particle_mask, &myCounts, &myStarts, true, 1 );
    read_var_from_file( particle_lat, part_latitude_name,  input_trajectory.c_str(),
            &particle_mask, &myCounts, &myStarts, true, 1 );

    //const int Ntraj = particle_traj.size();
    //const int Ntime_traj = particle_time.size();
    const int Ntime_traj = myCounts[0];
    const int Ntraj      = myCounts[1];
    if (wRank == 0) { fprintf(stdout, "Ntraj = %'d    Ntime_traj = %'d\n", Ntraj, Ntime_traj); }

    // List the fields to track along particle trajectories
    std::vector<const std::vector<double>*> fields_to_track;
    std::vector<std::string> names_of_tracked_fields;

    names_of_tracked_fields.push_back( "vorticity" );
    fields_to_track.push_back( &vort );
    const int vort_ind = 0;

    names_of_tracked_fields.push_back( "Pi" );
    fields_to_track.push_back( &Pi );
    //const int Pi_ind = 1;

    names_of_tracked_fields.push_back( "Lambda_rot" );
    fields_to_track.push_back( &Lambda_rot );
    //const int Lambda_rot_ind = 2;

    names_of_tracked_fields.push_back( "Lambda_nonlin" );
    fields_to_track.push_back( &Lambda_nonlin );
    //const int Lambda_nonlin_ind = 3;

    names_of_tracked_fields.push_back( "Lambda_full" );
    fields_to_track.push_back( &Lambda_full );
    //const int Lambda_full_ind = 4;

    names_of_tracked_fields.push_back( "fine_KE" );
    fields_to_track.push_back( &fine_KE );
    const int fine_KE_ind = 5;

    names_of_tracked_fields.push_back( "coarse_KE" );
    fields_to_track.push_back( &coarse_KE );
    const int coarse_KE_ind = 6;

    // Storage for tracked fields
    std::vector< std::vector< double > > field_trajectories(fields_to_track.size());

    for (size_t Ifield = 0; Ifield < fields_to_track.size(); ++Ifield) {
        field_trajectories.at(Ifield).resize(Ntraj * Ntime_traj, constants::fill_value);
    }

    particles_project_onto_trajectory(
        field_trajectories, particle_time, particle_lat, particle_lon,
        fields_to_track, time, latitude, longitude, particle_mask, mask, myCounts);

    //
    //// Initialize particle output file
    ////    and write what we have so far
    //
    if (wRank == 0) { fprintf(stdout, "Initializing output file\n"); fflush(stdout); }
    names_of_tracked_fields.push_back("ddt_fine_KE");
    names_of_tracked_fields.push_back("ddt_coarse_KE");
    names_of_tracked_fields.push_back("ddt_enstrophy");
    initialize_projected_particle_file( particle_time, particle_traj, names_of_tracked_fields, 
            output_name.c_str());


    if (wRank == 0) { fprintf(stdout, "Done projection, outputting partial results.\n"); fflush(stdout); }
    size_t starts[2], counts[2];
    starts[0] = size_t(myStarts[0]);
    counts[0] = size_t(myCounts[0]);

    starts[1] = size_t(myStarts[1]);
    counts[1] = size_t(myCounts[1]);

    MPI_Barrier(MPI_COMM_WORLD);

    write_field_to_output(particle_lon, "longitude", starts, counts, output_name.c_str(), &particle_mask);
    write_field_to_output(particle_lat, "latitude",  starts, counts, output_name.c_str(), &particle_mask);

    for (size_t Ifield = 0; Ifield < fields_to_track.size(); ++Ifield) {
        MPI_Barrier(MPI_COMM_WORLD);
        write_field_to_output(  field_trajectories.at(Ifield), 
                                names_of_tracked_fields.at(Ifield).c_str(),  
                                starts, counts, output_name.c_str(), &particle_mask);
    }

    // Compute KE and enstrophy for each trajectory
    if (wRank == 0) { fprintf(stdout, "Computing enstrophy\n"); fflush(stdout); }
    int index;
    size_t II;
    std::vector<double> enstrophy( particle_lon.size(), constants::fill_value );
    #pragma omp parallel \
    default(none) shared( field_trajectories, enstrophy, particle_mask ) private( II )
    {
        #pragma omp for collapse(1) schedule(static)
        for (II = 0; II < enstrophy.size(); ++II) {
            if ( particle_mask.at(II) ) {
                enstrophy.at(II) = 0.5 * constants::rho0 * ( pow( field_trajectories.at( vort_ind ).at(II), 2.) );
            }
        }
    }

    // Compute derivatives ( ddt(KE), ddt(enstrophy) )
    if (wRank == 0) { fprintf(stdout, "Computing KE and enstrophy time derivatives\n"); fflush(stdout); }
    std::vector<double> ddt_coarse_KE( Ntraj * Ntime_traj, constants::fill_value ), 
                        ddt_fine_KE(   Ntraj * Ntime_traj, constants::fill_value ),
                        ddt_enstrophy( Ntraj * Ntime_traj, constants::fill_value ),
                        diff_vector;
    int LB_ret, diff_index, Itraj, Itime, ind;
    double tmp_c_KE, tmp_f_KE, tmp_EN;

    #pragma omp parallel \
    default(none) \
    shared( stdout, particle_mask, particle_time, ddt_coarse_KE, ddt_fine_KE, ddt_enstrophy, \
            enstrophy, field_trajectories ) \
    private( Itime, Itraj, index, tmp_c_KE, tmp_f_KE, tmp_EN, diff_index, ind, diff_vector, LB_ret )
    {
        #pragma omp for collapse(2) schedule(static)
        for (Itraj = 0; Itraj < Ntraj; ++Itraj) {
            for (Itime = 0; Itime < Ntime_traj; ++Itime) {

                index = Index(0, 0, Itime,      Itraj,
                              1, 1, Ntime_traj, Ntraj);

                if ( particle_mask.at(index) ) {
                    // build the differentiation vector 
                    //    we're going to call it 'lat', since it's the second-to-last dim
                    get_diff_vector( diff_vector, LB_ret, particle_time, "lat", 
                                     0, 0, Itime,      Itraj,
                                     1, 1, Ntime_traj, Ntraj,
                                     particle_mask, 1, constants::DiffOrd );

                    tmp_c_KE = 0.;
                    tmp_f_KE = 0.;
                    tmp_EN = 0.;
                    for ( ind = LB_ret; ind < LB_ret + (int) diff_vector.size(); ++ind ) {
                        diff_index = Index(0, 0, ind,        Itraj,
                                           1, 1, Ntime_traj, Ntraj);
                        tmp_c_KE += field_trajectories.at( coarse_KE_ind ).at(diff_index)
                                            * diff_vector.at(ind - LB_ret);
                        tmp_f_KE += field_trajectories.at( fine_KE_ind ).at(diff_index)
                                            * diff_vector.at(ind - LB_ret);
                        tmp_EN   += enstrophy.at(diff_index)
                                            * diff_vector.at(ind - LB_ret);
                    }

                    ddt_coarse_KE.at(index) = tmp_c_KE;
                    ddt_fine_KE.at(  index) = tmp_f_KE;
                    ddt_enstrophy.at(index) = tmp_EN;
                }
            }
        }
    }

    //
    write_field_to_output(ddt_fine_KE,   "ddt_fine_KE",   starts, counts, output_name.c_str(), &particle_mask);
    write_field_to_output(ddt_coarse_KE, "ddt_coarse_KE", starts, counts, output_name.c_str(), &particle_mask);
    write_field_to_output(ddt_enstrophy, "ddt_enstrophy", starts, counts, output_name.c_str(), &particle_mask);

    /*
    // Compute correlations
    if (wRank == 0) { fprintf(stdout, "Computing correlations\n"); fflush(stdout); }
    double ddt_cKE_mean, ddt_fKE_mean, ddt_EN_mean, Pi_mean, Lambda_mean,
           numer_cKE_Pi, numer_cKE_La, numer_fKE_Pi, numer_fKE_La, numer_EN_La, 
           denom_cKE, denom_fKE, denom_EN, denom_Pi, denom_La;
    std::vector<double> correl_cKE_Pi( Ntraj, constants::fill_value ), 
                        correl_cKE_La( Ntraj, constants::fill_value ), 
                        correl_fKE_Pi( Ntraj, constants::fill_value ), 
                        correl_fKE_La( Ntraj, constants::fill_value ), 
                        correl_EN_La(  Ntraj, constants::fill_value ),

                        out_mask(      Ntraj, 1 ),

                        numer_cKE_Pi_vec( Ntraj, constants::fill_value ), 
                        numer_cKE_La_vec( Ntraj, constants::fill_value ), 
                        numer_fKE_Pi_vec( Ntraj, constants::fill_value ), 
                        numer_fKE_La_vec( Ntraj, constants::fill_value ),
                        numer_EN_La_vec(  Ntraj, constants::fill_value ),

                        denom_cKE_vec( Ntraj, constants::fill_value ),
                        denom_fKE_vec( Ntraj, constants::fill_value ),
                        denom_Pi_vec(  Ntraj, constants::fill_value ),
                        denom_EN_vec(  Ntraj, constants::fill_value ),
                        denom_La_vec(  Ntraj, constants::fill_value ),

                        mean_ddt_cKE_vec( Ntraj, constants::fill_value ),
                        mean_ddt_fKE_vec( Ntraj, constants::fill_value ),
                        mean_ddt_EN_vec(  Ntraj, constants::fill_value ),
                        mean_Pi_vec(      Ntraj, constants::fill_value ),
                        mean_Lambda_vec(  Ntraj, constants::fill_value );
    int Itime_min, Itime_max, count;

    #pragma omp parallel \
    default(none) \
    shared( particle_mask, ddt_coarse_KE, ddt_fine_KE, ddt_enstrophy, field_trajectories, \
            correl_cKE_Pi, correl_cKE_La, correl_fKE_Pi, correl_fKE_La, correl_EN_La, out_mask,\
            numer_cKE_Pi_vec, numer_cKE_La_vec, numer_fKE_Pi_vec, numer_fKE_La_vec, numer_EN_La_vec,\
            denom_cKE_vec, denom_fKE_vec, denom_EN_vec, denom_Pi_vec, denom_La_vec,\
            mean_ddt_cKE_vec, mean_ddt_fKE_vec, mean_ddt_EN_vec, mean_Pi_vec, mean_Lambda_vec,\
            stdout ) \
    private( Itraj, Itime, index, ddt_cKE_mean, ddt_fKE_mean, ddt_EN_mean, Pi_mean, Lambda_mean, \
             numer_cKE_Pi, numer_cKE_La, numer_fKE_Pi, numer_fKE_La, numer_EN_La, \
             denom_cKE, denom_fKE, denom_EN, denom_Pi, denom_La,\
             Itime_max, Itime_min, count )
    {
        #pragma omp for collapse(1) schedule(static)
        for ( Itraj = 0; Itraj < Ntraj; ++Itraj ) {

            // get the starting index (since the first chunk might be masked)
            //   i.e. iterate until we hit a mask == 1 value
            Itime = 0;
            while ( Itime < Ntime_traj - constants::DiffOrd ) {
                index = Index(0, 0, Itime,      Itraj,
                              1, 1, Ntime_traj, Ntraj);
                if ( particle_mask.at(index) == 1) { break; }
                Itime++;
            }
            Itime += constants::DiffOrd / 2;
            Itime_min = Itime;

            // get the upper bound in time for the trajectory
            while ( Itime < Ntime_traj - constants::DiffOrd ) {
                index = Index(0, 0, Itime,      Itraj,
                              1, 1, Ntime_traj, Ntraj);

                if ( particle_mask.at(index) == 0) { break; }

                Itime++;
            }
            Itime_max = Itime;

            //
            //// First, compute the mean(s)
            //
            ddt_cKE_mean = particles_trajectory_mean( ddt_coarse_KE, 
                                                        particles_mask, Itraj, Ntraj, Itime_max, Itime_min );
            ddt_fKE_mean = particles_trajectory_mean( ddt_fine_KE,   
                                                        particles_mask, Itraj, Ntraj, Itime_max, Itime_min );
            ddt_EN_mean  = particles_trajectory_mean( ddt_enstrophy, 
                                                        particles_mask, Itraj, Ntraj, Itime_max, Itime_min );

            Pi_mean      = particles_trajectory_mean( field_trajectories.at(Pi_ind), 
                                                        particles_mask, Itraj, Ntraj, Itime_max, Itime_min );
            Lambda_mean  = particles_trajectory_mean( field_trajectories.at(Lambda_ind), 
                                                        particles_mask, Itraj, Ntraj, Itime_max, Itime_min );

            mean_ddt_cKE_vec.at( Itraj ) = ddt_cKE_mean;
            mean_ddt_fKE_vec.at( Itraj ) = ddt_fKE_mean;
            mean_ddt_EN_vec.at(  Itraj ) = ddt_EN_mean;
            mean_Pi_vec.at(      Itraj ) = Pi_mean;
            mean_Lambda_vec.at(  Itraj ) = Lambda_mean;


            // Now accumulate the correlation
            numer_cKE_Pi = 0.;
            numer_cKE_La = 0.;
            numer_fKE_Pi = 0.;
            numer_fKE_La = 0.;
            numer_EN_La  = 0.;

            denom_cKE   = 0.;
            denom_fKE   = 0.;
            denom_EN    = 0.;
            denom_Pi    = 0.;
            denom_La    = 0.;

            for (Itime = Itime_min; Itime < Itime_max; ++Itime) {

                // First, compute the mean
                index = Index(0, 0, Itime,      Itraj,
                              1, 1, Ntime_traj, Ntraj);

                numer_cKE_Pi +=    ( ddt_coarse_KE.at(index)                       - ddt_cKE_mean ) 
                                 * ( field_trajectories.at( Pi_ind ).at(index)     - Pi_mean );
                numer_cKE_La +=    ( ddt_coarse_KE.at(index)                       - ddt_cKE_mean ) 
                                 * ( field_trajectories.at( Lambda_ind ).at(index) - Lambda_mean );
                numer_fKE_Pi +=    ( ddt_fine_KE.at(index)                         - ddt_fKE_mean ) 
                                 * ( field_trajectories.at( Pi_ind ).at(index)     - Pi_mean );
                numer_fKE_La +=    ( ddt_fine_KE.at(index)                         - ddt_fKE_mean ) 
                                 * ( field_trajectories.at( Lambda_ind ).at(index) - Lambda_mean );
                numer_EN_La  +=    ( ddt_enstrophy.at(index)                       - ddt_EN_mean ) 
                                 * ( field_trajectories.at( Lambda_ind ).at(index) - Lambda_mean );

                denom_cKE += pow( ddt_coarse_KE.at(index)                        - ddt_cKE_mean , 2.);
                denom_fKE += pow( ddt_fine_KE.at(index)                          - ddt_fKE_mean , 2.);
                denom_EN  += pow( ddt_enstrophy.at(index)                        - ddt_EN_mean , 2.);
                denom_Pi  += pow( field_trajectories.at( Pi_ind ).at(index)      - Pi_mean ,     2.);
                denom_La  += pow( field_trajectories.at( Lambda_ind ).at(index)  - Lambda_mean , 2.);

            }

            numer_cKE_Pi_vec.at(Itraj) = numer_cKE_Pi;
            numer_cKE_La_vec.at(Itraj) = numer_cKE_La;
            numer_fKE_Pi_vec.at(Itraj) = numer_fKE_Pi;
            numer_fKE_La_vec.at(Itraj) = numer_fKE_La;
            numer_EN_La_vec.at( Itraj) = numer_EN_La;

            denom_cKE = sqrt( denom_cKE );
            denom_fKE = sqrt( denom_fKE );
            denom_EN  = sqrt( denom_EN  );
            denom_Pi  = sqrt( denom_Pi  );
            denom_La  = sqrt( denom_La  );

            denom_cKE_vec.at(Itraj) = denom_cKE;
            denom_fKE_vec.at(Itraj) = denom_fKE;
            denom_EN_vec.at(Itraj)  = denom_EN;
            denom_Pi_vec.at(Itraj)  = denom_Pi;
            denom_La_vec.at( Itraj) = denom_La;

            correl_cKE_Pi.at(Itraj) = numer_cKE_Pi / ( denom_cKE * denom_Pi );
            correl_cKE_La.at(Itraj) = numer_cKE_La / ( denom_cKE * denom_La );
            correl_fKE_Pi.at(Itraj) = numer_fKE_Pi / ( denom_fKE * denom_Pi );
            correl_fKE_La.at(Itraj) = numer_fKE_La / ( denom_fKE * denom_La );
            correl_EN_La.at( Itraj) = numer_EN_La  / ( denom_EN  * denom_La );
            
        }
    }


    // Write correlations
    size_t correl_starts[1], correl_counts[1];
    correl_starts[0] = (size_t) myStarts[1];
    correl_counts[0] = (size_t) myCounts[1];

    write_field_to_output(  correl_cKE_Pi, "correl_coarse_KE_Pi",
                            correl_starts, correl_counts, output_name.c_str(), &out_mask);
    write_field_to_output(  correl_cKE_La, "correl_coarse_KE_La",
                            correl_starts, correl_counts, output_name.c_str(), &out_mask);
    write_field_to_output(  correl_fKE_Pi, "correl_fine_KE_Pi",
                            correl_starts, correl_counts, output_name.c_str(), &out_mask);
    write_field_to_output(  correl_fKE_La, "correl_fine_KE_La",
                            correl_starts, correl_counts, output_name.c_str(), &out_mask);
    write_field_to_output(  correl_EN_La, "correl_EN_La",
                            correl_starts, correl_counts, output_name.c_str(), &out_mask);

    //
    #if DEBUG >= 2
    write_field_to_output( mean_ddt_cKE_vec, "mean_ddt_cKE", 
                            correl_starts, correl_counts, output_name.c_str(), &out_mask);
    write_field_to_output( mean_ddt_fKE_vec, "mean_ddt_fKE",
                            correl_starts, correl_counts, output_name.c_str(), &out_mask);
    write_field_to_output( mean_ddt_EN_vec,  "mean_ddt_EN",
                            correl_starts, correl_counts, output_name.c_str(), &out_mask);
    write_field_to_output( mean_Pi_vec,      "mean_Pi",   
                            correl_starts, correl_counts, output_name.c_str(), &out_mask);
    write_field_to_output( mean_Lambda_vec,  "mean_Lambda",
                            correl_starts, correl_counts, output_name.c_str(), &out_mask);

    write_field_to_output( numer_cKE_Pi_vec, "numer_cKE_Pi",
                            correl_starts, correl_counts, output_name.c_str(), &out_mask);
    write_field_to_output( numer_cKE_La_vec, "numer_cKE_La",
                            correl_starts, correl_counts, output_name.c_str(), &out_mask);
    write_field_to_output( numer_fKE_Pi_vec, "numer_fKE_Pi",
                            correl_starts, correl_counts, output_name.c_str(), &out_mask);
    write_field_to_output( numer_fKE_La_vec, "numer_fKE_La",
                            correl_starts, correl_counts, output_name.c_str(), &out_mask);
    write_field_to_output( numer_EN_La_vec,  "numer_EN_La",
                            correl_starts, correl_counts, output_name.c_str(), &out_mask);

    write_field_to_output( denom_cKE_vec, "denom_cKE",
                            correl_starts, correl_counts, output_name.c_str(), &out_mask);
    write_field_to_output( denom_fKE_vec, "denom_fKE",
                            correl_starts, correl_counts, output_name.c_str(), &out_mask);
    write_field_to_output( denom_EN_vec,  "denom_EN",
                            correl_starts, correl_counts, output_name.c_str(), &out_mask);
    write_field_to_output( denom_Pi_vec,  "denom_Pi",
                            correl_starts, correl_counts, output_name.c_str(), &out_mask);
    write_field_to_output( denom_La_vec,  "denom_La",
                            correl_starts, correl_counts, output_name.c_str(), &out_mask);
    #endif
    */

    MPI_Finalize();
    return 0;
}
