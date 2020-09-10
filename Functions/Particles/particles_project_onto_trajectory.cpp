#include <math.h>
#include <algorithm>
#include <vector>
#include <mpi.h>
#include <omp.h>
#include "../../constants.hpp"
#include "../../functions.hpp"
#include "../../particles.hpp"

void particles_project_onto_trajectory(
        std::vector< std::vector<double> > & field_trajectories,
        const std::vector<double> & trajectory_time,
        const std::vector<double> & trajectory_lat,
        const std::vector<double> & trajectory_lon,
        const std::vector<const std::vector<double>*> & fields_to_track,
        const std::vector<double> & time,
        const std::vector<double> & lat,
        const std::vector<double> & lon,
        const std::vector<bool> & particle_mask,
        const std::vector<bool> & field_mask,
        const std::vector<int> & myCounts,
        const MPI_Comm comm
        ) {


    int wRank, wSize;
    MPI_Comm_rank( comm, &wRank );
    MPI_Comm_size( comm, &wSize );

    double traj_time, traj_lat, traj_lon,
           field_val, time_p;

    const int Ntime      = time.size(),
              Ntime_traj = myCounts[0],
              Ntraj      = myCounts[1];

    if (wRank == 0) {
        fprintf(stdout, "Beginning particle projections (%'d trajectories, %'d points per trajectory).\n",
                Ntraj, Ntime_traj);
        fflush(stdout);
    }

    int left, right, bottom, top,
        ref_ind = 0, Ip,
        print_count = 1;
    size_t index;
    bool do_print = false;

    for ( int Itime_traj = 0; Itime_traj < Ntime_traj; ++Itime_traj ) {

        do_print = Itime_traj + 1 > 0.05 * print_count * Ntime_traj ;

        if (do_print and (wRank == 0)) {
            fprintf(stdout, "  ... done %d percent\n", print_count * 5);
            fflush(stdout);
            print_count++;
        }

        traj_time = trajectory_time.at(Itime_traj);

        while (traj_time > time.at(ref_ind+1)) { ref_ind++; }

        time_p =    ( traj_time          - time.at(ref_ind) ) 
                  / ( time.at(ref_ind+1) - time.at(ref_ind) );

        size_t Ifield;
        #pragma omp parallel \
        default(none) \
        shared( time, lat, lon, particle_mask, field_mask, Itime_traj, traj_time, time_p,\
                ref_ind, trajectory_time, trajectory_lat, trajectory_lon,\
                field_trajectories, fields_to_track, stdout, do_print, wRank )\
        private( Ip, index, Ifield, traj_lat, traj_lon, field_val,\
                 left, right, bottom, top)
        {
            #pragma omp for collapse(1) schedule(static)
            for (Ip = 0; Ip < Ntraj; ++Ip) {

                index = Index(0, 0, Itime_traj, Ip,
                              1, 1, Ntime_traj, Ntraj);

                for (Ifield = 0; Ifield < fields_to_track.size(); ++Ifield) {
                    field_trajectories.at(Ifield).at(index) = constants::fill_value;
                }

                if ( particle_mask.at(index) ) {

                    traj_lon = trajectory_lon.at(index);
                    traj_lat = trajectory_lat.at(index);

                    particles_get_edges(left, right, bottom, top, traj_lat, traj_lon, lat, lon);
                    if ( (bottom >= 0) and (top >= 0) ) { 

                        for (Ifield = 0; Ifield < fields_to_track.size(); ++Ifield) {
                            field_val = particles_interp_from_edges(traj_lat, traj_lon, lat, lon, 
                                                                    fields_to_track.at(Ifield), field_mask,
                                                                    left, right, bottom, top, time_p, ref_ind, Ntime);
                            field_trajectories.at(Ifield).at(index) = field_val;

                            if (do_print and (wRank == 0) and (Ip == 0)) {
                                fprintf(stdout, "    field %'zu has value %'g\n", Ifield, field_val);
                                fflush(stdout);
                            }

                        }
                    }
                }
            }
        } 
    }
}
