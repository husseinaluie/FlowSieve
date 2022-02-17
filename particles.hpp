#ifndef PARTICLES_HPP
#define PARTICLES_HPP 1

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <string>
#include <mpi.h>
#include "constants.hpp"


void particles_evolve_trajectories(
        std::vector<double> & part_lon_hist,
        std::vector<double> & part_lat_hist,
        std::vector<double> & rev_part_lon_hist,
        std::vector<double> & rev_part_lat_hist,
        std::vector< std::vector<double> > & field_trajectories,
        std::vector< std::vector<double> > & rev_field_trajectories,
        const std::vector<double> & starting_lat,
        const std::vector<double> & starting_lon,
        const std::vector<double> & target_times,
        const double particle_lifespan,
        const std::vector<double> & vel_lon,
        const std::vector<double> & vel_lat,
        const std::vector<const std::vector<double>*> & fields_to_track,
        const std::vector<std::string> & names_of_tracked_fields,
        const std::vector<double> & time,
        const std::vector<double> & lat,
        const std::vector<double> & lon,
        const std::vector<bool> & mask,
        const MPI_Comm comm = MPI_COMM_WORLD
        );

void particles_get_edges(
        int & left,
        int & right,
        int & bottom,
        int & top,
        const double & ref_lat,
        const double & ref_lon,
        const std::vector<double> & lat,
        const std::vector<double> & lon
        );

double particles_interp_from_edges(
        double ref_lat,
        double ref_lon,
        const std::vector<double> & lat,
        const std::vector<double> & lon,
        const std::vector<double> * field,
        const std::vector<bool> & mask,
        const int left,
        const int right,
        const int bottom,
        const int top,
        const double time_p,
        const int Itime,
        const int Ntime
        );

void particles_initial_positions(
        std::vector<double> & starting_lat, 
        std::vector<double> & starting_lon, 
        const int Npts,
        const std::vector<double> & latitude, 
        const std::vector<double> & longitude, 
        const std::vector<bool> & mask,
        const MPI_Comm comm = MPI_COMM_WORLD
        );

void particles_fore_back_difference(
        std::vector<double> & traj_dists,
        const std::vector<double> & fore_lon_hist,
        const std::vector<double> & fore_lat_hist,
        const std::vector<double> & back_lon_hist,
        const std::vector<double> & back_lat_hist,
        const MPI_Comm comm = MPI_COMM_WORLD
        );

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
        const MPI_Comm comm = MPI_COMM_WORLD
        );

#endif
