#ifndef PARTICLES_HPP
#define PARTICLES_HPP 1

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <string>
#include "constants.hpp"


void particles_evolve_trajectories(
        std::vector<double> & part_lon_hist,
        std::vector<double> & part_lat_hist,
        const std::vector<double> & starting_lat,
        const std::vector<double> & starting_lon,
        const std::vector<double> & target_times,
        const std::vector<double> & vel_lon,
        const std::vector<double> & vel_lat,
        const std::vector<double> & time,
        const std::vector<double> & lat,
        const std::vector<double> & lon,
        const std::vector<double> & mask
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
        const std::vector<double> lat,
        const std::vector<double> lon,
        const std::vector<double> field,
        const std::vector<double> mask,
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
        const std::vector<double> & mask
        );

#endif
