#include <cassert>
#include <ctime>
#include <math.h>
#include <algorithm>
#include <vector>
#include <mpi.h>
#include <omp.h>
#include <random>
#include "../../constants.hpp"
#include "../../functions.hpp"
#include "../../particles.hpp"
#include "../../netcdf_io.hpp"

bool check_if_recycling( 
        const double dt, 
        const double time, 
        const double next_recycle_time,
        const double particle_lifespan
        ) {
    // Check if this particle is going to recycle
    //      i.e. if it gets reset to a new random location
    // A non-positive lifespan means no recycling
    if ( particle_lifespan <= 0 ) { return false; }
    if ( time >= next_recycle_time ) { return true; }
    return false; // back-up catch-all
}

void recycle_position( 
        double & lon0,
        double & lat0,
        const std::vector<double> & longitude,
        const std::vector<double> & latitude,
        const double rand1, // [0,1]
        const double rand2  // [0,1]
        ) {


    // reset longitude
    lon0 = rand1 * longitude.front() + (1.-rand1) * longitude.back();
    
    // reset latitude [with cosine weighting]
    double sin_LB = sin( latitude.front() );
    double sin_UB = sin( latitude.back() );
    lat0 = asin( rand2 * sin_UB + (1.-rand2) * sin_LB );
}

void move_on_sphere(
        double & x_dest,
        double & y_dest,
        const double x_init,
        const double y_init,
        const double u,
        const double v
        ) {

    x_dest = x_init + u / ( constants::R_earth * cos(y_init) );
    y_dest = y_init + v / constants::R_earth;

    // Account for moving across the poles
    if (y_dest >  M_PI/2) { y_dest -= y_dest - ( M_PI/2); x_dest += M_PI; }
    if (y_dest < -M_PI/2) { y_dest -= y_dest - (-M_PI/2); x_dest += M_PI; }

    // Adjust lon for periodicity
    if (x_dest >  M_PI) { x_dest -= 2 * M_PI; }
    if (x_dest < -M_PI) { x_dest += 2 * M_PI; }

}

double get_at_point( 
        const double t, 
        const double lon0, 
        const double lat0,
        const std::vector<double> & var0,
        const std::vector<double> & var1,
        const double t0,
        const double t1,
        const std::vector<double> & lat,
        const std::vector<double> & lon,
        const std::vector<bool> & mask
        ) {

    int left, right, bottom, top;
    particles_get_edges(left, right, bottom, top, lat0, lon0, lat, lon);

    double f0, f1;
    f0 = particles_interp_from_edges( lat0, lon0, lat, lon, &var0, mask, left, right, bottom, top );
    f1 = particles_interp_from_edges( lat0, lon0, lat, lon, &var1, mask, left, right, bottom, top );

    double t_p = ( t - t0 ) / (t1 - t0);

    double var_val = f0 * (1 - t_p) + f1 * t_p;

    return var_val;
}

void SimplecticEuler_update(
        double & x_new,
        double & y_new,
        double & dt_new,
        const double t0,
        const double x0,
        const double y0,
        const double dt_seed,
        const double dt_max,
        const double dt_min,
        const std::vector<double> & u_lon_0,
        const std::vector<double> & u_lon_1,
        const std::vector<double> & u_lat_0,
        const std::vector<double> & u_lat_1,
        const double tb0,
        const double tb1,
        const std::vector<double> & lat,
        const std::vector<double> & lon,
        const std::vector<bool> & mask
        ) {

    double dt = dt_seed;

    double u0 = get_at_point( t0, x0, y0, u_lon_0, u_lon_1, tb0, tb1, lat, lon, mask );
    double x1, y1;
    move_on_sphere( x1, y1, x0, y0, dt * u0, 0. );

    double v0 = get_at_point( t0, x1, y0, u_lat_0, u_lat_1, tb0, tb1, lat, lon, mask );
    double x2, y2;
    move_on_sphere( x2, y2, x1, y1, 0., dt * v0 );

    dt_new = dt;
    x_new = x2;
    y_new = y2;
}

void RKF45_update( // Fehlbers's RK 4(5)
        double & x_new,
        double & y_new,
        double & dt_new,
        const double t0,
        const double x0,
        const double y0,
        const double dt_seed,
        const double dt_max,
        const double dt_min,
        const std::vector<double> & u_lon_0,
        const std::vector<double> & u_lon_1,
        const std::vector<double> & u_lat_0,
        const std::vector<double> & u_lat_1,
        const double tb0,
        const double tb1,
        const std::vector<double> & lat,
        const std::vector<double> & lon,
        const std::vector<bool> & mask
        ) {

    // Table 5.1 (page 177 of "Solving Ordinary Differential Equations, 2nd Ed.", Hairer, Norsett, Wanner

    // Butchers tableau
    const double 
        A10 = 1./4,
        A20 = 3./32,        A21 = 9./32,
        A30 = 1932./2197,   A31 = -7200./2197,  A32 = 7296./2197,
        A40 = 439./216,     A41 = -8.,          A42 = 3680./513,    A43 = -845./4104,
        A50 = -8./27,       A51 = 2.,           A52 = -3544./2565,  A53 = 1859./4104,   A54 = -11./40;

    const double    C0 = 0, 
                    C1 = 1./4, 
                    C2 = 3./8,
                    C3 = 12./13,
                    C4 = 1.,
                    C5 = 1./2;

    const double B0 = 25./216, B1 = 0., B2 = 1408./2565, B3 = 2197./4104, B4 = -1./5, B5 = 0.;
    const double B0hat = 16./135, B1hat = 0., B2hat = 6656./12825, B3hat = 28561./56430, B4hat = -9./50, B5hat = 2./55;

    // and a pile of working variables
    double t, u0, u1, u2, u3, u4, u5, 
              v0, v1, v2, v3, v4, v5, 
                  x1, x2, x3, x4, x5, x_final, 
                  y1, y2, y3, y4, y5, y_final;
    double dt = dt_seed, dt_tmp = dt_seed;

    const double vel_tol = 1e-6, // in m/s, will then be scaled by dt
                 eps     = 1.;   // if displacement error is > vel_tol * dt, reduce step size
    double TE = 10 * eps;

    while ( TE > eps ) {

        dt = dt_tmp;

        // Note: the u,v values have dt multiplied in, so they are all displacements in m
    
        // First increment
        t = t0 + C0 * dt;
        u0 = dt * get_at_point( t, x0, y0, u_lon_0, u_lon_1, tb0, tb1, lat, lon, mask );
        v0 = dt * get_at_point( t, x0, y0, u_lat_0, u_lat_1, tb0, tb1, lat, lon, mask );

        // Second increment
        t = t0 + C1 * dt;
        move_on_sphere( x1, y1, x0, y0, A10 * u0, A10 * v0 );
        u1 = dt * get_at_point( t, x1, y1, u_lon_0, u_lon_1, tb0, tb1, lat, lon, mask );
        v1 = dt * get_at_point( t, x1, y1, u_lat_0, u_lat_1, tb0, tb1, lat, lon, mask );

        // Third increment
        t = t0 + C2 * dt;
        move_on_sphere( x2, y2, x0, y0, A20 * u0 + A21 * u1, A20 * v0 + A21 * v1 );
        u2 = dt * get_at_point( t, x2, y2, u_lon_0, u_lon_1, tb0, tb1, lat, lon, mask );
        v2 = dt * get_at_point( t, x2, y2, u_lat_0, u_lat_1, tb0, tb1, lat, lon, mask );

        // Fourth increment
        t = t0 + C3 * dt;
        move_on_sphere( x3, y3, x0, y0, A30 * u0 + A31 * u1 + A32 * u2, 
                                        A30 * v0 + A31 * v1 + A32 * v2 );
        u3 = dt * get_at_point( t, x3, y3, u_lon_0, u_lon_1, tb0, tb1, lat, lon, mask );
        v3 = dt * get_at_point( t, x3, y3, u_lat_0, u_lat_1, tb0, tb1, lat, lon, mask );

        // Fifth increment
        t = t0 + C4 * dt;
        move_on_sphere( x4, y4, x0, y0, A40 * u0 + A41 * u1 + A42 * u2 + A43 * u3, 
                                        A40 * v0 + A41 * v1 + A42 * v2 + A43 * v3);
        u4 = dt * get_at_point( t, x4, y4, u_lon_0, u_lon_1, tb0, tb1, lat, lon, mask );
        v4 = dt * get_at_point( t, x4, y4, u_lat_0, u_lat_1, tb0, tb1, lat, lon, mask );

        // Sixth increment
        t = t0 + C5 * dt;
        move_on_sphere( x5, y5, x0, y0, A50 * u0 + A51 * u1 + A52 * u2 + A53 * u3 + A54 * u4, 
                                        A50 * v0 + A51 * v1 + A52 * v2 + A53 * v3 + A54 * v4);
        u5 = dt * get_at_point( t, x5, y5, u_lon_0, u_lon_1, tb0, tb1, lat, lon, mask );
        v5 = dt * get_at_point( t, x5, y5, u_lat_0, u_lat_1, tb0, tb1, lat, lon, mask );


        // Net displacement vector
        double  u_final_LO = B0 * u0 + B1 * u1 + B2 * u2 + B3 * u3 + B4 * u4 * B5 * u5,
                v_final_LO = B0 * v0 + B1 * v1 + B2 * v2 + B3 * v3 + B4 * v4 * B5 * v5;
        double  u_final_HO = B0hat * u0 + B1hat * u1 + B2hat * u2 + B3hat * u3 + B4hat * u4 * B5hat * u5,
                v_final_HO = B0hat * v0 + B1hat * v1 + B2hat * v2 + B3hat * v3 + B4hat * v4 * B5hat * v5;

        // The error estimate is the different in the two displacement vectors
        double displacement_error = sqrt( pow( u_final_LO - u_final_HO, 2.) + pow( v_final_LO - v_final_HO, 2. ) );
        TE = displacement_error / ( vel_tol * dt );

        // For numerical reasons, we advance forward with the low-order solution to preserve our error estimates
        move_on_sphere( x_final, y_final, x0, y0, u_final_LO, v_final_LO );

        // If the error was too large, then we'll need to redo the step with a smaller dt.
        // If the error was smaller than our tolerance, than the step is adaptively increased
        if ( TE == 0 ) {
            dt_tmp = dt_max;
        } else {
            dt_tmp = std::fmax( std::fmin( 0.9 * dt * pow( eps / TE, 1./5 ), dt_max ), dt_min );
        }

        if (dt_tmp <= dt_min) { break; } // if we're at the smallest dt, just stop
    }

    // After success, store the desired values
    dt_new = dt;
    x_new = x_final;
    y_new = y_final;
}

void DOPRI5_update( // Dormand-Prince 5(4)
        double & x_new,
        double & y_new,
        double & dt_new,
        const double t0,
        const double x0,
        const double y0,
        const double dt_seed,
        const double dt_max,
        const double dt_min,
        const std::vector<double> & u_lon_0,
        const std::vector<double> & u_lon_1,
        const std::vector<double> & u_lat_0,
        const std::vector<double> & u_lat_1,
        const double tb0,
        const double tb1,
        const std::vector<double> & lat,
        const std::vector<double> & lon,
        const std::vector<bool> & mask
        ) {

    // Table 5.2 (page 178 of "Solving Ordinary Differential Equations, 2nd Ed.", Hairer, Norsett, Wanner)


    // Butcher tableau
    const double 
        A10 = 1./5,
        A20 = 3./40,        A21 = 9./40,
        A30 = 44./45,       A31 = -56./15,      A32 = 32./9,
        A40 = 19372./6561,  A41 = -25360./2187, A42 = 64448./6561,  A43 = -212./729,
        A50 = 9017./3168,   A51 = -355./33,     A52 = 46732./5247,  A53 = 49./176,   A54 = -5103./18656,
        A60 = 35./384,      A61 = 0.,           A62 = 500./1113,    A63 = 125./192,  A64 = -2187./6784,   A65 = 11./84;

    const double    C0 = 0, 
                    C1 = 1./5, 
                    C2 = 3./10,
                    C3 = 4./5,
                    C4 = 8./9,
                    C5 = 1.,
                    C6 = 1.;

    const double B0    = 35./384,     B1    = 0., B2    = 500./1113,   B3    = 125./192, B4    = -2187./6784,    
                 B5    = 11./84,      B6    = 0.;
    const double B0hat = 5179./57600, B1hat = 0., B2hat = 7571./16695, B3hat = 393./640, B4hat = -92097./339200, 
                 B5hat = 187./2100,   B6hat = 1./40;

    // and a pile of working variables
    double t, u0, u1, u2, u3, u4, u5, u6,
              v0, v1, v2, v3, v4, v5, v6,
                  x1, x2, x3, x4, x5, x6, x_final, 
                  y1, y2, y3, y4, y5, y6, y_final;
    double dt = dt_seed, dt_tmp = dt_seed;

    const double vel_tol = 1e-6, // in m/s, will then be scaled by dt to give displacement tolerance
                 eps     = 1.;   // if displacement error is > vel_tol * dt, reduce step size
    double TE = 10 * eps;

    while ( TE > eps ) {

        dt = dt_tmp;
    
        // First increment
        t = t0 + C0 * dt;
        u0 = dt * get_at_point( t, x0, y0, u_lon_0, u_lon_1, tb0, tb1, lat, lon, mask );
        v0 = dt * get_at_point( t, x0, y0, u_lat_0, u_lat_1, tb0, tb1, lat, lon, mask );

        // Second increment
        t = t0 + C1 * dt;
        move_on_sphere( x1, y1, x0, y0, A10 * u0, A10 * v0 );
        u1 = dt * get_at_point( t, x1, y1, u_lon_0, u_lon_1, tb0, tb1, lat, lon, mask );
        v1 = dt * get_at_point( t, x1, y1, u_lat_0, u_lat_1, tb0, tb1, lat, lon, mask );

        // Third increment
        t = t0 + C2 * dt;
        move_on_sphere( x2, y2, x0, y0, A20 * u0 + A21 * u1, A20 * v0 + A21 * v1 );
        u2 = dt * get_at_point( t, x2, y2, u_lon_0, u_lon_1, tb0, tb1, lat, lon, mask );
        v2 = dt * get_at_point( t, x2, y2, u_lat_0, u_lat_1, tb0, tb1, lat, lon, mask );

        // Fourth increment
        t = t0 + C3 * dt;
        move_on_sphere( x3, y3, x0, y0, A30 * u0 + A31 * u1 + A32 * u2, 
                                        A30 * v0 + A31 * v1 + A32 * v2 );
        u3 = dt * get_at_point( t, x3, y3, u_lon_0, u_lon_1, tb0, tb1, lat, lon, mask );
        v3 = dt * get_at_point( t, x3, y3, u_lat_0, u_lat_1, tb0, tb1, lat, lon, mask );

        // Fifth increment
        t = t0 + C4 * dt;
        move_on_sphere( x4, y4, x0, y0, A40 * u0 + A41 * u1 + A42 * u2 + A43 * u3, 
                                        A40 * v0 + A41 * v1 + A42 * v2 + A43 * v3);
        u4 = dt * get_at_point( t, x4, y4, u_lon_0, u_lon_1, tb0, tb1, lat, lon, mask );
        v4 = dt * get_at_point( t, x4, y4, u_lat_0, u_lat_1, tb0, tb1, lat, lon, mask );

        // Sixth increment
        t = t0 + C5 * dt;
        move_on_sphere( x5, y5, x0, y0, A50 * u0 + A51 * u1 + A52 * u2 + A53 * u3 + A54 * u4, 
                                        A50 * v0 + A51 * v1 + A52 * v2 + A53 * v3 + A54 * v4);
        u5 = dt * get_at_point( t, x5, y5, u_lon_0, u_lon_1, tb0, tb1, lat, lon, mask );
        v5 = dt * get_at_point( t, x5, y5, u_lat_0, u_lat_1, tb0, tb1, lat, lon, mask );

        // Seventh increment
        t = t0 + C6 * dt;
        move_on_sphere( x6, y6, x0, y0, A60 * u0 + A61 * u1 + A62 * u2 + A63 * u3 + A64 * u4 + A65 * u5, 
                                        A60 * v0 + A61 * v1 + A62 * v2 + A63 * v3 + A64 * v4 + A65 * v5);
        u6 = dt * get_at_point( t, x6, y6, u_lon_0, u_lon_1, tb0, tb1, lat, lon, mask );
        v6 = dt * get_at_point( t, x6, y6, u_lat_0, u_lat_1, tb0, tb1, lat, lon, mask );


        // Net displacement vector
        double  u_final_LO = B0 * u0 + B1 * u1 + B2 * u2 + B3 * u3 + B4 * u4 * B5 * u5 + B6 * u6,
                v_final_LO = B0 * v0 + B1 * v1 + B2 * v2 + B3 * v3 + B4 * v4 * B5 * v5 + B6 * v6;
        double  u_final_HO = B0hat * u0 + B1hat * u1 + B2hat * u2 + B3hat * u3 + B4hat * u4 * B5hat * u5 + B6hat * u6,
                v_final_HO = B0hat * v0 + B1hat * v1 + B2hat * v2 + B3hat * v3 + B4hat * v4 * B5hat * v5 + B6hat * v6;

        // The error estimate is the different in the two displacement vectors
        double displacement_error = sqrt( pow( u_final_LO - u_final_HO, 2.) + pow( v_final_LO - v_final_HO, 2. ) );
        TE = displacement_error / ( vel_tol * dt );

        // For numerical reasons, we advance forward with the low-order solution to preserve our error estimates
        move_on_sphere( x_final, y_final, x0, y0, u_final_LO, v_final_LO );

        // If the error was too large, then we'll need to redo the step with a smaller dt.
        // If the error was smaller than our tolerance, than the step is adaptively increased
        if ( TE == 0 ) {
            dt_tmp = dt_max;
        } else {
            dt_tmp = std::fmax( std::fmin( 0.9 * dt * pow( eps / TE, 1./5 ), dt_max ), dt_min );
        }

        if (dt_tmp <= dt_min) { break; } // if we're at the smallest dt, just stop
    }

    // After success, store the desired values
    dt_new = dt;
    x_new = x_final;
    y_new = y_final;
}

void particles_evolve_trajectories(
        std::vector<double> & part_lon_hist,
        std::vector<double> & part_lat_hist,
        std::vector< std::vector<double> > & field_trajectories,
        const std::vector<double> & starting_lat,
        const std::vector<double> & starting_lon,
        const std::vector<double> & target_times,
        const double particle_lifespan,
        std::vector<double> & vel_lon_0,
        std::vector<double> & vel_lon_1,
        std::vector<double> & vel_lat_0,
        std::vector<double> & vel_lat_1,
        const std::string zonal_vel_name,
        const std::string merid_vel_name,
        const std::string input_fname,
        const std::vector<const std::vector<double>*> & fields_to_track,
        const std::vector<std::string> & names_of_tracked_fields,
        const std::vector<double> & time,
        const std::vector<double> & lat,
        const std::vector<double> & lon,
        const std::vector<bool> & mask,
        const MPI_Comm comm
        ) {

    int wRank, wSize;
    MPI_Comm_rank( comm, &wRank );
    MPI_Comm_size( comm, &wSize );

    double t_part, lon0, lat0, lon_new, lat_new,
           dx_loc, dy_loc, dt, dt_max, dt_min, dt_new,
           up_0, up_1, vp_0, vp_1,
           vel_lon_part, vel_lat_part, field_val,
           test_val, time_p, time_block_0, time_block_1,
           data_load_time;

    const double dlon = lon.at(1) - lon.at(0),
                 dlat = lat.at(1) - lat.at(0),
                 cfl_max  = 1e-1, // this is for the upper-bound of dt
                 cfl_min  = 1e-5, // this is for the lower-bound of dt
                 U0   = 2.,
                 dt_target = target_times.at(1) - target_times.at(0);

    const unsigned int  Nlat   = lat.size(),
                        Nlon   = lon.size(),
                        Ntime  = time.size(),
                        Ndepth = 1,
                        Nparts = starting_lat.size(),
                        Nouts  = target_times.size();

    int left, right, bottom, top;
    bool do_recycle;

    std::vector<double> lon_tmps = starting_lon, 
                        lat_tmps = starting_lat,
                        recycle_times(Nparts,0);

    unsigned int out_ind=0, prev_out_ind, step_iter, Ip;
    size_t index, next_load_index = 1;

    int perc_base = 5;
    int perc = 0, perc_count=0;

    std::random_device rd;  // seed generator for random
    std::mt19937_64 gen(rd());  // replace rd() with a number to fix the seed
    std::uniform_real_distribution<> get_rand(0, 1);

    if ( particle_lifespan > 0 ) {
        for (Ip = 0; Ip < Nparts; Ip++) {
            if ( constants::PARTICLE_RECYCLE_TYPE == constants::ParticleRecycleType::FixedInterval ) {
                recycle_times[Ip] = time.at(0) + particle_lifespan;
            } else if (constants::PARTICLE_RECYCLE_TYPE == constants::ParticleRecycleType::Stochastic) {
                recycle_times[Ip] = time.at(0) - log( get_rand(gen) ) * particle_lifespan; 
            }
        }
    }

    // Step through time 'blocks' i.e. get to the time of the next velocity data
    //  each particle is using adaptive stepping, so just loop through them getting there
    // We've already loaded in the first two times, so just flag the target time and continue.
    prev_out_ind = 0;
    for ( next_load_index = 1; next_load_index < (unsigned int) std::fmax(2.,(double)Ntime); next_load_index++ ) {

        #if DEBUG >= 0
        if ( (wRank == 0) and (Ntime > 1) ) {
            // Every perc_base percent, print a dot, but only the first thread
            while ( ((double)(next_load_index+1) / Ntime) * 100 >= perc ) {
                perc_count++;
                if (perc_count % 5 == 0) { fprintf(stdout, "|"); }
                else                     { fprintf(stdout, "."); }
                fflush(stdout);
                perc += perc_base;
            }
        }
        #endif

        data_load_time = (Ntime == 1) ? target_times.back() : time.at(next_load_index);

        if (next_load_index > 1) {
            // Swap time 1 to time 0
            vel_lon_1.swap(vel_lon_0);
            vel_lat_1.swap(vel_lat_0);

            // Load in next time-block as time 1
            read_var_from_file_at_time( vel_lon_1, next_load_index, zonal_vel_name, input_fname );
            read_var_from_file_at_time( vel_lat_1, next_load_index, merid_vel_name, input_fname );
        }
        #pragma omp parallel \
        default(none) \
        shared( lat, lon, vel_lon_0, vel_lon_1, vel_lat_0, vel_lat_1, mask,\
                target_times, time, part_lon_hist, part_lat_hist, recycle_times,\
                field_trajectories, fields_to_track, lon_tmps, lat_tmps,\
                wRank, wSize, stdout)\
        private(Ip, index, t_part, step_iter, lon0, lat0, \
                dt_max, dt_min, dt_new, lon_new, lat_new, \
                up_0, up_1, vp_0, vp_1, \
                dx_loc, dy_loc, dt, time_p, vel_lon_part, vel_lat_part, field_val,\
                do_recycle, time_block_0, time_block_1, \
                left, right, bottom, top) \
        firstprivate( Nparts, Nouts, Ntime, dlon, dlat, dt_target, data_load_time, \
                      particle_lifespan, next_load_index, prev_out_ind ) \
        reduction( max:out_ind )
        {

            // Make a thread-local random generator
            std::random_device rd;  // seed generator for random
            std::mt19937_64 gen(rd());  // replace rd() with a number to fix the seed
            std::uniform_real_distribution<> get_rand(0, 1);

            #pragma omp for collapse(1) schedule(guided)
            for (Ip = 0; Ip < Nparts; ++Ip) {

                do_recycle = false;
                dt = 0.;

                //
                time_block_0 = time.at(next_load_index-1);
                time_block_1 = data_load_time;

                // Particle time
                t_part    = time_block_0;   //
                out_ind   = prev_out_ind+1; // Index for the next output
                step_iter = 0;              //

                // Particle position
                lon0 = lon_tmps[Ip];
                lat0 = lat_tmps[Ip];

                // Check if initial positions are fill_value, and recycle if they are
                if ( (lon0 == constants::fill_value) or (lat0 == constants::fill_value) ) {
                    recycle_position( lon0, lat0, lon, lat, get_rand(gen), get_rand(gen) );
                }

                // Seed values for velocities (only used for dt)
                particles_get_edges(left, right, bottom, top, lat0, lon0, lat, lon);
                vel_lon_part = particles_interp_from_edges( lat0, lon0, lat, lon, &vel_lon_0, 
                        mask, left, right, bottom, top );

                vel_lat_part = particles_interp_from_edges( lat0, lon0, lat, lon, &vel_lat_0, 
                        mask, left, right, bottom, top );

                // Get initial values for tracked fields
                if ( next_load_index == 1 ) {
                    index = Index(0,       0,      0,      Ip,
                                  Ntime,   Ndepth, Nouts,  Nparts);
                    field_trajectories.at(0).at(index) = vel_lon_part;
                    field_trajectories.at(1).at(index) = vel_lat_part;
                }

                while (t_part < time_block_1) {

                    // Get local dt
                    //   we'll use the previous velocities, which should
                    //   be fine, since it doesn't change very quickly
                    dx_loc = dlon * constants::R_earth * cos(lat0);
                    dy_loc = dlat * constants::R_earth;

                    dt_max = cfl_max * std::min( dx_loc / std::max(vel_lon_part, 1e-3), 
                                                 dy_loc / std::max(vel_lat_part, 1e-3) );
                    dt_min = cfl_min * std::min( dx_loc / std::max(vel_lon_part, 1e-3), 
                                                 dy_loc / std::max(vel_lat_part, 1e-3) );
                    if (dt == 0) { dt = dt_min; } // start with dt_min, and adapt from there

                    // Adjust dt to avoid stepping too far passed an output
                    if ( ( out_ind < Nouts ) and ( (t_part + dt) > target_times.at(out_ind) + 0.01*dt_target ) )
                    {
                        dt_max = target_times.at(out_ind) - t_part;
                        dt = dt_max;
                        dt_min = dt_max;
                    }

                    //
                    //// RKF45 Time-stepping
                    //
                    /*
                    RKF45_update( lon_new, lat_new, dt_new,
                            t_part, lon0, lat0, dt, dt_max, dt_min,
                            vel_lon_0, vel_lon_1, vel_lat_0, vel_lat_1,
                            time_block_0, time_block_1, lat, lon, mask
                            );
                    lon0 = lon_new;
                    lat0 = lat_new;
                    dt = dt_new;
                    t_part += dt;
                    */

                    //
                    //// Simplectide Euler time-stepping
                    //
                    /*
                    SimplecticEuler_update( lon_new, lat_new, dt_new,
                            t_part, lon0, lat0, dt_min, dt_min, dt_min,
                            vel_lon_0, vel_lon_1, vel_lat_0, vel_lat_1,
                            time_block_0, time_block_1, lat, lon, mask
                            );
                    lon0 = lon_new;
                    lat0 = lat_new;
                    dt = dt_new;
                    t_part += dt;
                    */

                    //
                    //// DOPRI5 Time-stepping
                    //
                    DOPRI5_update( lon_new, lat_new, dt_new,
                            t_part, lon0, lat0, dt, dt_max, dt_min,
                            vel_lon_0, vel_lon_1, vel_lat_0, vel_lat_1,
                            time_block_0, time_block_1, lat, lon, mask
                            );
                    lon0 = lon_new;
                    lat0 = lat_new;
                    dt = dt_new;
                    t_part += dt;

                    /*
                    // If our particle went out of bounds, just recycle now.
                    if ( (lat0 <= lat.front()) or (lat0 >= lat.back()) ) {
                        recycle_position( lon0, lat0, lon, lat, get_rand(gen), get_rand(gen) );

                        // We also need to flag the recycle in the previous output
                        if (out_ind > 0) { 
                            index = Index(0,       0,      out_ind-1, Ip,
                                          Ntime,   Ndepth, Nouts,     Nparts);

                            part_lon_hist.at(index) = constants::fill_value;
                            part_lat_hist.at(index) = constants::fill_value;
                        }
                    }
                    */

                    // Keep track of if we're going to recycle at the next output
                    //  we're syncing recycling with outputs so that we can flag
                    //  them as singleton nan values
                    do_recycle = do_recycle or check_if_recycling( dt, t_part, recycle_times[Ip], particle_lifespan );


                    // Track, if at right time
                    if ( (t_part >= target_times.at(out_ind)) and (out_ind < Nouts) ){

                        index = Index(0,       0,      out_ind, Ip,
                                      Ntime,   Ndepth, Nouts,   Nparts);

                        if ( do_recycle ) {
                            // Flag the recycle with fill_value
                            part_lon_hist.at(index) = constants::fill_value;
                            part_lat_hist.at(index) = constants::fill_value;
                            if ( constants::PARTICLE_RECYCLE_TYPE == constants::ParticleRecycleType::FixedInterval ) {
                                recycle_times[Ip] += particle_lifespan; 
                            } else if (constants::PARTICLE_RECYCLE_TYPE == constants::ParticleRecycleType::Stochastic) {
                                recycle_times[Ip] += -log( get_rand(gen) ) * particle_lifespan; 
                            }

                            // And recycle
                            recycle_position( lon0, lat0, lon, lat, get_rand(gen), get_rand(gen) );

                            do_recycle = false;
                        } else {
                            part_lon_hist.at(index) = lon0;
                            part_lat_hist.at(index) = lat0;
                        }

                        up_0 = get_at_point( t_part, lon0, lat0, 
                                vel_lon_0, vel_lon_1, time_block_0, time_block_1, lat, lon, mask );
                        vp_0 = get_at_point( t_part, lon0, lat0, 
                                vel_lat_0, vel_lat_1, time_block_0, time_block_1, lat, lon, mask );

                        field_trajectories.at(0).at(index) = up_0;
                        field_trajectories.at(1).at(index) = vp_0;

                        out_ind++;
                    }

                    // If we've somehow gone too far in time, stop.
                    if ( (t_part >= target_times.back()) or (out_ind >= Nouts) ) { break; }

                    step_iter++;
                } // close inner time loop
                lon_tmps[Ip] = lon0;
                lat_tmps[Ip] = lat0;
            } // close particle loop
        } // close pragma 
        prev_out_ind = out_ind-1;
        #if DEBUG >= 1
        if (wRank == 0) {
            fprintf( stdout, "Finished time-block %zu. Previous out index is %d.\n", 
                    next_load_index, prev_out_ind );
            fflush( stdout );
        }
        #endif
    } // close time-block loop
}
