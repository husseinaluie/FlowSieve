#include <math.h>
#include <algorithm>
#include <vector>
#include <omp.h>
#include "../../constants.hpp"
#include "../../functions.hpp"
#include "../../particles.hpp"

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
        ) {

    double t_part, lon0, lat0,
           dx_loc, dy_loc, dt,
           vel_lon_part, vel_lat_part,
           time_p;

    const double dlon = lon.at(1) - lon.at(0),
                 dlat = lat.at(1) - lat.at(0),
                 cfl  = 0.005,
                 U0   = 2.;

    const int Nlat   = lat.size(),
              Nlon   = lon.size(),
              Ntime  = time.size(),
              Ndepth = 1,
              Nparts = starting_lat.size(),
              Nouts  = target_times.size();

    int left, right, bottom, top;

    std::vector<double> vel_lon_sub( Nlat * Nlon ), 
                        vel_lat_sub( Nlat * Nlon ),
                        mask_sub(    Nlat * Nlon );

    int out_ind, step_iter, ref_ind, 
        Itime, index, Ip;

    #pragma omp parallel \
    default(none) \
    shared( lat, lon, vel_lon, vel_lat, mask, stdout,\
            starting_lon, starting_lat,\
            target_times, time, part_lon_hist, part_lat_hist)\
    private(Ip, Itime, index, \
            t_part, out_ind, step_iter, ref_ind, lon0, lat0, \
            dx_loc, dy_loc, dt, time_p, vel_lon_part, vel_lat_part,\
            left, right, bottom, top)
    {
        #pragma omp for collapse(1) schedule(dynamic)
        for (Ip = 0; Ip < Nparts; ++Ip) {
    
            // Particle time
            t_part    = time.at(0);
            out_ind   = 1;
            step_iter = 0;
            ref_ind   = 0;

            //
            Itime = 0;

            // Particle position
            lon0 = starting_lon.at(Ip);
            lat0 = starting_lat.at(Ip);

            // Seed values for velocities (only used for dt)
            vel_lon_part = U0;
            vel_lat_part = U0;

            fprintf(stdout, "Particle %03d of %03d\n", Ip+1, Nparts);
            fflush(stdout);

            while (t_part < target_times.back()) {

                if (t_part >= time.at(Itime)) { Itime++; }

                // Get local dt
                //   we'll use the previous velocities, which should
                //   be fine, since it doesn't change very quickly
                dx_loc = dlon * constants::R_earth * cos(lat0);
                dy_loc = dlat * constants::R_earth;

                dt = cfl * std::min( dx_loc / std::max(vel_lon_part, 1e-3), 
                                     dy_loc / std::max(vel_lat_part, 1e-3) );
                //fprintf(stdout, "  (dt, dx, dy)       = (%g, %g, %g)\n", dt, dx_loc, dy_loc);
                //fprintf(stdout, "      (u_lon, u_lat) = (%g, %g)\n", vel_lon_part, vel_lat_part);
                //fflush(stdout);

                // Subset velocities by time
                time_p =    ( t_part             - time.at(ref_ind) ) 
                          / ( time.at(ref_ind+1) - time.at(ref_ind) );

                //
                //// Get u_lon at position at advance lon position
                //
                particles_get_edges(left, right, bottom, top, lat0, lon0, lat, lon);
                if ( (bottom < 0) or (top < 0) ) { break; }
                vel_lon_part = particles_interp_from_edges(lat0, lon0, lat, lon, vel_lon, 
                        mask, left, right, bottom, top, time_p, Itime, Ntime);

                // convert to radial velocity and step in space
                lon0 += dt * vel_lon_part / (constants::R_earth * cos(lat0));
                if (lon0 >  M_PI) { lon0 -= 2 * M_PI; }
                if (lon0 < -M_PI) { lon0 += 2 * M_PI; }

                //
                //// Get u_lat at position at advance lat position
                //
                particles_get_edges(left, right, bottom, top, lat0, lon0, lat, lon);
                if ( (bottom < 0) or (top < 0) ) { break; }
                vel_lat_part = particles_interp_from_edges(lat0, lon0, lat, lon, vel_lat, 
                        mask, left, right, bottom, top, time_p, Itime, Ntime);

                // convert to radial velocity and step in space
                lat0 += dt * vel_lat_part / constants::R_earth;

                // Update time
                t_part += dt;
                if (t_part > time.at(ref_ind+1)) { ref_ind++; }

                // Track, if at right time
                if (      (t_part < target_times.back()) 
                        and (t_part + dt >= target_times[out_ind]) 
                   ){
                    index = Index(0,       0,      out_ind, Ip,
                                  Ntime,   Ndepth, Nouts,   Nparts);
                    part_lon_hist.at(index) = lon0;
                    part_lat_hist.at(index) = lat0;

                    out_ind++;
                }

                step_iter++;
            }
        } 
    }
}
