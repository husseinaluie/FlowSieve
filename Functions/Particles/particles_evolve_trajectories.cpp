#include <math.h>
#include <algorithm>
#include <vector>
#include <mpi.h>
#include <omp.h>
#include "../../constants.hpp"
#include "../../functions.hpp"
#include "../../particles.hpp"

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
        const std::vector<double> & vel_lon,
        const std::vector<double> & vel_lat,
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

    double t_part, lon0, lat0,
           dx_loc, dy_loc, dt,
           vel_lon_part, vel_lat_part, field_val,
           time_p;

    const double dlon = lon.at(1) - lon.at(0),
                 dlat = lat.at(1) - lat.at(0),
                 cfl  = 1e-5,
                 U0   = 2.,
                 dt_target = target_times.at(1) - target_times.at(0);

    const unsigned int  Nlat   = lat.size(),
                        Nlon   = lon.size(),
                        Ntime  = time.size(),
                        Ndepth = 1,
                        Nparts = starting_lat.size(),
                        Nouts  = target_times.size();

    int left, right, bottom, top;

    std::vector<double> vel_lon_sub( Nlat * Nlon ), 
                        vel_lat_sub( Nlat * Nlon );
    std::vector<bool>   mask_sub(    Nlat * Nlon );

    unsigned int out_ind, step_iter, ref_ind, Ip;
    size_t index;

    #pragma omp parallel \
    default(none) \
    shared( lat, lon, vel_lon, vel_lat, mask, stdout,\
            starting_lon, starting_lat,\
            target_times, time, part_lon_hist, part_lat_hist,\
            rev_part_lon_hist, rev_part_lat_hist,\
            field_trajectories, rev_field_trajectories, fields_to_track,\
            wRank, wSize)\
    private(Ip, index, \
            t_part, out_ind, step_iter, ref_ind, lon0, lat0, \
            dx_loc, dy_loc, dt, time_p, vel_lon_part, vel_lat_part, field_val,\
            left, right, bottom, top)
    {
        #pragma omp for collapse(1) schedule(dynamic)
        for (Ip = 0; Ip < Nparts; ++Ip) {
    
            // Particle time
            t_part    = time.at(0);     //
            out_ind   = 0;              //
            step_iter = 0;              //
            ref_ind   = 0;              // time index in source velocity 

            // Particle position
            lon0 = starting_lon.at(Ip);
            lat0 = starting_lat.at(Ip);

            // Get initial values for tracked fields
            index = Index(0,       0,      out_ind, Ip,
                          Ntime,   Ndepth, Nouts,   Nparts);
            part_lon_hist.at(index) = lon0;
            part_lat_hist.at(index) = lat0;

            time_p =    ( t_part             - time.at(ref_ind) ) 
                      / ( time.at(ref_ind+1) - time.at(ref_ind) );

            particles_get_edges(left, right, bottom, top, lat0, lon0, lat, lon);
            for (size_t Ifield = 0; Ifield < fields_to_track.size(); ++Ifield) {
                field_val = particles_interp_from_edges(lat0, lon0, lat, lon, 
                        fields_to_track.at(Ifield), mask,
                        left, right, bottom, top, time_p, ref_ind, Ntime);
                field_trajectories.at(Ifield).at(index) = field_val;
            }
            out_ind++;

            // Seed values for velocities (only used for dt)
            vel_lon_part = U0;
            vel_lat_part = U0;

            #if DEBUG >= 2
            fprintf(stdout, "Particle %03d of %03d (rank %d of %d)\n", 
                    Ip+1 + Nparts * wRank, Nparts * wSize, wRank + 1, wSize);
            fflush(stdout);
            #endif

            while (t_part < target_times.back()) {

                // Get local dt
                //   we'll use the previous velocities, which should
                //   be fine, since it doesn't change very quickly
                dx_loc = dlon * constants::R_earth * cos(lat0);
                dy_loc = dlat * constants::R_earth;

                dt = cfl * std::min( dx_loc / std::max(vel_lon_part, 1e-3), 
                                     dy_loc / std::max(vel_lat_part, 1e-3) );

                if ( ( (size_t)out_ind < target_times.size() )
                     and ( (t_part + dt) - target_times.at(out_ind) > (dt_target / 50.) ) 
                   )
                {
                    dt = target_times.at(out_ind) - t_part;
                }

                // Subset velocities by time
                time_p =    ( t_part             - time.at(ref_ind) ) 
                          / ( time.at(ref_ind+1) - time.at(ref_ind) );

                //
                //// Get u_lon at position at advance lon position
                //
                particles_get_edges(left, right, bottom, top, lat0, lon0, lat, lon);
                if ( (bottom < 0) or (top < 0) ) { break; }
                vel_lon_part = particles_interp_from_edges(lat0, lon0, lat, lon, &vel_lon, 
                        mask, left, right, bottom, top, time_p, ref_ind, Ntime);
                if ( fabs(vel_lon_part) > 100. ) { break; }

                // convert to radial velocity and step in space
                lon0 += dt * vel_lon_part / (constants::R_earth * cos(lat0));
                if (lon0 >  M_PI) { lon0 -= 2 * M_PI; }
                if (lon0 < -M_PI) { lon0 += 2 * M_PI; }

                //
                //// Get u_lat at position at advance lat position
                //
                particles_get_edges(left, right, bottom, top, lat0, lon0, lat, lon);
                if ( (bottom < 0) or (top < 0) ) { break; }
                vel_lat_part = particles_interp_from_edges(lat0, lon0, lat, lon, &vel_lat, 
                        mask, left, right, bottom, top, time_p, ref_ind, Ntime);
                if ( fabs(vel_lat_part) > 100. ) { break; }

                // convert to radial velocity and step in space
                lat0 += dt * vel_lat_part / constants::R_earth;

                // Update time
                t_part += dt;

                // Track, if at right time
                if (      (t_part < target_times.back()) 
                      and (t_part >= target_times.at(out_ind)) 
                   ){

                    index = Index(0,       0,      out_ind, Ip,
                                  Ntime,   Ndepth, Nouts,   Nparts);
                    part_lon_hist.at(index) = lon0;
                    part_lat_hist.at(index) = lat0;

                    particles_get_edges(left, right, bottom, top, lat0, lon0, lat, lon);

                    for (size_t Ifield = 0; Ifield < fields_to_track.size(); ++Ifield) {
                        field_val = particles_interp_from_edges(lat0, lon0, lat, lon, 
                                fields_to_track.at(Ifield), mask,
                                left, right, bottom, top, time_p, ref_ind, Ntime);
                        field_trajectories.at(Ifield).at(index) = field_val;
                    }

                    out_ind++;
                }

                if ( (t_part >= target_times.back()) or (out_ind >= Nouts) ) { break; }

                if (t_part > time.at(ref_ind+1)) { ref_ind++; }
                step_iter++;
            }
            #if DEBUG >= 2
            fprintf(stdout, "  time %g / %g (%g %% - %d)\n", 
                    t_part, target_times.back(), 100 * t_part / target_times.back(), step_iter);
            fflush(stdout);
            #endif


            /*
            #if DEBUG >= 1
            //
            //// Reversing the trajectory
            //
            out_ind--;

            // Restart particle position
            //    to evolve backwards
            lon0 = starting_lon.at(Ip);
            lat0 = starting_lat.at(Ip);

            while ( (t_part > 0) and (out_ind >= 0) and (ref_ind >= 0) ) {

                // Get local dt
                //   we'll use the previous velocities, which should
                //   be fine, since it doesn't change very quickly
                dx_loc = dlon * constants::R_earth * cos(lat0);
                dy_loc = dlat * constants::R_earth;

                dt = cfl * std::min( dx_loc / std::max(vel_lon_part, 1e-3), 
                                     dy_loc / std::max(vel_lat_part, 1e-3) );
                dt *= -1.; // negate for reverse-stepping

                if ( target_times[out_ind] - (t_part + dt) > (dt_target / 50.) ) {
                    dt = target_times[out_ind] - t_part;
                }

                // Subset velocities by time
                time_p =    ( t_part             - time.at(ref_ind) ) 
                          / ( time.at(ref_ind+1) - time.at(ref_ind) );

                //
                //// Get u_lon at position at advance lon position
                //
                particles_get_edges(left, right, bottom, top, lat0, lon0, lat, lon);
                if ( (bottom < 0) or (top < 0) ) { break; }
                vel_lon_part = particles_interp_from_edges(lat0, lon0, lat, lon, &vel_lon, 
                        mask, left, right, bottom, top, time_p, ref_ind, Ntime);

                // convert to radial velocity and step in space
                lon0 += dt * vel_lon_part / (constants::R_earth * cos(lat0));
                if (lon0 >  M_PI) { lon0 -= 2 * M_PI; }
                if (lon0 < -M_PI) { lon0 += 2 * M_PI; }

                //
                //// Get u_lat at position at advance lat position
                //
                particles_get_edges(left, right, bottom, top, lat0, lon0, lat, lon);
                if ( (bottom < 0) or (top < 0) ) { break; }
                vel_lat_part = particles_interp_from_edges(lat0, lon0, lat, lon, &vel_lat, 
                        mask, left, right, bottom, top, time_p, ref_ind, Ntime);

                // convert to radial velocity and step in space
                lat0 += dt * vel_lat_part / constants::R_earth;

                // Update time
                t_part += dt;

                // Track, if at right time
                if (t_part <= target_times[out_ind]) {
                    //fprintf(stdout, "  t, out_ind = (%.3g, %d)\n", t_part, out_ind);
                    index = Index(0,       0,      out_ind, Ip,
                                  Ntime,   Ndepth, Nouts,   Nparts);

                    rev_part_lon_hist.at(index) = lon0;
                    rev_part_lat_hist.at(index) = lat0;

                    particles_get_edges(left, right, bottom, top, lat0, lon0, lat, lon);
                    for (size_t Ifield = 0; Ifield < fields_to_track.size(); ++Ifield) {
                        field_val = particles_interp_from_edges(lat0, lon0, lat, lon, 
                                fields_to_track.at(Ifield), mask,
                                left, right, bottom, top, time_p, ref_ind, Ntime);
                        rev_field_trajectories.at(Ifield).at(index) = field_val;
                    }

                    out_ind--;
                }

                if (t_part < time.at(ref_ind)) { ref_ind--; }
                step_iter++;
            }
            #endif
            */
        } 
    }
}
