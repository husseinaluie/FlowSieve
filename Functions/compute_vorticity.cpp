#include <vector>
#include <omp.h>
#include "../functions.hpp"
#include "../constants.hpp"

/*!
 * \brief Wrapper for computing vorticity
 *
 * This wrapper applies compute_vorticity_at_point() at each point in the
 * grid. If the compiler flag COMP_VORT is set to false, then this procedure
 * is skipped.
 *
 *  @param[in,out]      vort_r,vort_lon,vort_lat    where to store computed vorticity components (array)
 *  @param[in,out]      vel_div                     where to store computed velocity divergence (array)
 *  @param[in,out]      OkuboWeiss                  where to store computed OkuboWeiss (array)
 *  @param[in]          source_data                 dataset class instance containing data (Psi, Phi, etc)
 *  @param[in]          u_r,u_lon,u_lat             velocity components
 *  @param[in]          mask                        2D array to distinguish land from water
 *  @param[in]          comm                        MPI communicator object
 *
 */
void compute_vorticity(
        std::vector<double> & vort_r,
        std::vector<double> & vort_lon,
        std::vector<double> & vort_lat,
        std::vector<double> & vel_div,
        std::vector<double> & OkuboWeiss,
        std::vector<double> & cyclonic_energy,
        std::vector<double> & anticyclonic_energy,
        std::vector<double> & divergent_strain_energy,
        std::vector<double> & traceless_strain_energy,
        const dataset & source_data,
        const std::vector<double> & u_r,
        const std::vector<double> & u_lon,
        const std::vector<double> & u_lat,
        const MPI_Comm comm
        ) {

    const std::vector<bool> &mask = source_data.mask;

    const int   Ntime   = source_data.Ntime,
                Ndepth  = source_data.Ndepth,
                Nlat    = source_data.Nlat,
                Nlon    = source_data.Nlon;

    double vort_r_tmp, vort_lon_tmp, vort_lat_tmp, div_tmp, OkuboWeiss_tmp,
           cyclonic_energy_tmp, anticyclonic_energy_tmp, 
           divergent_strain_energy_tmp, traceless_strain_energy_tmp;
    int Itime, Idepth, Ilat, Ilon;
    size_t index; 
    const size_t Npts = u_lon.size();

    // If any of the 'output' arrays are size zero 
    // don't do them (this is essentially how to 'turn off' outputs)
    const bool do_vort_r   = vort_r.size() > 0;
    const bool do_vort_lon = vort_lon.size() > 0;
    const bool do_vort_lat = vort_lat.size() > 0;

    const bool do_vel_div = vel_div.size() > 0;

    const bool do_OkuboWeiss = OkuboWeiss.size() > 0;

    const bool  do_cyclonic_energy = cyclonic_energy.size() > 0,
                do_anticyclonic_energy = anticyclonic_energy.size() > 0,
                do_divergent_strain_energy = divergent_strain_energy.size() > 0,
                do_traceless_strain_energy = traceless_strain_energy.size() > 0;

    #if DEBUG >= 2
    int wRank, wSize;
    MPI_Comm_rank( comm, &wRank );
    MPI_Comm_size( comm, &wSize );

    if (wRank == 0) { fprintf(stdout, "  Starting vorticity computation.\n"); }
    fflush(stdout);
    #endif

    #pragma omp parallel \
    default(none) \
    shared( source_data, mask, u_r, u_lon, u_lat, \
            vort_r, vort_lon, vort_lat, vel_div, OkuboWeiss, \
            cyclonic_energy, anticyclonic_energy, \
            divergent_strain_energy, traceless_strain_energy ) \
    private( Itime, Idepth, Ilat, Ilon, index, vort_r_tmp, vort_lon_tmp, vort_lat_tmp, \
             div_tmp, OkuboWeiss_tmp, cyclonic_energy_tmp, \
            anticyclonic_energy_tmp, \
            divergent_strain_energy_tmp, traceless_strain_energy_tmp )\
    firstprivate( Npts, Nlon, Nlat, Ndepth, Ntime, do_vort_r, do_vort_lon, do_vort_lat, do_vel_div,\
                  do_OkuboWeiss, do_cyclonic_energy, do_anticyclonic_energy, \
                  do_divergent_strain_energy, do_traceless_strain_energy )
    {
        #pragma omp for collapse(1) schedule(guided)
        for (index = 0; index < Npts; index++) {

            vort_r_tmp   = 0.; 
            vort_lon_tmp = 0.; 
            vort_lat_tmp = 0.; 

            div_tmp = 0.; 

            OkuboWeiss_tmp = 0.; 

            if ( mask.at(index) ) { // Skip land areas

                Index1to4(index, Itime, Idepth, Ilat, Ilon,
                                 Ntime, Ndepth, Nlat, Nlon);

                compute_vorticity_at_point(
                        vort_r_tmp, vort_lon_tmp, vort_lat_tmp, div_tmp, OkuboWeiss_tmp, 
                        cyclonic_energy_tmp, anticyclonic_energy_tmp, 
                        divergent_strain_energy_tmp, traceless_strain_energy_tmp,
                        source_data, u_r, u_lon, u_lat,        
                        Itime, Idepth, Ilat, Ilon);
            }

            if (do_vort_r)   { vort_r.at(  index) = vort_r_tmp; }
            if (do_vort_lon) { vort_lon.at(index) = vort_lon_tmp; }
            if (do_vort_lat) { vort_lat.at(index) = vort_lat_tmp; }

            if (do_vel_div) { vel_div.at(index) = div_tmp; }

            if (do_OkuboWeiss) { OkuboWeiss.at(index) = OkuboWeiss_tmp; }

            if (do_cyclonic_energy) { cyclonic_energy.at(index) = cyclonic_energy_tmp; }
            if (do_anticyclonic_energy) { anticyclonic_energy.at(index) = anticyclonic_energy_tmp; }
            if (do_divergent_strain_energy) { divergent_strain_energy.at(index) = divergent_strain_energy_tmp; }
            if (do_traceless_strain_energy) { traceless_strain_energy.at(index) = traceless_strain_energy_tmp; }


        } // end index loop
    } // end pragma
    #if DEBUG >= 2
    if (wRank == 0) { fprintf(stdout, "     ... done.\n"); }
    #endif
} // end compute_vorticity
