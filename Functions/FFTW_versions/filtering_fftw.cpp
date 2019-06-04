#include <math.h>
#include <algorithm>
#include <vector>
#include <omp.h>
#include <mpi.h>
#include "../../functions.hpp"
#include "../../fft_based.hpp"
#include "../../netcdf_io.hpp"
#include "../../constants.hpp"
#include <fftw3.h>
#include <assert.h>

void filtering_fftw(
        const std::vector<double> & full_u_z,       /**< [in] Full u_r velocity array */
        const std::vector<double> & full_u_x,       /**< [in] Full u_lon velocity array */
        const std::vector<double> & full_u_y,       /**< [in] Full u_lat velocity array */
        const std::vector<double> & full_rho,       /**< [in] Full density array */
        const std::vector<double> & full_p,         /**< [in] Full pressure array */
        const std::vector<double> & scales,         /**< [in] Array of filtering scales */
        const std::vector<double> & dAreas,         /**< [in] Array of cell areas (2D) (compute_areas()) */
        const std::vector<double> & time,           /**< [in] Time dimension (1D) */
        const std::vector<double> & depth,          /**< [in] Depth dimension (1D) */
        const std::vector<double> & longitude,      /**< [in] Longitude dimension (1D) */
        const std::vector<double> & latitude,       /**< [in] Latitude dimension (1D) */
        const std::vector<double> & mask,           /**< [in] Array to distinguish between land and water cells (2D) */
        const std::vector<int>    & myCounts,       /**< [in] Array of dimension sizes */
        const std::vector<int>    & myStarts,       /**< [in] Array of dimension sizes */
        const MPI_Comm comm                         /**< [in] MPI Communicator */
        ) {

    assert(constants::CARTESIAN);
    assert(constants::PERIODIC_X);
    assert(constants::PERIODIC_Y);

    int wRank, wSize;
    MPI_Comm_rank( comm, &wRank );
    MPI_Comm_size( comm, &wSize );

    // Get dimension sizes
    const int Nscales = scales.size();
    const int Ntime   = myCounts.at(0);
    const int Ndepth  = myCounts.at(1);
    const int Nlat    = myCounts.at(2);
    const int Nlon    = myCounts.at(3);

    const double Llat = Nlat * (latitude.at( 1) - latitude.at( 0));
    const double Llon = Nlon * (longitude.at(1) - longitude.at(0));

    const unsigned int num_pts = Ntime * Ndepth * Nlat * Nlon;
    char filename [50];
    
    const int ndims = 4;
    size_t starts[ndims] = {
        size_t(myStarts.at(0)), size_t(myStarts.at(1)), 
        size_t(myStarts.at(2)), size_t(myStarts.at(3))};
    size_t counts[ndims] = {
        size_t(Ntime), size_t(Ndepth), 
        size_t(Nlat), size_t(Nlon)};
    std::vector<std::string> vars_to_write;

    #if DEBUG >= 1
    if (wRank == 0) { fprintf(stdout, "Converting to Cartesian velocities.\n"); }
    #endif

    std::vector<double> local_kernel(Nlat * Nlon);

    std::vector<double> coarse_u_x(num_pts);
    std::vector<double> coarse_u_y(num_pts);
    std::vector<double> coarse_u_z(num_pts);
    vars_to_write.push_back("coarse_u_r");
    vars_to_write.push_back("coarse_u_lon");
    vars_to_write.push_back("coarse_u_lat");

    std::vector<double> full_KE(num_pts);

    size_t index;
    #pragma omp parallel default(none) private(index) \
        shared(full_KE, full_u_x, full_u_y, full_u_z)
    {
        #pragma omp for collapse(1) schedule(static)
        for (index = 0; index < num_pts; ++index) {
            full_KE.at(index) = 
                0.5 * ( 
                          pow(full_u_x.at(index), 2) 
                        + pow(full_u_y.at(index), 2) 
                        + pow(full_u_z.at(index), 2) 
                      );
        } 
    } // done pragma

    double scale;

    std::vector<double> fine_u_x(num_pts);
    std::vector<double> fine_u_y(num_pts);
    std::vector<double> fine_u_z(num_pts);
    vars_to_write.push_back("fine_u_r");
    vars_to_write.push_back("fine_u_lon");
    vars_to_write.push_back("fine_u_lat");

    std::vector<double> div_J(num_pts);
    vars_to_write.push_back("div_Jtransport");

    std::vector<double> fine_KE(num_pts);
    std::vector<double> coarse_KE(num_pts);
    vars_to_write.push_back("fine_KE");
    vars_to_write.push_back("coarse_KE");

    // If computing energy transfers, we'll need some more arrays
    std::vector<double> coarse_uxux(num_pts);
    std::vector<double> coarse_uxuy(num_pts);
    std::vector<double> coarse_uxuz(num_pts);
    std::vector<double> coarse_uyuy(num_pts);
    std::vector<double> coarse_uyuz(num_pts);
    std::vector<double> coarse_uzuz(num_pts);
    vars_to_write.push_back("coarse_uxux");
    vars_to_write.push_back("coarse_uxuy");
    vars_to_write.push_back("coarse_uxuz");
    vars_to_write.push_back("coarse_uyuy");
    vars_to_write.push_back("coarse_uyuz");
    vars_to_write.push_back("coarse_uzuz");

    /*
    std::vector<double> tau_xx(num_pts);
    std::vector<double> tau_xy(num_pts);
    std::vector<double> tau_xz(num_pts);
    std::vector<double> tau_yy(num_pts);
    std::vector<double> tau_yz(num_pts);
    std::vector<double> tau_zz(num_pts);
    vars_to_write.push_back("tau_xx");
    vars_to_write.push_back("tau_xy");
    vars_to_write.push_back("tau_xz");
    vars_to_write.push_back("tau_yy");
    vars_to_write.push_back("tau_yz");
    vars_to_write.push_back("tau_zz");

    double s_xx, s_xy, s_xz, s_yy, s_yz, s_zz;
    double ux_loc, uy_loc, uz_loc, pi_tmp;
    std::vector<double> S_xx(num_pts);
    std::vector<double> S_xy(num_pts);
    std::vector<double> S_xz(num_pts);
    std::vector<double> S_yy(num_pts);
    std::vector<double> S_yz(num_pts);
    std::vector<double> S_zz(num_pts);
    vars_to_write.push_back("S_xx");
    vars_to_write.push_back("S_xy");
    vars_to_write.push_back("S_xz");
    vars_to_write.push_back("S_yy");
    vars_to_write.push_back("S_yz");
    vars_to_write.push_back("S_zz");
    */

    std::vector<double> div(num_pts);
    vars_to_write.push_back("full_vel_div");
    vars_to_write.push_back("coarse_vel_div");

    // Also an array for the transfer itself
    std::vector<double> energy_transfer(num_pts);
    vars_to_write.push_back("energy_transfer");

    // Placeholder
    std::vector<double> coarse_p;

    // Initialize fftw plans
    const int rank = 2; 
    const int howmany = Ntime;
    const int idist = Nlat * Nlon; 
    const int odist = Nlat * (Nlon/2 + 1);
    const int istride = 1; 
    const int ostride = 1;
    int n[] = {Nlat, Nlon};
    int *inembed = NULL, *onembed = NULL;
    const int Nk = Ntime * Ndepth * Nlat * (Nlon/2 + 1);

    #if DEBUG >= 1
    if (wRank == 0) { fprintf(stdout, "Preparing FFTW plans.\n"); }
    #endif
    double *fft_inp = fftw_alloc_real(Ntime * Ndepth * Nlat * Nlon);

    #if DEBUG >= 1
    if (wRank == 0) { fprintf(stdout, "   ... complex storage\n"); }
    #endif
    fftw_complex *fft_out = fftw_alloc_complex(Nk);

    #if DEBUG >= 1
    if (wRank == 0) { fprintf(stdout, "   ... forward plan\n"); }
    #endif
    fftw_plan fft = fftw_plan_many_dft_r2c(
            rank, n, howmany,
            fft_inp, inembed, istride, idist,
            fft_out, onembed, ostride, odist,
            FFTW_PATIENT);

    #if DEBUG >= 1
    if (wRank == 0) { fprintf(stdout, "   ... backward plan\n"); }
    #endif
    fftw_plan ifft = fftw_plan_many_dft_c2r(
            rank, n, howmany,
            fft_out, onembed, ostride, odist,
            fft_inp, inembed, istride, idist,
            FFTW_PATIENT);

    //
    //// Begin the main filtering loop
    //
    #if DEBUG>=1
    if (wRank == 0) { fprintf(stdout, "Beginning main filtering loop.\n\n"); }
    #endif
    for (int Iscale = 0; Iscale < Nscales; Iscale++) {

        // Create the output file
        snprintf(filename, 50, "filter_%.6gkm.nc", scales.at(Iscale)/1e3);
        initialize_output_file(time, depth, longitude, latitude, 
                mask, vars_to_write, filename, scales.at(Iscale));

        #if DEBUG >= 0
        if (wRank == 0) { 
            fprintf(stdout, "Scale %d of %d (%.5g km)\n", 
                Iscale+1, Nscales, scales.at(Iscale)/1e3); 
        }
        #endif

        scale = scales.at(Iscale);

        // Filter velocity fields
        fft_filter(coarse_u_x, full_u_x, fft, ifft, fft_inp, fft_out, scale,
                Ntime, Ndepth, Nlat, Nlon, Llat, Llon);
        fft_filter(coarse_u_y, full_u_y, fft, ifft, fft_inp, fft_out, scale,
                Ntime, Ndepth, Nlat, Nlon, Llat, Llon);
        fft_filter(coarse_u_z, full_u_z, fft, ifft, fft_inp, fft_out, scale,
                Ntime, Ndepth, Nlat, Nlon, Llat, Llon);

        // Also filter KE
        fft_filter(coarse_KE, full_KE, fft, ifft, fft_inp, fft_out, scale,
                Ntime, Ndepth, Nlat, Nlon, Llat, Llon);

        for (size_t index = 0; index < num_pts; ++index) {
            fine_u_x.at(index) = full_u_x.at(index) - coarse_u_x.at(index);
            fine_u_y.at(index) = full_u_y.at(index) - coarse_u_y.at(index);
            fine_u_z.at(index) = full_u_z.at(index) - coarse_u_z.at(index);
            fine_KE.at(index) = full_KE.at(index) - coarse_KE.at(index);
        }

        // Filter Quadratics (for tau)
        fft_filter_product(coarse_uxux, full_u_x, full_u_x, fft, ifft, fft_inp, fft_out,
                scale, Ntime, Ndepth, Nlat, Nlon, Llat, Llon);
        fft_filter_product(coarse_uxuy, full_u_x, full_u_y, fft, ifft, fft_inp, fft_out,
                scale, Ntime, Ndepth, Nlat, Nlon, Llat, Llon);
        fft_filter_product(coarse_uxuz, full_u_x, full_u_z, fft, ifft, fft_inp, fft_out,
                scale, Ntime, Ndepth, Nlat, Nlon, Llat, Llon);
        fft_filter_product(coarse_uyuy, full_u_y, full_u_y, fft, ifft, fft_inp, fft_out,
                scale, Ntime, Ndepth, Nlat, Nlon, Llat, Llon);
        fft_filter_product(coarse_uyuz, full_u_y, full_u_z, fft, ifft, fft_inp, fft_out,
                scale, Ntime, Ndepth, Nlat, Nlon, Llat, Llon);
        fft_filter_product(coarse_uzuz, full_u_z, full_u_z, fft, ifft, fft_inp, fft_out,
                scale, Ntime, Ndepth, Nlat, Nlon, Llat, Llon);

        #if DEBUG >= 1
        fprintf(stdout, "  = Rank %d finished filtering =\n", wRank);
        #endif

        // Write to file
        write_field_to_output(fine_u_x,   "fine_u_lon",   starts, counts, filename);
        write_field_to_output(fine_u_y,   "fine_u_lat",   starts, counts, filename);
        write_field_to_output(fine_u_z,   "fine_u_r",     starts, counts, filename);
        write_field_to_output(coarse_u_x, "coarse_u_lon", starts, counts, filename);
        write_field_to_output(coarse_u_y, "coarse_u_lat", starts, counts, filename);
        write_field_to_output(coarse_u_z, "coarse_u_r",   starts, counts, filename);

        write_field_to_output(coarse_uxux, "coarse_uxux",   starts, counts, filename);
        write_field_to_output(coarse_uxuy, "coarse_uxuy",   starts, counts, filename);
        write_field_to_output(coarse_uxuz, "coarse_uxuz",   starts, counts, filename);
        write_field_to_output(coarse_uyuy, "coarse_uyuy",   starts, counts, filename);
        write_field_to_output(coarse_uyuz, "coarse_uyuz",   starts, counts, filename);
        write_field_to_output(coarse_uzuz, "coarse_uzuz",   starts, counts, filename);

        write_field_to_output(fine_KE,   "fine_KE",   starts, counts, filename);
        write_field_to_output(coarse_KE, "coarse_KE", starts, counts, filename);

        // Compute the energy transfer through the filter scale
        compute_energy_transfer_through_scale(
                energy_transfer, 
                coarse_u_x,  coarse_u_y,  coarse_u_z,
                coarse_uxux, coarse_uxuy, coarse_uxuz,
                coarse_uyuy, coarse_uyuz, coarse_uzuz,
                Ntime, Ndepth, Nlat, Nlon,
                longitude, latitude, mask);

        /*
        for (int Itime = 0; Itime < Ntime; Itime++) {
            for (int Idepth = 0; Idepth < Ndepth; Idepth++) {
                for (int Ilat = 0; Ilat < Nlat; Ilat++) {
                    for (int Ilon = 0; Ilon < Nlon; Ilon++) {
                        index = Index(Itime, Idepth, Ilat, Ilon,
                                      Ntime, Ndepth, Nlat, Nlon);

                        ux_loc = coarse_u_x.at(index);
                        uy_loc = coarse_u_y.at(index);
                        uz_loc = coarse_u_z.at(index);

                        // Compute subfilter-scale stress
                        tau_xx.at(index) = coarse_uxux.at(index) - ux_loc * ux_loc;
                        tau_xy.at(index) = coarse_uxuy.at(index) - ux_loc * uy_loc;
                        tau_xz.at(index) = coarse_uxuz.at(index) - ux_loc * uz_loc;
                        tau_yy.at(index) = coarse_uyuy.at(index) - uy_loc * uy_loc;
                        tau_yz.at(index) = coarse_uyuz.at(index) - uy_loc * uz_loc;
                        tau_zz.at(index) = coarse_uzuz.at(index) - uz_loc * uz_loc;

                        // Compute large-scale strain
                        compute_largescale_strain(
                                s_xx, s_xy, s_xz, s_yy, s_yz, s_zz, 
                                coarse_u_x, coarse_u_y, coarse_u_z,
                                Itime, Idepth, Ilat, Ilon,
                                Ntime, Ndepth, Nlat, Nlon,
                                longitude, latitude, mask);
                
                        S_xx.at(index) = s_xx;
                        S_xy.at(index) = s_xy;
                        S_xz.at(index) = s_xz;
                        S_yy.at(index) = s_yy;
                        S_yz.at(index) = s_yz;
                        S_zz.at(index) = s_zz;

                        pi_tmp = -constants::rho0 * 
                            (   s_xx * tau_xx.at(index) 
                              + s_yy * tau_yy.at(index) 
                              + s_zz * tau_zz.at(index)
                              + 2 * s_xy * tau_xy.at(index) 
                              + 2 * s_xz * tau_xz.at(index) 
                              + 2 * s_yz * tau_yz.at(index) 
                            );

                        energy_transfer.at(index) = pi_tmp;
                    }
                }
            }
        }
        write_field_to_output(tau_xx, "tau_xx", starts, counts, filename);
        write_field_to_output(tau_xy, "tau_xy", starts, counts, filename);
        write_field_to_output(tau_xz, "tau_xz", starts, counts, filename);
        write_field_to_output(tau_yy, "tau_yy", starts, counts, filename);
        write_field_to_output(tau_yz, "tau_yz", starts, counts, filename);
        write_field_to_output(tau_zz, "tau_zz", starts, counts, filename);

        write_field_to_output(S_xx, "S_xx", starts, counts, filename);
        write_field_to_output(S_xy, "S_xy", starts, counts, filename);
        write_field_to_output(S_xz, "S_xz", starts, counts, filename);
        write_field_to_output(S_yy, "S_yy", starts, counts, filename);
        write_field_to_output(S_yz, "S_yz", starts, counts, filename);
        write_field_to_output(S_zz, "S_zz", starts, counts, filename);
        */

        write_field_to_output(energy_transfer, "energy_transfer", starts, counts, filename);

        // Compute the divergence of the coarse field and full field (for comparison)
        compute_div_vel(div, coarse_u_x, coarse_u_y, coarse_u_z, longitude, latitude,
                Ntime, Ndepth, Nlat, Nlon, mask);
        write_field_to_output(div, "coarse_vel_div", starts, counts, filename);

        compute_div_vel(div, full_u_x, full_u_y, full_u_z, longitude, latitude,
                Ntime, Ndepth, Nlat, Nlon, mask);
        write_field_to_output(div, "full_vel_div", starts, counts, filename);

        #if DEBUG >= 1
        if (wRank == 0) { fprintf(stdout, "Starting compute_div_transport\n"); }
        #endif

        compute_div_transport(
                div_J,
                coarse_u_x,  coarse_u_y,  coarse_u_z,
                coarse_uxux, coarse_uxuy, coarse_uxuz,
                coarse_uyuy, coarse_uyuz, coarse_uzuz,
                coarse_p, longitude, latitude,
                Ntime, Ndepth, Nlat, Nlon,
                mask);
        write_field_to_output(div_J, "div_Jtransport", starts, counts, filename);

        #if DEBUG >= 0
        // Flushing stdout is necessary for SLURM outputs.
        fflush(stdout);
        #endif

    }  // end for(scale) block

    // Cleanup some fftw stuff
    fftw_cleanup();
    fftw_destroy_plan(fft);
    fftw_destroy_plan(ifft);
    fftw_free(fft_inp);
    fftw_free(fft_out);

} // end filtering
