#include <math.h>
#include <algorithm>
#include <vector>
#include <omp.h>
#include <mpi.h>
#include "../../functions.hpp"
#include "../../fft_based.hpp"
#include "../../netcdf_io.hpp"
#include "../../constants.hpp"
#include "../../postprocess.hpp"
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
        const std::vector<bool> & mask,             /**< [in] Array to distinguish between land and water cells (2D) */
        const std::vector<int>    & myCounts,       /**< [in] Array of dimension sizes */
        const std::vector<int>    & myStarts,       /**< [in] Array of dimension sizes */
        const MPI_Comm comm                         /**< [in] MPI Communicator */
        ) {

    static_assert( constants::CARTESIAN  );
    static_assert( constants::PERIODIC_X );
    static_assert( constants::PERIODIC_Y );

    int wRank, wSize;
    MPI_Comm_rank( comm, &wRank );
    MPI_Comm_size( comm, &wSize );

    size_t index;
    double scale;

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

    std::vector<double> coarse_u_x(num_pts), coarse_u_y(num_pts), coarse_u_z(num_pts);
    vars_to_write.push_back("coarse_u_x");
    vars_to_write.push_back("coarse_u_y");
    vars_to_write.push_back("coarse_u_z");

    std::vector<double> coarse_KE(num_pts);
    vars_to_write.push_back("coarse_KE");

    std::vector<double> div;
    if (not(constants::MINIMAL_OUTPUT)) {
        div.resize(num_pts);
        vars_to_write.push_back("full_vel_div");
        vars_to_write.push_back("coarse_vel_div");
    }

    // If computing energy transfers, we'll need some more arrays
        std::vector<double> coarse_uxux, coarse_uxuy, coarse_uxuz,
                                         coarse_uyuy, coarse_uyuz, 
                                                      coarse_uzuz;   
        std::vector<double> energy_transfer;
    if (constants::COMP_TRANSFERS) {
        coarse_uxux.resize(num_pts);
        coarse_uxuy.resize(num_pts);
        coarse_uxuz.resize(num_pts);

        coarse_uyuy.resize(num_pts);
        coarse_uyuz.resize(num_pts);

        coarse_uzuz.resize(num_pts);

        // Also an array for the transfer itself
        energy_transfer.resize(num_pts);
        vars_to_write.push_back("energy_transfer");
    }

    // Preset some post-processing variables
    std::vector<const std::vector<double>*> postprocess_fields;
    std::vector<std::string> postprocess_names;

    postprocess_names.push_back("coarse_KE");
    postprocess_fields.push_back(&coarse_KE);

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
        if (not(constants::NO_FULL_OUTPUTS)) {
            snprintf(filename, 50, "filter_%.6gkm.nc", scales.at(Iscale)/1e3);
            initialize_output_file(time, depth, longitude, latitude, 
                    dAreas, vars_to_write, filename, scales.at(Iscale));
        }

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

        for (index = 0; index < num_pts; ++index) {
            coarse_KE.at(index) = 0.5 * constants::rho0 * (
                                    pow(coarse_u_x.at(index), 2.0)
                                  + pow(coarse_u_y.at(index), 2.0)
                                  + pow(coarse_u_z.at(index), 2.0)
                                );
        }

        if (constants::COMP_TRANSFERS) {
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
        }

        #if DEBUG >= 1
        fprintf(stdout, "  = Rank %d finished filtering =\n", wRank);
        fflush(stdout);
        #endif

        // Write to file
        if (not(constants::NO_FULL_OUTPUTS)) {
            write_field_to_output(coarse_u_x, "coarse_u_x", starts, counts, filename);
            write_field_to_output(coarse_u_y, "coarse_u_y", starts, counts, filename);
            write_field_to_output(coarse_u_z, "coarse_u_z",   starts, counts, filename);

            write_field_to_output(coarse_KE, "coarse_KE", starts, counts, filename);
        }

        if (constants::COMP_TRANSFERS) {
            // Compute the energy transfer through the filter scale
            compute_energy_transfer_through_scale(
                    energy_transfer, 
                    coarse_u_x,  coarse_u_y,  coarse_u_z,
                    coarse_uxux, coarse_uxuy, coarse_uxuz,
                    coarse_uyuy, coarse_uyuz, coarse_uzuz,
                    Ntime, Ndepth, Nlat, Nlon,
                    longitude, latitude, mask);

            if (not(constants::NO_FULL_OUTPUTS)) {
                write_field_to_output(energy_transfer, "energy_transfer", starts, counts, filename);
            }
        }


        if (not(constants::MINIMAL_OUTPUT)) {
            compute_div_vel(div, coarse_u_x, coarse_u_y, coarse_u_z, longitude, latitude,
                    Ntime, Ndepth, Nlat, Nlon, mask);
            if (not(constants::NO_FULL_OUTPUTS)) {
                write_field_to_output(div, "coarse_vel_div", starts, counts, filename);
            }

            compute_div_vel(div, full_u_x, full_u_y, full_u_z, longitude, latitude,
                    Ntime, Ndepth, Nlat, Nlon, mask);
            if (not(constants::NO_FULL_OUTPUTS)) {
                write_field_to_output(div, "full_vel_div", starts, counts, filename);
            }
        }

        if (constants::APPLY_POSTPROCESS) {
            MPI_Barrier(MPI_COMM_WORLD);

            if (wRank == 0) { fprintf(stdout, "Beginning post-process routines\n"); }
            fflush(stdout);

            Apply_Postprocess_Routines(
                    postprocess_fields, postprocess_names,
                    time, depth, latitude, longitude,
                    mask, dAreas,
                    myCounts, myStarts,
                    scales.at(Iscale));
        }


    }  // end for(scale) block

    // Cleanup some fftw stuff
    fftw_cleanup();
    fftw_destroy_plan(fft);
    fftw_destroy_plan(ifft);
    fftw_free(fft_inp);
    fftw_free(fft_out);

} // end filtering
