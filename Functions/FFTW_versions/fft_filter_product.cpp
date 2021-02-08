#include <math.h>
#include <algorithm>
#include <vector>
#include <omp.h>
#include <mpi.h>
#include <fftw3.h>
#include "../../constants.hpp"
#include "../../functions.hpp"
#include "../../fft_based.hpp"

void fft_filter_product(
        std::vector<double> & coarse_field,
        const std::vector<double> & full_field1,
        const std::vector<double> & full_field2,
        const fftw_plan  fft,
        const fftw_plan ifft,
        double * fft_inp,
        fftw_complex * fft_out,
        const double scale,
        const int Ntime,
        const int Ndepth,
        const int Ny,
        const int Nx,
        const double Ly,
        const double Lx
        ) {

    int wRank, wSize;
    MPI_Comm_rank( MPI_COMM_WORLD, &wRank );
    MPI_Comm_size( MPI_COMM_WORLD, &wSize );

    size_t index;
    int Itime, Idepth, Iky, Ikx;
    double ky, kx, kk;
    const double dky = 2. * M_PI / Ly;
    const double dkx = 2. * M_PI / Lx;
    const int Nky = Ny;
    const int Nkx = Nx/2 + 1;

    #if DEBUG >= 1
    if (wRank == 0) { fprintf(stdout, "  Copying full field into fft input.\n"); }
    #endif
    // Copy full_field into fft_inp
    for (index = 0; index < full_field1.size(); index++) {
        fft_inp[index] = full_field1.at(index) * full_field2.at(index);
    }

    #if DEBUG >= 1
    if (wRank == 0) { fprintf(stdout, "  Forward transforming.\n"); }
    #endif
    // Forward transform
    fftw_execute(fft);

    #if DEBUG >= 1
    if (wRank == 0) { fprintf(stdout, "  Applying filter.\n"); }
    #endif
    // Apply filter
    for (Iky = 0; Iky < Nky; Iky++) {

        // y wavenumber
        if (Iky < Nky/2) { ky = dky * Iky; }
        else {             ky = dky * (Nky - Iky); }

        for (Ikx = 0; Ikx < Nkx; Ikx++) {

            // x wavenumber ( this one is halved, so no if/else )
            kx = dkx * Ikx;

            // Compute full
            kk = sqrt(kx*kx + ky*ky);

            // Remove anything higher frequency than filter scale
            if (kk >= 2. * M_PI / scale) {
                for (Itime = 0; Itime < Ntime; Itime++) {
                    for (Idepth = 0; Idepth < Ndepth; Idepth++) {
                        index = ( ( Itime * Ndepth + Idepth) * Nky + Iky ) * Nkx + Ikx;
                        fft_out[index][0] = 0.;
                        fft_out[index][1] = 0.;
                    }
                }
            }
        }
    }

    #if DEBUG >= 1
    if (wRank == 0) { fprintf(stdout, "  Back transforming.\n"); }
    #endif
    // Transform back
    fftw_execute(ifft);

    #if DEBUG >= 1
    if (wRank == 0) { fprintf(stdout, "  Copying fft output into coarse field.\n"); }
    #endif
    // Copy fft_inp into coarse_field
    for (index = 0; index < full_field1.size(); index++) {
        coarse_field.at(index) = fft_inp[index] / (Nx * Ny);
    }

    #if DEBUG >= 1
    if (wRank == 0) { fprintf(stdout, "  \n"); }
    #endif
}
