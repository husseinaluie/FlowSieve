
#ifndef FFT_BASED_HPP
#define FFT_BASED_HPP 1

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <mpi.h>
#include <fftw3.h>
#include "constants.hpp"

//
//// FFTW-based scripts
//

void filtering_fftw(
        const std::vector<double> & u_r, 
        const std::vector<double> & u_lon, 
        const std::vector<double> & u_lat,
        const std::vector<double> & rho,
        const std::vector<double> & p,
        const std::vector<double> & scales, 
        const std::vector<double> & dAreas, 
        const std::vector<double> & time, 
        const std::vector<double> & depth,
        const std::vector<double> & longitude, 
        const std::vector<double> & latitude,
        const std::vector<double> & mask,
        const std::vector<int>    & myCounts,
        const std::vector<int>    & myStarts,
        const MPI_Comm comm = MPI_COMM_WORLD);

void fft_filter(
        std::vector<double> & coarse_field,
        const std::vector<double> & full_field,
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
        const double Lx);

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
        const double Lx);

#endif
