#include "../constants.hpp"
#include "../functions.hpp"
#include "../preprocess.hpp"
#include <algorithm>
#include <vector>
#include <omp.h>
//#include <math.h>
#include <cmath>

// Pragma-magic to allow reduction over vector operator
//      thanks to: https://stackoverflow.com/questions/43168661/openmp-and-reduction-on-stdvector
#pragma omp declare reduction(vec_double_plus : std::vector<double> : \
                              std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<double>())) \
                    initializer(omp_priv = decltype(omp_orig)(omp_orig.size()))

#pragma omp declare reduction(vec_long_double_plus : std::vector<long double> : \
                              std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<long double>())) \
                    initializer(omp_priv = decltype(omp_orig)(omp_orig.size()))


void SphericalHarmonicSolver(
        std::vector<double> & reconstructed_field,
        const dataset & source_data,
        const std::vector<double> & Laplacian_of_field
        ) {

    const int max_L = 128;
    const size_t Npts = source_data.mask.size();
    std::vector<long double> sph_coeffs_Real( max_L * 2 * max_L, 0 );
    std::vector<long double> sph_coeffs_Imag( max_L * 2 * max_L, 0 );

    size_t Ipt;
    int L_degree, M_order;

    #pragma omp parallel \
    default(none) \
    shared( source_data, Laplacian_of_field ) \
    private( Ipt, L_degree, M_order ) \
    firstprivate( Npts ) \
    reduction( vec_long_double_plus:sph_coeffs_Real,sph_coeffs_Imag )
    {
        // Loop over space to compute the inner produce of Lap_of_field with Y_L^M
        // Loop over the different spherical hamonic degrees and ordrs
        // We skip L = 0 since that's the mean
        #pragma omp for collapse(2) schedule(guided)
        for ( L_degree = 1; L_degree < max_L; L_degree++ ) {
            for ( Ipt = 0; Ipt < Npts; Ipt++ ) {
                // Divergence/Vorticity is zero of land, so avoid the computational cost
                if ( not(source_data.mask[Ipt]) ) { continue; }
                for ( int M_order = -L_degree; M_order <= L_degree; M_order++ ) {
                    long double lon = source_data.longitude[Ipt];
                    long double lat = (M_PI/2) - source_data.latitude[Ipt]; //convert to polar angle

                    long double Y_L_M = std::sph_legendre( L_degree, M_order, lat );

                    long double dA = source_data.areas[Ipt];

                    long double alpha_L_M_Re = dA * Y_L_M * Laplacian_of_field[Ipt] * cos( M_order * lon );
                    long double alpha_L_M_Im = dA * Y_L_M * Laplacian_of_field[Ipt] * sin( M_order * lon );

                    sph_coeffs_Real[ L_degree * (2*max_L) + (M_order + max_L) ] += alpha_L_M_Re;
                    sph_coeffs_Imag[ L_degree * (2*max_L) + (M_order + max_L) ] += alpha_L_M_Im;
                }
            }
        }
    }


    // Now that we have the coefficients, use them to build the new field
    reconstructed_field.resize(Npts, 0);

    #pragma omp parallel \
    default(none) \
    shared( source_data, reconstructed_field, sph_coeffs_Real, sph_coeffs_Imag ) \
    private( Ipt, L_degree, M_order ) \
    firstprivate( Npts )
    {
        // Loop over the different spherical hamonic degrees and ordrs
        #pragma omp for collapse(1) schedule(static)
        for ( Ipt = 0; Ipt < Npts; Ipt++ ) {
            long double lon = source_data.longitude[Ipt];
            long double lat = (M_PI/2) - source_data.latitude[Ipt]; //convert to polar angle

            for ( L_degree = 1; L_degree < max_L; L_degree++ ) {
                for ( M_order = -L_degree; M_order <= L_degree; M_order++ ) {
                    long double Y_L_M = std::sph_legendre( L_degree, M_order, lat );

                    long double alpha_L_M_Re = sph_coeffs_Real[ L_degree * (2*max_L) + (M_order + max_L) ];
                    long double alpha_L_M_Im = sph_coeffs_Imag[ L_degree * (2*max_L) + (M_order + max_L) ];

                    alpha_L_M_Re = - alpha_L_M_Re / ( L_degree * (L_degree+1) );
                    alpha_L_M_Im = - alpha_L_M_Im / ( L_degree * (L_degree+1) );

                    // Our input field is real, so only keep the real parts
                    long double real_from_real = alpha_L_M_Re * Y_L_M * cos( M_order * lon );
                    long double real_from_imag = alpha_L_M_Im * Y_L_M * sin( M_order * lon );

                    //reconstructed_field[Ipt] += real_from_real - real_from_imag;
                    reconstructed_field[Ipt] += (double) (real_from_real + real_from_imag);
                }
            }
        }
    }


    // Finally, adjust so that Psi[0] = 0, Phi[0] = 0
    #pragma omp parallel \
    default(none) \
    shared( reconstructed_field ) \
    private( Ipt ) \
    firstprivate( Npts )
    {
        #pragma omp for collapse(1) schedule(static)
        for ( Ipt = 1; Ipt < Npts; Ipt++ ) {
            reconstructed_field[Ipt] = reconstructed_field[Ipt] - reconstructed_field[0];
        }
    }
    reconstructed_field[0] = 0.;

}
