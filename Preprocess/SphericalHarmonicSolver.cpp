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
        std::vector<double> & Psi,
        std::vector<double> & Phi,
        const dataset & source_data,
        const std::vector<double> & vorticity,
        const std::vector<double> & divergence,
        const int max_L
        ) {

    //const int max_L = 16;
    const size_t Npts = source_data.mask.size();
    std::vector<long double> vort_coeffs_Real( max_L * (2 * max_L + 1), 0 ),
                             vort_coeffs_Imag( max_L * (2 * max_L + 1), 0 ),
                             div_coeffs_Real(  max_L * (2 * max_L + 1), 0 ),
                             div_coeffs_Imag(  max_L * (2 * max_L + 1), 0 );

    size_t Ipt;
    int L_degree;

    long double surface_area = 0.;
    #pragma omp parallel default(none) shared( source_data ) private(Ipt) \
    reduction( +:surface_area ) firstprivate( Npts )
    {
        #pragma omp for collapse(1) schedule(static)
        for ( Ipt = 0; Ipt < Npts; Ipt++ ) {
            surface_area += source_data.areas[Ipt];
        }
    }
    surface_area /= pow( (long double) constants::R_earth, 2. );
    fprintf( stdout, "1 - Numerical Surface Area / True Surface Area = %.3e\n", (double) (1 - surface_area / (4 * M_PI)) );


    #pragma omp parallel \
    default(none) \
    shared( source_data, vorticity, divergence ) \
    private( Ipt, L_degree ) \
    firstprivate( Npts, max_L ) \
    reduction( vec_long_double_plus:vort_coeffs_Real,vort_coeffs_Imag,div_coeffs_Real,div_coeffs_Imag )
    {
        // Loop over space to compute the inner produce of Lap_of_field with Y_L^M
        // Loop over the different spherical hamonic degrees and ordrs
        // We skip L = 0 since that's the mean
        #pragma omp for collapse(2) schedule(dynamic,16)
        for ( L_degree = 1; L_degree < max_L; L_degree++ ) {
            for ( Ipt = 0; Ipt < Npts; Ipt++ ) {

                // Divergence/Vorticity is zero on land, so avoid the computational cost
                if ( (vorticity[Ipt] == 0) and (divergence[Ipt] == 0) ) { continue; }

                long double lon, lat;
                if ( constants::GRID_TYPE == constants::GridType::MeshGrid ) {
                    int Itime, Idepth, Ilat, Ilon;
                    source_data.index1to4_local(Ipt, Itime, Idepth, Ilat, Ilon);
                    lon = source_data.longitude[Ilon];
                    lat = (M_PI/2) - source_data.latitude[Ilat];
                } else {
                    lon = source_data.longitude[Ipt];
                    lat = (M_PI/2) - source_data.latitude[Ipt]; //convert to polar angle
                }

                for ( int M_order = -L_degree; M_order <= L_degree; M_order++ ) {

                    long double Y_L_M = std::sph_legendre( L_degree, M_order, lat );

                    long double dA = source_data.areas[Ipt];

                    long double cos_M_lon = cos( M_order * lon );
                    long double sin_M_lon = sin( M_order * lon );

                    long double vort_L_M_Re = dA * Y_L_M * vorticity[Ipt] * cos_M_lon;
                    long double vort_L_M_Im = dA * Y_L_M * vorticity[Ipt] * sin_M_lon;

                    long double div_L_M_Re = dA * Y_L_M * divergence[Ipt] * cos_M_lon;
                    long double div_L_M_Im = dA * Y_L_M * divergence[Ipt] * sin_M_lon;

                    const double m_factor = ( M_order == 0 ) ? 0.5 : 1.;

                    vort_coeffs_Real[ L_degree * (2*max_L+1) + (M_order + max_L) ] += vort_L_M_Re * m_factor;
                    vort_coeffs_Imag[ L_degree * (2*max_L+1) + (M_order + max_L) ] += vort_L_M_Im * m_factor;

                    div_coeffs_Real[ L_degree * (2*max_L+1) + (M_order + max_L) ] += div_L_M_Re * m_factor;
                    div_coeffs_Imag[ L_degree * (2*max_L+1) + (M_order + max_L) ] += div_L_M_Im * m_factor;
                }
            }
        }
    }


    // Now that we have the coefficients, use them to build the new field
    std::vector<long double> long_Psi(Npts, 0.), long_Phi(Npts, 0.);

    #pragma omp parallel \
    default(none) \
    shared( source_data, long_Psi, long_Phi, vort_coeffs_Real, vort_coeffs_Imag, div_coeffs_Real, div_coeffs_Imag ) \
    private( Ipt, L_degree ) \
    firstprivate( Npts, surface_area, max_L )
    {
        // Loop over the different spherical hamonic degrees and ordrs
        #pragma omp for collapse(1) schedule(static)
        for ( Ipt = 0; Ipt < Npts; Ipt++ ) {

            long double lon, lat;
            if ( constants::GRID_TYPE == constants::GridType::MeshGrid ) {
                int Itime, Idepth, Ilat, Ilon;
                source_data.index1to4_local(Ipt, Itime, Idepth, Ilat, Ilon);
                lon = source_data.longitude[Ilon];
                lat = (M_PI/2) - source_data.latitude[Ilat];
            } else {
                lon = source_data.longitude[Ipt];
                lat = (M_PI/2) - source_data.latitude[Ipt]; //convert to polar angle
            }

            for ( L_degree = 1; L_degree < max_L; L_degree++ ) {
                for ( int M_order = -L_degree; M_order <= L_degree; M_order++ ) {
                    long double Y_L_M = std::sph_legendre( L_degree, M_order, lat );

                    long double vort_L_M_Re = vort_coeffs_Real[ L_degree * (2*max_L+1) + (M_order + max_L) ];
                    long double vort_L_M_Im = vort_coeffs_Imag[ L_degree * (2*max_L+1) + (M_order + max_L) ];

                    long double div_L_M_Re = div_coeffs_Real[ L_degree * (2*max_L+1) + (M_order + max_L) ];
                    long double div_L_M_Im = div_coeffs_Imag[ L_degree * (2*max_L+1) + (M_order + max_L) ];

                    vort_L_M_Re = - vort_L_M_Re / ( L_degree * (L_degree+1) );
                    vort_L_M_Im = - vort_L_M_Im / ( L_degree * (L_degree+1) );

                    div_L_M_Re = - div_L_M_Re / ( L_degree * (L_degree+1) );
                    div_L_M_Im = - div_L_M_Im / ( L_degree * (L_degree+1) );

                    // Our input field is real, so only keep the real parts
                    long double cos_M_lon = cos( M_order * lon );
                    long double sin_M_lon = sin( M_order * lon );

                    long double Psi_from_real = vort_L_M_Re * Y_L_M * cos_M_lon;
                    long double Psi_from_imag = vort_L_M_Im * Y_L_M * sin_M_lon;

                    #ifdef NDEBUG
                    long_Psi[Ipt] += (Psi_from_real + Psi_from_imag) * sqrt( surface_area );
                    #else
                    long_Psi.at(Ipt) += (Psi_from_real + Psi_from_imag) * sqrt( surface_area );
                    #endif

                    long double Phi_from_real = div_L_M_Re * Y_L_M * cos_M_lon;
                    long double Phi_from_imag = div_L_M_Im * Y_L_M * sin_M_lon;

                    #ifdef NDEBUG
                    long_Phi[Ipt] += (Phi_from_real + Phi_from_imag) * sqrt( surface_area );
                    #else
                    long_Phi.at(Ipt) += (Phi_from_real + Phi_from_imag) * sqrt( surface_area );
                    #endif
                }
            }
        }
    }

    // Finally, adjust so that Psi[0] = 0, Phi[0] = 0
    #pragma omp parallel \
    default(none) shared( long_Psi, long_Phi, Psi, Phi ) \
    private( Ipt ) firstprivate( Npts )
    {
        #pragma omp for collapse(1) schedule(static)
        for ( Ipt = 1; Ipt < Npts; Ipt++ ) {
            Psi[Ipt] = long_Psi[Ipt] - long_Psi[0];
            Phi[Ipt] = long_Phi[Ipt] - long_Phi[0];

        }
    }
    Psi[0] = 0.;
    Phi[0] = 0.;

}
