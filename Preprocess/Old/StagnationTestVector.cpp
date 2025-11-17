#include "../constants.hpp"
#include "../functions.hpp"
#include "../preprocess.hpp"
#include <algorithm>
#include <vector>
#include <omp.h>
#include <math.h>
#include <cassert>
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>

bool HelmholtzDataClass::StagnationTestVector( 
        ) const { 

    double norm_old_Psi = 0., norm_new_Psi = 0., norm_diff_Psi = 0.,
           norm_old_Phi = 0., norm_new_Phi = 0., norm_diff_Phi = 0.;
    const size_t Npts = Ncol;

    size_t Ipt;
    #pragma omp parallel default(none) private(Ipt) \
    firstprivate( Npts ) reduction( +:norm_old_Psi,norm_new_Psi,norm_diff_Psi,norm_old_Phi,norm_new_Phi,norm_diff_Phi )
    {
    #pragma omp for collapse(1) schedule(static)
        for ( Ipt = 0; Ipt < Npts; Ipt++ ) {
            norm_old_Psi  += pow( x0[Ipt],             2.);
            //norm_new_Psi  += pow( x0[Ipt] + soln[Ipt], 2.);
            //norm_diff_Psi += pow(           soln[Ipt], 2.);
            norm_new_Psi  += pow(           soln[Ipt], 2.);
            norm_diff_Psi += pow( x0[Ipt] - soln[Ipt], 2.);

            norm_old_Phi  += pow( x0[Ipt+Npts],                  2.);
            //norm_new_Phi  += pow( x0[Ipt+Npts] + soln[Ipt+Npts], 2.);
            //norm_diff_Phi += pow(                soln[Ipt+Npts], 2.);
            norm_new_Phi  += pow(                soln[Ipt+Npts], 2.);
            norm_diff_Phi += pow( x0[Ipt+Npts] - soln[Ipt+Npts], 2.);
        }
    }

    double stagnation_factor_Psi = sqrt(norm_diff_Psi) / ( sqrt(norm_new_Psi) + sqrt(norm_old_Psi) );
    double stagnation_factor_Phi = sqrt(norm_diff_Phi) / ( sqrt(norm_new_Phi) + sqrt(norm_old_Phi) );
    // By triangle inequality, stagnation_factors are in [0,1]
    // 0 indicates no change, 1 indicates largest possible change

    #if DEBUG >= 2
    fprintf( stdout, " Stagnation factors for (Psi,Phi) are (%.3e, %.3e)\n",
          stagnation_factor_Psi, stagnation_factor_Phi );
    #endif

    bool is_stagnated = ( stagnation_factor_Psi < stagnation_tolerance )
        and ( stagnation_factor_Phi < stagnation_tolerance );

    return is_stagnated;
}
