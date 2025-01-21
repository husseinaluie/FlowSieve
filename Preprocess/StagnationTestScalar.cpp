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

bool HelmholtzDataClass::StagnationTestScalar( 
        ) const { 

    double norm_old = 0., norm_new=0., norm_diff=0.;
    const size_t Npts = x0.size();

    size_t Ipt;
    #pragma omp parallel default(none) private(Ipt) \
    firstprivate( Npts ) reduction( +:norm_old,norm_new,norm_diff )
    {
    #pragma omp for collapse(1) schedule(static)
        for ( Ipt = 0; Ipt < Npts; Ipt++ ) {
            norm_old  += pow( x0[Ipt],             2.);
            norm_new  += pow( x0[Ipt] + soln[Ipt], 2.);
            norm_diff += pow(           soln[Ipt], 2.);
        }
    }

    double stagnation_factor = sqrt(norm_diff) / ( sqrt(norm_new) + sqrt(norm_old) );
    // By triangle inequality, stagnation_factor is in [0,1]
    // 0 indicates no change, 1 indicates largest possible change

    #if DEBUG >= 2
    fprintf( stdout, " Stagnation factors is %.3e\n",
          stagnation_factor );
    #endif

    return stagnation_factor < stagnation_tolerance;
}
