#include "../constants.hpp"
#include "../functions.hpp"
#include "../preprocess.hpp"
#include <algorithm>
#include <vector>
#include <omp.h>
#include <math.h>
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>

// Set the seed for the solver
void HelmholtzDataClass::Set_Seed(
        const dataset & data
        ) {

    x0.resize(2*Ncol);
    //x0[0] = 0.;
    //x0[Ncol] = 0.;

    size_t ii, index;
    const size_t Npts = data.mask.size();

    #pragma omp parallel default(none) \
    shared( pt_maps_to, num_mapped_before_col, data, x0 ) \
    private( ii, index ) \
    firstprivate( Npts, Ncol )
    {
        #pragma omp for collapse(1) schedule(static)
        for ( ii = 0; ii < Npts; ++ii ) {
            /*
            if ( (index > 0) and (pt_maps_to[ii] == ii) ) {
                x0[     index] = data.variables.at("Psi")[ii];
                x0[Ncol+index] = data.variables.at("Phi")[ii];
            }
            */
            if ( pt_maps_to[ii] != ii ) { continue; }
            index = ii - num_mapped_before_col[ii];
            x0[     index] = data.variables.at("Psi")[ii];
            x0[Ncol+index] = data.variables.at("Phi")[ii];
        }
    }

}
