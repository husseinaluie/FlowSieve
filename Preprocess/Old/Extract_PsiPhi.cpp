#include "../constants.hpp"
#include "../functions.hpp"
#include "../preprocess.hpp"
#include <algorithm>
#include <vector>
#include <omp.h>
#include <math.h>


// Extract Psi/Phi from the solver grid onto the physical grid
void HelmholtzDataClass::Extract_PsiPhi(
        dataset & data
        ) {


    // Zero out first
    std::fill( data.variables.at("Psi").begin(),
               data.variables.at("Psi").end(),
               0. );
    std::fill( data.variables.at("Phi").begin(),
               data.variables.at("Phi").end(),
               0. );

    // We'll extract from x0, since it stores the accumulated solution

    size_t ii, index;
    const size_t Npts = data.mask.size();
    #pragma omp parallel default(none) \
    shared( data, pt_maps_to, num_mapped_before_col, soln, x0 ) \
    private( ii, index ) \
    firstprivate( Npts ) 
    {
        #pragma omp for collapse(1) schedule(static)
        for ( ii = 0; ii < Npts; ++ii ) {
            if ( pt_maps_to[ii] == 0 ) { continue; }
            index = pt_maps_to[ii] - num_mapped_before_col[pt_maps_to[ii]];
            data.variables.at("Psi")[ii] = x0[       index];
            data.variables.at("Phi")[ii] = x0[Ncol + index];
        }
    }

}
