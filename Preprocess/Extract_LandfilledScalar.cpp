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

// Compute projection errors
void HelmholtzDataClass::Extract_LandfilledScalar(
        dataset & data       
        ) {

    const size_t Npts = data.mask.size();
    size_t Ipt, Ipt_mapped;
    double w;

    fprintf( stdout, "Npts = %'zu, len(x0) = %'zu, max(Ipt_mapped) = %'zu\n",
          Npts, x0.size(), num_land_before.back() );
    /*
    fprintf( stdout, "len(land-filled-scalar) = %'zu, len(source-scalar) = %'zu\n",
          data.variables.at("land_filled_scalar").size(),
          data.variables.at("scalar").size()
          );
    */

    #pragma omp parallel default(none) \
    shared( data ) \
    private( Ipt, Ipt_mapped, w ) \
    firstprivate( Npts )
    {
        #pragma omp for collapse(1) schedule(static)
        for ( Ipt = 0; Ipt < Npts; Ipt++ ) {

            if ( data.mask[Ipt] ) {
                // If it's water, keep the original value
                data.variables.at("land_filled_scalar")[Ipt] = 
                    data.variables.at("scalar")[Ipt];
            } else {
                //w = weight_err ? data.areas.at(Ipt) : 1.;
                w = 1.;
                Ipt_mapped = num_land_before[Ipt];
                if (Ipt_mapped >= x0.size()) {
                    data.variables.at("land_filled_scalar")[Ipt] = -1e5;
                } else {
                    data.variables.at("land_filled_scalar")[Ipt] = x0[Ipt_mapped] / w;
                }
            }
        }
    }
}
