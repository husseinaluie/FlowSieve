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
void HelmholtzDataClass::Extract_ProjectedVars(
        dataset & data       
        ) {

    // Get the product A * x (which should be ~b)
    // note that x0 stores the accumulated solution
    Eigen::VectorXd Ax = LHS * x0;

    const size_t Npts = data.mask.size();
    size_t Ipt, Ipt_mapped;
    double w;

    #pragma omp parallel default(none) \
    shared( data, Ax ) \
    private( Ipt, Ipt_mapped, w ) \
    firstprivate( Npts )
    {
        #pragma omp for collapse(1) schedule(static)
        for ( Ipt = 0; Ipt < Npts; Ipt++ ) {

            w = weight_err ? data.areas.at(Ipt) : 1.;
            if ( all_land_neighbours[Ipt] ) {
                // If it's landlocked, then we mapped it away, so just zero
                if ( use_vel ) {
                    data.variables.at("projected_u_lon")[Ipt] = 0.;
                    data.variables.at("projected_u_lat")[Ipt] = 0.;
                }
                if ( use_vort_div ) {
                    data.variables.at("projected_vort" )[Ipt] = 0.;
                    data.variables.at("projected_div"  )[Ipt] = 0.;
                }
                continue;
            }
            //Ipt_mapped = pt_maps_to[Ipt] - num_mapped_before_row[pt_maps_to[Ipt]];
            // This doesn't get the right value for coastal

            Ipt_mapped = Ipt - num_mapped_before_row[Ipt];

            if ( use_vel ) {
                data.variables.at("projected_u_lon")[Ipt] = Ax[0*Nrow + Ipt_mapped] / w;
                data.variables.at("projected_u_lat")[Ipt] = Ax[1*Nrow + Ipt_mapped] / w;
            }
            if ( use_vort_div ) {
                const size_t base = use_vel ? 2*Nrow : 0;
                data.variables.at("projected_vort" )[Ipt] = Ax[base +        Ipt_mapped] / (w*Tikhov);
                data.variables.at("projected_div"  )[Ipt] = Ax[base + Nrow + Ipt_mapped] / (w*Tikhov);
            }
        }
    }
}
