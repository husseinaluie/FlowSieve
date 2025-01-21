#include "../constants.hpp"
#include "../functions.hpp"
#include "../preprocess.hpp"
#include <algorithm>
#include <vector>
#include <omp.h>
#include <math.h>
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>

// Build the RHS of the problem
void HelmholtzDataClass::Build_RHS(
        const dataset & data
        ) {

    size_t Ipt, index;
    double w;
    size_t Npts = data.mask.size();
    RHS.resize( (use_vort_div and use_vel) ? 4*Nrow : 2*Nrow );
    #pragma omp parallel default(none) \
    private( Ipt, index ) \
    firstprivate( Nrow )
    {
        #pragma omp for collapse(1) schedule(static)
        for ( Ipt = 0; Ipt < ( (use_vort_div and use_vel) ? 4*Nrow : 2*Nrow); ++Ipt) {
            RHS[Ipt] = 0.;
        }
    }

    #pragma omp parallel default(none) \
    shared( data ) \
    private( Ipt, index, w ) \
    firstprivate( Npts )
    {
        #pragma omp for collapse(1) schedule(static)
        for ( Ipt = 0; Ipt < Npts; ++Ipt) {
            if ( pt_maps_to[Ipt] != Ipt ) { continue; } // skip mapped points
            if ( Ipt == 0 ) { continue; }
            if ( not(data.mask[Ipt]) ) { continue; } // if it's land, skip [i.e. set RHS to zero]
            w = weight_err ? data.areas.at(Ipt) : 1.;
            index = Ipt - num_mapped_before_row[Ipt];
            if ( use_vel ) {
                RHS[         index] = data.variables.at("u_lon").at(Ipt) * w;
                RHS[  Nrow + index] = data.variables.at("u_lat").at(Ipt) * w;
            }
            if ( use_vort_div ) {
                const size_t base = use_vel ? 2*Nrow : 0;
                RHS[base +        index] = data.variables.at("vort" ).at(Ipt) * w * Tikhov;
                RHS[base + Nrow + index] = data.variables.at("div"  ).at(Ipt) * w * Tikhov;
            }
        }
    }

    #if DEBUG >= 2
    // Get norms of the RHS segments
    double uo_norm = 0, vo_norm = 0, vort_norm = 0, div_norm = 0;
    #pragma omp parallel default(none) \
    private( Ipt ) \
    reduction( +:uo_norm,vo_norm,vort_norm,div_norm )
    {
        #pragma omp for collapse(1) schedule(static)
        for ( Ipt = 0; Ipt < Nrow; ++Ipt) {
            if ( use_vel ) {
                uo_norm   += pow(RHS[         Ipt], 2.);
                vo_norm   += pow(RHS[  Nrow + Ipt], 2.);
            }
            if ( use_vort_div ) {
                const size_t base = use_vel ? 2*Nrow : 0;
                vort_norm += pow(RHS[base +        Ipt], 2.);
                div_norm  += pow(RHS[base + Nrow + Ipt], 2.);
            }
        }
    }
    if ( use_vel ) {
        uo_norm   = sqrt( uo_norm   / Nrow );
        vo_norm   = sqrt( vo_norm   / Nrow );
    }
    if ( use_vort_div ) {
        vort_norm = sqrt( vort_norm / Nrow );
        div_norm  = sqrt( div_norm  / Nrow );
    }

    if ( use_vort_div and use_vel ) {
        fprintf( stdout, "RHS segment norms are %.2e, %.2e, %.2e, and %.2e\n",
                uo_norm, vo_norm, vort_norm, div_norm );
    } else if ( use_vort_div ) {
        fprintf( stdout, "RHS segment norms are %.2e, and %.2e\n",
                vort_norm, div_norm );
    } else {
        fprintf( stdout, "RHS segment norms are %.2e and %.2e\n",
                uo_norm, vo_norm );
    }
    #endif

}
