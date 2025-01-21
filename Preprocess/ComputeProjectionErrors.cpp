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
void HelmholtzDataClass::ComputeProjectionErrors(
        const dataset & data       
        ) {

    // Get the product A * x (which should be ~b)
    Eigen::VectorXd Ax = LHS * x0;

    size_t Ipt, Ipt_mapped;

    #if DEBUG >= 2
    // Get L2 error w.r.t. RHS for comparison
    double RHS_u_err = 0, RHS_v_err = 0, RHS_vort_err = 0, RHS_div_err = 0,
           RHS_u_ref = 0, RHS_v_ref = 0, RHS_vort_ref = 0, RHS_div_ref = 0;

    #pragma omp parallel default(none) \
    shared( data, Ax ) \
    private( Ipt ) \
    firstprivate( Ncol ) \
    reduction( +:RHS_u_err,RHS_v_err,RHS_vort_err,RHS_div_err,RHS_u_ref,RHS_v_ref,RHS_vort_ref,RHS_div_ref )
    {
        #pragma omp for collapse(1) schedule(static)
        for ( Ipt = 0; Ipt < Nrow; Ipt++ ) {

            // Velocity Errors
            if ( use_vel ) {
                RHS_u_err    += pow(  Ax[0*Nrow + Ipt] - RHS[0*Nrow + Ipt], 2);
                RHS_v_err    += pow(  Ax[1*Nrow + Ipt] - RHS[1*Nrow + Ipt], 2);

                RHS_u_ref    += pow(  RHS[0*Nrow + Ipt], 2);
                RHS_v_ref    += pow(  RHS[1*Nrow + Ipt], 2);
            }

            if ( use_vort_div ) {
                const size_t base = use_vel ? 2*Nrow : 0;
                RHS_vort_err += pow(  Ax[base +        Ipt] - RHS[base +        Ipt], 2);
                RHS_div_err  += pow(  Ax[base + Nrow + Ipt] - RHS[base + Nrow + Ipt], 2);
                RHS_vort_ref += pow(  RHS[base +        Ipt], 2);
                RHS_div_ref  += pow(  RHS[base + Nrow + Ipt], 2);
            }

        }
    }
    fprintf(stdout, "(u,v,vort,div) -> (%.2e, %.2e, %.2e, %.2e) / (%.2e, %.2e, %.2e, %.2e)\n",
            RHS_u_err, RHS_v_err, RHS_vort_err, RHS_div_err, 
            RHS_u_ref, RHS_v_ref, RHS_vort_ref, RHS_div_ref );
    #endif

    // Now loop over space and get the L2 and Linf errors for
    // vel, vort, div
    double vel_2_err  = 0, vel_inf_err  = 0,
           vort_2_err = 0, vort_inf_err = 0,
           div_2_err  = 0, div_inf_err  = 0;
    double vel_2_norm  = 0, vel_inf_norm  = 0,
           vort_2_norm = 0, vort_inf_norm = 0,
           div_2_norm  = 0, div_inf_norm  = 0;


    const size_t Npts = data.mask.size();
    double area, total_area = 0, incr, w;

    #pragma omp parallel default(none) \
    shared( data, pt_maps_to, num_mapped_before_row, Ax ) \
    private( Ipt, incr ) \
    firstprivate( Npts, stdout ) \
    reduction( +:total_area,vel_2_err,vort_2_err,div_2_err,vel_2_norm,vort_2_norm,div_2_norm ) \
    reduction( max:vel_inf_err,vort_inf_err,div_inf_err,vel_inf_norm,vort_inf_norm,div_inf_norm )
    {
        #pragma omp for collapse(1) schedule(static)
        for ( Ipt = 0; Ipt < Npts; Ipt++ ) {

            if (pt_maps_to[Ipt] == 0) { continue; }

            double area = data.areas[Ipt];
            double w = weight_err ? data.areas.at(Ipt) : 1.;
            size_t Ipt_mapped = Ipt - num_mapped_before_row[Ipt];

            total_area += area;

            if ( use_vel ) {
                // Velocity Errors
                incr = sqrt(   pow(  Ax[0*Nrow + Ipt_mapped] / w - data.variables.at("u_lon")[Ipt], 2)
                        + pow(  Ax[1*Nrow + Ipt_mapped] / w - data.variables.at("u_lat")[Ipt], 2)
                        );
                vel_inf_err = std::fmax( incr, vel_inf_err );
                vel_2_err += area * pow( incr, 2 ); 
            }

            if ( use_vort_div ) {
                const size_t base = use_vel ? 2*Nrow : 0;

                // Vorticity errors
                incr = std::fabs( Ax[base + Ipt_mapped] / (w*Tikhov) - data.variables.at("vort")[Ipt] );
                vort_inf_err = std::fmax( incr, vort_inf_err );
                vort_2_err += area * pow( incr, 2 );

                // Divergence errors
                incr = std::fabs( Ax[base + Nrow + Ipt_mapped] / (w*Tikhov) - data.variables.at("div")[Ipt] );
                div_inf_err = std::fmax( incr, div_inf_err );
                div_2_err += area * pow( incr, 2 );
            }

            // Reference norms
            if ( use_vel ) {
                incr = sqrt(   pow(  data.variables.at("u_lon")[Ipt], 2)
                        + pow(  data.variables.at("u_lat")[Ipt], 2)
                        );
                vel_inf_norm = std::fmax( incr, vel_inf_norm );
                vel_2_norm += area * pow( incr, 2 ); 
            }

            if ( use_vort_div ) {
                incr = std::fabs( data.variables.at("vort")[Ipt] );
                vort_inf_norm = std::fmax( incr, vort_inf_norm );
                vort_2_norm += area * pow( incr, 2 );

                incr = std::fabs( data.variables.at("div")[Ipt] );
                div_inf_norm = std::fmax( incr, div_inf_norm );
                div_2_norm += area * pow( incr, 2 );
            }
        }
    }

    vel_2_err  = sqrt( vel_2_err  / total_area );
    vort_2_err = sqrt( vort_2_err / total_area );
    div_2_err  = sqrt( div_2_err  / total_area );

    // Save the convergence records
    vel_2_errors.push_back(  vel_2_err  );
    vort_2_errors.push_back( vort_2_err );
    div_2_errors.push_back(  div_2_err  );

    vel_inf_errors.push_back(  vel_inf_err  );
    vort_inf_errors.push_back( vort_inf_err );
    div_inf_errors.push_back(  div_inf_err  );

    // Save the convergence reference values
    vel_2_norm  = sqrt( vel_2_norm  / total_area );
    vort_2_norm = sqrt( vort_2_norm / total_area );
    div_2_norm  = sqrt( div_2_norm  / total_area );

    vel_2_norms.push_back(  vel_2_norm  );
    vort_2_norms.push_back( vort_2_norm );
    div_2_norms.push_back(  div_2_norm  );

    vel_inf_norms.push_back(  vel_inf_norm  );
    vort_inf_norms.push_back( vort_inf_norm );
    div_inf_norms.push_back(  div_inf_norm  );

}
