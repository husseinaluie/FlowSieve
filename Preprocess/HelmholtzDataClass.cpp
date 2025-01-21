#include "../constants.hpp"
#include "../functions.hpp"
#include "../preprocess.hpp"
#include <algorithm>
#include <vector>
#include <omp.h>
#include <math.h>

// Class Constructor
HelmholtzDataClass::HelmholtzDataClass() {
};

// Set the solver
void HelmholtzDataClass::InitializeSolver() {
    solver.compute( LHS );
    if ( solver.info() == Eigen::NumericalIssue ) {
        fprintf( stderr, "The provided data did not satisfy the prerequisites..\n" );
        throw std::runtime_error("Failed to initialize the solver.");
    } else if ( solver.info() == Eigen::NoConvergence ) {
        fprintf( stderr, "Iterative procedure did not converge.\n" );
        throw std::runtime_error("Failed to initialize the solver.");
    } else if ( solver.info() == Eigen::InvalidInput ) {
        fprintf( stderr, "The inputs are invalid, or the algorithm has been improperly called.\n" );
        throw std::runtime_error("Failed to initialize the solver.");
    } else if ( solver.info() != Eigen::Success ) {
        fprintf( stderr, "Eigen decomposition failed in an unknown way.\n" );
        throw std::runtime_error("Failed to initialize the solver.");
    }
};

// Apply the solver with guess x0
void HelmholtzDataClass::Solve() {
    //Eigen::VectorXd RHS_tmp = RHS - LHS * x0;
    //soln = solver.solve( RHS_tmp );
    soln = solver.solveWithGuess( RHS, x0 );
};

// In the case of small grids, do a direct solve
void HelmholtzDataClass::DirectSolve() {
    direct_solver.compute(LHS);
    if ( solver.info() == Eigen::NumericalIssue ) {
        fprintf( stderr, "The provided data did not satisfy the prerequisites..\n" );
        throw std::runtime_error("Failed to initialize the solver.");
    } else if ( solver.info() == Eigen::NoConvergence ) {
        fprintf( stderr, "Iterative procedure did not converge.\n" );
        throw std::runtime_error("Failed to initialize the solver.");
    } else if ( solver.info() == Eigen::InvalidInput ) {
        fprintf( stderr, "The inputs are invalid, or the algorithm has been improperly called.\n" );
        throw std::runtime_error("Failed to initialize the solver.");
    } else if ( solver.info() != Eigen::Success ) {
        fprintf( stderr, "Eigen decomposition failed in an unknown way.\n" );
        throw std::runtime_error("Failed to initialize the solver.");
    }
    soln = direct_solver.solve( RHS );
};


bool HelmholtzDataClass::IsConverged() const {
    
    const double vel_2_REL =  use_vel ? vel_2_errors.back()  / vel_2_norms.back() : 1.;
    if ( not(use_vort_div) ) {
        return vel_2_REL < tolerance;
    }

    const double vort_2_REL = use_vort_div ? vort_2_errors.back() / vort_2_norms.back() : 1.;
    const double div_2_REL  = use_vort_div ? div_2_errors.back()  / div_2_norms.back() : 1.;
    if ( not(use_vel) ) {
        return ( vort_2_REL < tolerance ) and ( div_2_REL < tolerance );
    }

    const bool all_below_tolerance = 
        ( vel_2_REL < tolerance ) and ( vort_2_REL < tolerance ) and ( div_2_REL < tolerance );
                            
    return all_below_tolerance;
    // Add conditions for inf norm too
}

// Clear
void HelmholtzDataClass::clear() {
    all_land_neighbours.clear();
    pt_maps_to.clear();
    num_mapped_before_col.clear();
    num_mapped_before_row.clear();
    num_land_before.clear();

    island_reps.clear();
    coastal_boundaries.clear();

    num_coastal = 0;
    num_all_land = 0;
    Ncol = 0;
    Nrow = 0;
    Npts_mapped = 0;

    // Eigen
    LHS.resize( 0,0 );
    RHS.resize( 0 );
    x0.resize( 0 );
    soln.resize( 0 );

    // Convergence tracking
    vel_2_errors.clear();
    vort_2_errors.clear();
    div_2_errors.clear();
    vel_inf_errors.clear();
    vort_inf_errors.clear();
    div_inf_errors.clear();
    vel_2_norms.clear();
    vort_2_norms.clear();
    div_2_norms.clear();
    vel_inf_norms.clear();
    vort_inf_norms.clear();
    div_inf_norms.clear();
};
