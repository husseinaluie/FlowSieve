#include <vector>
#include <algorithm>
#include <math.h>
#include <cassert>
#include "../../differentiation_tools.hpp"
#include "../../constants.hpp"

#include <Eigen/Dense>

void non_uniform_diff_vector(
        std::vector<double> & diff_array,
        const std::vector<double> & grid,
        const int Iref,
        const int LB,
        const int UB,
        const int diff_ord) {

    // For now, only works for diff_ord == 2
    /*
    if ( diff_ord != 2 ) {
        throw std::runtime_error("non-uniform derivatives currently only implemented for 2nd order accuracy");
    }
    */

    const size_t Ndpts = UB-LB+1;
    const int Nconds = constants::DiffOrd + 0; //2;

    double scale_factor = 1.;
    
    // Calculate S1, Sx, Sxx, etc, and store into A (LHS)
    Eigen::Matrix<double, Nconds, Nconds> A;
    for ( int i = 0; i < Nconds; i++ ) {
        for ( int j = 0; j <= i; j++ ) {

            double Sij_value = 0.;

            for ( size_t In = LB; In <= UB; In++ ) {

                double dA = 1.;
                double x = grid.at(In) - grid.at(Iref);
                
                Sij_value += pow( x, i+j ) * dA;

            }

            A(i,j) = Sij_value;
            A(j,i) = Sij_value;
        }
    }

    // Precompute the LDLt solver
    Eigen::LDLT<Eigen::MatrixXd> ldlt;
    ldlt.compute(A);

    // Build the RHS and solve at each point
    diff_array.resize(Ndpts, 0.);
    for ( size_t In = LB; In <= UB; In++ ) {
        Eigen::Vector<double, Nconds> RHS, sol;

        double dA = 1.;
        double x = grid.at(In) - grid.at(Iref);

        for ( int i = 0; i < Nconds; i++ ) {
            RHS[i] = pow(x, i) * dA;
        }

        //sol = A.completeOrthogonalDecomposition().solve(RHS);
        //sol = A.ldlt().solve(RHS);
        sol = ldlt.solve(RHS);
        diff_array.at(In-LB) = sol[1];
    }
}
