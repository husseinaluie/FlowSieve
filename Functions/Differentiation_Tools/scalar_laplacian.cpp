#include "../../constants.hpp"
#include "../../functions.hpp"
#include "../../differentiation_tools.hpp"
#include <algorithm>
#include <vector>
#include <omp.h>
#include <math.h>

void scalar_laplacian(  
        std::vector<double> & Lap_scalar,
        const std::vector<double> & scalar,
        const dataset & source_data,
        const std::vector<short int> & mask,
        std::vector<double>* work_grad_comp_1,
        std::vector<double>* work_grad_comp_2
    ) {

    const size_t Npts = scalar.size();

    // back-up storage in case not user provided [not allocated unless needed]
    std::vector<double> local_grad_comp_1;
    std::vector<double> local_grad_comp_2;

    // temporary work variables to store gradient component
    //  if they were provided by the user, us the pre-allocated version
    std::vector<double>& grad_comp_1 =
        work_grad_comp_1 ? *work_grad_comp_1 : local_grad_comp_1;

    std::vector<double>& grad_comp_2 =
        work_grad_comp_2 ? *work_grad_comp_2 : local_grad_comp_2;

    if (!work_grad_comp_1) { grad_comp_1.resize(Npts, 0.); }
    if (!work_grad_comp_2) { grad_comp_2.resize(Npts, 0.); }

    if ( grad_comp_1.size() != Npts ) {
        throw std::runtime_error("grad_comp_1 storage array has the wrong size.");
    }
    if ( grad_comp_2.size() != Npts ) {
        throw std::runtime_error("grad_comp_2 storage array has the wrong size.");
    }

    // Compute gradient of scalar
    grad_scalar( grad_comp_1, grad_comp_2, scalar, source_data, mask );

    // Compute divergence of gradient
    div_vec( Lap_scalar, grad_comp_1, grad_comp_2, source_data, mask );

}
