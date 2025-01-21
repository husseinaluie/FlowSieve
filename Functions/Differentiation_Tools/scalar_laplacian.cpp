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
        const std::vector<short int> & mask
    ) {

    const size_t Npts = scalar.size();

    // temporary work variables to store gradient component
    std::vector<double> grad_comp_1(Npts, 0), grad_comp_2(Npts, 0);

    // Compute gradient of scalar
    grad_scalar( grad_comp_1, grad_comp_2, scalar, source_data, mask );

    // Compute divergence of gradient
    div_vec( Lap_scalar, grad_comp_1, grad_comp_2, source_data, mask );

}
