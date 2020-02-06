#include <math.h>
#include <algorithm>
#include <vector>
#include <omp.h>
#include "../../functions.hpp"
#include "../../functions_sw.hpp"

void Compute_misc2_SW(
        std::vector<double> & out_array, 
        const std::vector<double> & u_tilde, 
        const std::vector<double> & v_tilde, 
        const std::vector<double> & tau_hpx, 
        const std::vector<double> & tau_hpy, 
        const std::vector<double> & time,
        const std::vector<double> & depth,
        const std::vector<double> & latitude,
        const std::vector<double> & longitude,
        const std::vector<double> & mask,
        const int & Ntime,
        const int & Ndepth,
        const int & Nlat,
        const int & Nlon
        ){

    const int num_pts = u_tilde.size();
    int index;

    //
    #pragma omp parallel \
    default(none) \
    shared(u_tilde, v_tilde, tau_hpx, tau_hpy, out_array) \
    private(index)
    {
        #pragma omp for collapse(1) schedule(static)
        for (index = 0; index < num_pts; ++index) {

            out_array.at(index) = 
                  u_tilde.at(index) * tau_hpx.at(index) 
                + v_tilde.at(index) * tau_hpy.at(index);
        }
    }
}

