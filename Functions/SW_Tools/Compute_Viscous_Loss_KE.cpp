#include <math.h>
#include <algorithm>
#include <vector>
#include <omp.h>
#include "../../functions.hpp"
#include "../../functions_sw.hpp"

void Compute_Viscous_Loss_KE(
        std::vector<double> & out_array, 
        const std::vector<double> & h_bar, 
        const std::vector<double> & u_tilde, 
        const std::vector<double> & v_tilde, 
        const std::vector<double> & lap_u_tilde, 
        const std::vector<double> & lap_v_tilde, 
        const double nu,
        const std::vector<double> & rho, 
        const std::vector<double> & time,
        const std::vector<double> & depth,
        const std::vector<double> & latitude,
        const std::vector<double> & longitude,
        const int & Ntime,
        const int & Ndepth,
        const int & Nlat,
        const int & Nlon
        ){

    const int num_pts = Ntime * Ndepth * Nlat * Nlon;

    int Itime, Idepth, Ilat, Ilon, index;

    //
    #pragma omp parallel \
    default(none) \
    shared(h_bar, u_tilde, v_tilde, lap_u_tilde, lap_v_tilde,\
            rho, out_array, Ntime, Ndepth, Nlat, Nlon) \
    private(index, Itime, Idepth, Ilat, Ilon)
    {

        #pragma omp for collapse(1) schedule(static)
        for (index = 0; index < num_pts; index++) {

            Index1to4(index, Itime, Idepth, Ilat, Ilon,
                             Ntime, Ndepth, Nlat, Nlon);

            out_array.at(index) = nu * rho.at(Idepth) * h_bar.at(index)
                                    * (    u_tilde.at(index) * lap_u_tilde.at(index)
                                        +  v_tilde.at(index) * lap_v_tilde.at(index) 
                                      );

        }
    }
}
