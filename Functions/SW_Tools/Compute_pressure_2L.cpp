#include <math.h>
#include <algorithm>
#include <vector>
#include <omp.h>
#include "../../functions.hpp"
#include "../../functions_sw.hpp"

void Compute_pressure_2L(
        std::vector<double> & pressure, 
        const std::vector<double> & h, 
        const std::vector<double> & rho,
        const int & Ntime,
        const int & Ndepth,
        const int & Nlat,
        const int & Nlon
        ){

    const int num_pts = Ntime * Ndepth * Nlat * Nlon;

    int Itime, Idepth, Ilat, Ilon, index, ind_UL, ind_LL;

    #pragma omp parallel \
    default(none) \
    shared(h, pressure, rho, Ntime, Ndepth, Nlat, Nlon) \
    private(Itime, Idepth, Ilat, Ilon, index, ind_UL, ind_LL)
    {
        #pragma omp for collapse(1) schedule(static)
        for (index = 0; index < num_pts; index++) {

            Index1to4(index, Itime, Idepth, Ilat, Ilon,
                             Ntime, Ndepth, Nlat, Nlon);

            ind_UL = Index(Itime, 0,      Ilat, Ilon, 
                           Ntime, Ndepth, Nlat, Nlon);
            ind_LL = Index(Itime, 1,      Ilat, Ilon, 
                           Ntime, Ndepth, Nlat, Nlon);

            if (Idepth == 0) {
                pressure.at(ind_UL) += constants::g * rho.at(0) * h.at(index);
                pressure.at(ind_LL) += constants::g * rho.at(0) * h.at(index);
            } else if (Idepth == 1) {
                pressure.at(ind_UL) += constants::g * rho.at(0) * h.at(index);
                pressure.at(ind_LL) += constants::g * rho.at(1) * h.at(index);
            }
        }
    }
}
