#include <math.h>
#include <algorithm>
#include <vector>
#include <omp.h>
#include "../../functions.hpp"
#include "../../functions_sw.hpp"

void Compute_Viscous_Loss_PE(
        std::vector<double> & out_array, 
        const std::vector<double> & h_bar, 
        const std::vector<double> & lap_h_bar, 
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

    int Itime, Ilat, Ilon, ind_UL, ind_LL;

    //
    #pragma omp parallel \
    default(none) \
    shared(h_bar, lap_h_bar, rho, out_array,\
            Ntime, Ndepth, Nlat, Nlon) \
    private(Itime, Ilat, Ilon, ind_UL, ind_LL)
    {
        #pragma omp for collapse(3) schedule(static)
        for (Itime = 0; Itime < Ntime; Itime++) {
            for (Ilat = 0; Ilat < Nlat; Ilat++) {
                for (Ilon = 0; Ilon < Nlon; Ilon++) {

                    ind_UL = Index(Itime, 0,      Ilat, Ilon,
                                   Ntime, Ndepth, Nlat, Nlon);
                    ind_LL = Index(Itime, 1,      Ilat, Ilon,
                                   Ntime, Ndepth, Nlat, Nlon);

                    out_array.at(ind_UL) = 
                        constants::g * rho.at(0) * nu * (
                                ( h_bar.at(ind_UL) + h_bar.at(ind_LL) ) 
                                    * lap_h_bar.at(ind_UL)        
                                +
                                h_bar.at(ind_UL) * lap_h_bar.at(ind_LL)
                            );

                    out_array.at(ind_LL) = 
                        constants::g * rho.at(1) * nu 
                            * h_bar.at(ind_LL) * lap_h_bar.at(ind_LL);

                }
            }
        }
    }
}
