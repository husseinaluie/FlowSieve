#include <math.h>
#include <algorithm>
#include <vector>
#include <omp.h>
#include "../../differentiation_tools.hpp"
#include "../../functions.hpp"
#include "../../functions_sw.hpp"

void Compute_full_bc_term_PE_parts(
        std::vector<double> & out_array, 
        const std::vector<double> & h_bar, 
        const std::vector<double> & u_tilde, 
        const std::vector<double> & u_bar, 
        const std::vector<double> & v_tilde, 
        const std::vector<double> & v_bar, 
        const std::vector<double> & rhos, 
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

    int Itime, Ilat, Ilon, ind_UL, ind_LL;

    //
    std::vector<const std::vector<double>*> deriv_fields;
    deriv_fields.push_back(&h_bar);

    std::vector<double*> x_derivs_UL, y_derivs_UL, z_derivs_UL,
                         x_derivs_LL, y_derivs_LL, z_derivs_LL;
    double dhdx_UL, dhdy_UL, dhdx_LL, dhdy_LL;
    
    //
    #pragma omp parallel \
    default(none) \
    shared( u_bar, v_bar, h_bar, deriv_fields, out_array, \
            u_tilde, v_tilde, rhos, \
            mask, longitude, latitude, Ntime, Ndepth, Nlat, Nlon ) \
    private( Itime, Ilat, Ilon, ind_UL, ind_LL, \
             x_derivs_UL, y_derivs_UL, z_derivs_UL,\
             x_derivs_LL, y_derivs_LL, z_derivs_LL,\
             dhdx_UL, dhdy_UL, dhdx_LL, dhdy_LL)
    {

        // UL derivs
        x_derivs_UL.clear();
        y_derivs_UL.clear();
        z_derivs_UL.clear();

        x_derivs_UL.push_back(&dhdx_UL);
        y_derivs_UL.push_back(&dhdy_UL);
        z_derivs_UL.push_back(NULL);

        // LL derivs
        x_derivs_LL.clear();
        y_derivs_LL.clear();
        z_derivs_LL.clear();

        x_derivs_LL.push_back(&dhdx_LL);
        y_derivs_LL.push_back(&dhdy_LL);
        z_derivs_LL.push_back(NULL);

        #pragma omp for collapse(3) schedule(static)
        for (Itime = 0; Itime < Ntime; Itime++) {
            for (Ilat = 0; Ilat < Nlat; Ilat++) {
                for (Ilon = 0; Ilon < Nlon; Ilon++) {

                    ind_UL = Index(Itime, 0,      Ilat, Ilon,
                                   Ntime, Ndepth, Nlat, Nlon);

                    ind_LL = Index(Itime, 1,      Ilat, Ilon,
                                   Ntime, Ndepth, Nlat, Nlon);

                    // We'll need some derivatives
                    Cart_derivatives_at_point(
                            x_derivs_UL, y_derivs_UL,
                            z_derivs_UL, deriv_fields,
                            latitude, longitude,
                            Itime, 0,      Ilat, Ilon,
                            Ntime, Ndepth, Nlat, Nlon,
                            mask);

                    Cart_derivatives_at_point(
                            x_derivs_LL, y_derivs_LL,
                            z_derivs_LL, deriv_fields,
                            latitude, longitude,
                            Itime, 1,      Ilat, Ilon,
                            Ntime, Ndepth, Nlat, Nlon,
                            mask);
    
                    out_array.at(ind_LL) = 
                        constants::g * rhos.at(0) * (
                                dhdx_UL * h_bar.at(ind_LL) 
                                        * ( u_tilde.at(ind_LL) - u_bar.at(ind_LL) )

                              + dhdy_UL * h_bar.at(ind_LL) 
                                        * ( v_tilde.at(ind_LL) - v_bar.at(ind_LL) )
                            );
    
                    out_array.at(ind_UL) = 
                        constants::g * rhos.at(0) * (
                                dhdx_LL * h_bar.at(ind_UL) 
                                        * ( u_tilde.at(ind_UL) - u_bar.at(ind_UL) )

                              + dhdy_LL * h_bar.at(ind_UL) 
                                        * ( v_tilde.at(ind_UL) - v_bar.at(ind_UL) )
                            );

                }
            }
        }
    }
}

