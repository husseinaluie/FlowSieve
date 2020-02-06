#include <math.h>
#include <algorithm>
#include <vector>
#include <omp.h>
#include "../../differentiation_tools.hpp"
#include "../../functions.hpp"
#include "../../functions_sw.hpp"

void Compute_misc_conversion_SW(
        std::vector<double> & out_array, 
        const std::vector<double> & h_bar, 
        const std::vector<double> & u_bar, 
        const std::vector<double> & v_bar, 
        const std::vector<double> & rhos,
        const double alpha,
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
    deriv_fields.push_back(&u_bar);
    deriv_fields.push_back(&v_bar);

    std::vector<double*> x_UL_deriv_vals, y_UL_deriv_vals, z_UL_deriv_vals,
                         x_LL_deriv_vals, y_LL_deriv_vals, z_LL_deriv_vals;
    double dhdx_UL, dhdy_UL, dudx_UL, dudy_UL, dvdx_UL, dvdy_UL,
           dhdx_LL, dhdy_LL, dudx_LL, dudy_LL, dvdx_LL, dvdy_LL;
    
    //
    #pragma omp parallel \
    default(none) \
    shared(u_bar, v_bar, h_bar, rhos, deriv_fields, out_array, \
            mask, longitude, latitude, Ntime, Ndepth, Nlat, Nlon ) \
    private(Itime, Ilat, Ilon, ind_UL, ind_LL, \
            x_UL_deriv_vals, y_UL_deriv_vals, z_UL_deriv_vals,\
            x_LL_deriv_vals, y_LL_deriv_vals, z_LL_deriv_vals,\
            dhdx_UL, dhdy_UL, dudx_UL, dudy_UL, dvdx_UL, dvdy_UL,\
            dhdx_LL, dhdy_LL, dudx_LL, dudy_LL, dvdx_LL, dvdy_LL)
    {

        // Derivatives in upper layer
        x_UL_deriv_vals.clear();
        x_UL_deriv_vals.push_back(&dhdx_UL);
        x_UL_deriv_vals.push_back(&dudx_UL);
        x_UL_deriv_vals.push_back(&dvdx_UL);

        y_UL_deriv_vals.clear();
        y_UL_deriv_vals.push_back(&dhdy_UL);
        y_UL_deriv_vals.push_back(&dudy_UL);
        y_UL_deriv_vals.push_back(&dvdy_UL);

        z_UL_deriv_vals.clear();
        z_UL_deriv_vals.push_back(NULL);
        z_UL_deriv_vals.push_back(NULL);
        z_UL_deriv_vals.push_back(NULL);

        // Derivatives in lower layer
        x_LL_deriv_vals.clear();
        x_LL_deriv_vals.push_back(&dhdx_LL);
        x_LL_deriv_vals.push_back(&dudx_LL);
        x_LL_deriv_vals.push_back(&dvdx_LL);

        y_LL_deriv_vals.clear();
        y_LL_deriv_vals.push_back(&dhdy_LL);
        y_LL_deriv_vals.push_back(&dudy_LL);
        y_LL_deriv_vals.push_back(&dvdy_LL);

        z_LL_deriv_vals.clear();
        z_LL_deriv_vals.push_back(NULL);
        z_LL_deriv_vals.push_back(NULL);
        z_LL_deriv_vals.push_back(NULL);

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
                            x_UL_deriv_vals, y_UL_deriv_vals,
                            z_UL_deriv_vals, deriv_fields,
                            latitude, longitude,
                            Itime, 0,      Ilat, Ilon,
                            Ntime, Ndepth, Nlat, Nlon,
                            mask);

                    Cart_derivatives_at_point(
                            x_LL_deriv_vals, y_LL_deriv_vals,
                            z_LL_deriv_vals, deriv_fields,
                            latitude, longitude,
                            Itime, 1,      Ilat, Ilon,
                            Ntime, Ndepth, Nlat, Nlon,
                            mask);

                    out_array.at(ind_UL) = 
                        alpha * constants::g * rhos.at(0) * 
                            (   dudx_UL * dhdx_UL * dhdx_LL
                              + 0.5 * ( dudy_UL + dvdx_UL ) 
                                    * ( dhdx_UL * dhdy_LL + dhdy_UL * dhdx_LL )
                              + dvdy_UL * dhdy_UL * dhdy_LL
                            );

                    out_array.at(ind_LL) = 
                        alpha * constants::g * rhos.at(0) * 
                            (   dudx_LL * dhdx_UL * dhdx_LL
                              + 0.5 * ( dudy_LL + dvdx_LL ) 
                                    * ( dhdx_UL * dhdy_LL + dhdy_UL * dhdx_LL )
                              + dvdy_LL * dhdy_UL * dhdy_LL
                            );

                }
            }
        }
    }
}
