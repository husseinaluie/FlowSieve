#include <math.h>
#include <algorithm>
#include <vector>
#include <omp.h>
#include "../../differentiation_tools.hpp"
#include "../../functions.hpp"
#include "../../functions_sw.hpp"

void Compute_misc1_SW(
        std::vector<double> & out_array, 
        const std::vector<double> & h_bar, 
        const std::vector<double> & u_bar, 
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

    std::vector<double*> x_deriv_vals, y_deriv_vals, z_deriv_vals;
    double dhdx, dhdy;
    
    //
    #pragma omp parallel \
    default(none) \
    shared(u_bar, v_bar, h_bar, rhos, deriv_fields, out_array, \
            mask, longitude, latitude, Ntime, Ndepth, Nlat, Nlon ) \
    private(Itime, Ilat, Ilon, ind_UL, ind_LL, \
            x_deriv_vals, y_deriv_vals, z_deriv_vals,\
            dhdx, dhdy)
    {

        x_deriv_vals.clear();
        x_deriv_vals.push_back(&dhdx);

        y_deriv_vals.clear();
        y_deriv_vals.push_back(&dhdy);

        z_deriv_vals.clear();
        z_deriv_vals.push_back(NULL);

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
                            x_deriv_vals, y_deriv_vals,
                            z_deriv_vals, deriv_fields,
                            latitude, longitude,
                            Itime, 0,      Ilat, Ilon,
                            Ntime, Ndepth, Nlat, Nlon,
                            mask);

                    out_array.at(ind_UL) = 
                        + constants::g * rhos.at(0) * h_bar.at(ind_LL) *
                        (dhdx * u_bar.at(ind_UL) + dhdy * v_bar.at(ind_UL) );

                    out_array.at(ind_LL) = 
                        - constants::g * rhos.at(0) * h_bar.at(ind_LL) *
                        (dhdx * u_bar.at(ind_LL) + dhdy * v_bar.at(ind_LL) );

                }
            }
        }
    }
}
