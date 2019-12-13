#include <math.h>
#include <algorithm>
#include <vector>
#include <omp.h>
#include "../../differentiation_tools.hpp"
#include "../../functions.hpp"
#include "../../functions_sw.hpp"

void Compute_Alt_PE_Transport(
        std::vector<double> & out_array, 
        const std::vector<double> & u_tilde,
        const std::vector<double> & v_tilde,
        const std::vector<double> & h_bar,
        const std::vector<double> & rho,
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

    const int num_pts = Ntime * Ndepth * Nlat * Nlon;

    int Itime, Idepth, Ilat, Ilon, index, ind_UL, ind_LL;

    std::vector<double> temp_array_u(num_pts);
    std::vector<double> temp_array_v(num_pts);
    #pragma omp parallel \
    default(none) \
    shared(u_tilde, v_tilde, h_bar, rho, \
            temp_array_u, temp_array_v, Ntime, Ndepth, Nlat, Nlon ) \
    private(Itime, Ilat, Ilon, index, ind_UL, ind_LL)
    {
        #pragma omp for collapse(3) schedule(static)
        for (Itime = 0; Itime < Ntime; ++Itime) {
            for (Ilat = 0; Ilat < Nlat; ++Ilat) {
                for (Ilon = 0; Ilon < Nlon; ++Ilon) {

                    ind_UL = Index(Itime, 0,      Ilat, Ilon,
                                   Ntime, Ndepth, Nlat, Nlon);
                    ind_LL = Index(Itime, 1,      Ilat, Ilon,
                                   Ntime, Ndepth, Nlat, Nlon);

                    temp_array_u.at(ind_UL) = 
                        constants::g * rho.at(0) 
                            * h_bar.at(ind_UL) * h_bar.at(ind_LL)
                            * ( u_tilde.at(ind_UL) + u_tilde.at(ind_LL) );

                    temp_array_v.at(ind_UL) = 
                        constants::g * rho.at(0) 
                            * h_bar.at(ind_UL) * h_bar.at(ind_LL)
                            * ( v_tilde.at(ind_UL) + v_tilde.at(ind_LL) );
                }
            }
        }
    }

    //
    std::vector<const std::vector<double>*> deriv_fields;
    deriv_fields.push_back(&temp_array_u);
    deriv_fields.push_back(&temp_array_v);

    std::vector<double*> x_deriv_vals, y_deriv_vals, z_deriv_vals;
    double dtudx, dtvdy;
    
    //
    #pragma omp parallel \
    default(none) \
    shared(temp_array_u, temp_array_v, deriv_fields, out_array, \
            mask, longitude, latitude, Ntime, Ndepth, Nlat, Nlon ) \
    private(Itime, Idepth, Ilat, Ilon, index,\
            x_deriv_vals, y_deriv_vals, z_deriv_vals,\
            dtudx, dtvdy)
    {

        x_deriv_vals.clear();
        x_deriv_vals.push_back(&dtudx);
        x_deriv_vals.push_back(NULL);

        y_deriv_vals.clear();
        y_deriv_vals.push_back(NULL);
        y_deriv_vals.push_back(&dtvdy);

        z_deriv_vals.clear();
        z_deriv_vals.push_back(NULL);
        z_deriv_vals.push_back(NULL);

        #pragma omp for collapse(1) schedule(static)
        for (index = 0; index < num_pts; ++index) {

            Index1to4(index, Itime, Idepth, Ilat, Ilon,
                             Ntime, Ndepth, Nlat, Nlon);

            // We'll need some derivatives
            Cart_derivatives_at_point(
                    x_deriv_vals, y_deriv_vals,
                    z_deriv_vals, deriv_fields,
                    latitude, longitude,
                    Itime, Idepth, Ilat, Ilon,
                    Ntime, Ndepth, Nlat, Nlon,
                    mask);

            // (u_i * field)_,i
            out_array.at(index) = dtudx + dtvdy;

        }
    }
}
