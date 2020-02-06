#include <math.h>
#include <algorithm>
#include <vector>
#include <omp.h>
#include "../../differentiation_tools.hpp"
#include "../../functions.hpp"
#include "../../functions_sw.hpp"

void Compute_Pi_SW(
        std::vector<double> & out_array, 
        const std::vector<double> & h_bar, 
        const std::vector<double> & u_tilde,
        const std::vector<double> & v_tilde,
        const std::vector<double> & uu_tilde,
        const std::vector<double> & uv_tilde,
        const std::vector<double> & vv_tilde,
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

    int Itime, Idepth, Ilat, Ilon, index;

    //
    std::vector<const std::vector<double>*> deriv_fields;
    deriv_fields.push_back(&u_tilde);
    deriv_fields.push_back(&v_tilde);

    std::vector<double*> x_deriv_vals, y_deriv_vals, z_deriv_vals;
    double dudx, dvdx, dudy, dvdy;
    double u_t, v_t, uu_t, uv_t, vv_t;
    
    //
    #pragma omp parallel \
    default(none) \
    shared(u_tilde, v_tilde, h_bar, uu_tilde, uv_tilde, \
            vv_tilde, deriv_fields, out_array, \
            mask, longitude, latitude, rho, \
            Ntime, Ndepth, Nlat, Nlon ) \
    private(Itime, Idepth, Ilat, Ilon, index, \
            x_deriv_vals, y_deriv_vals, z_deriv_vals, \
            dudx, dvdx, dudy, dvdy, u_t, v_t, uu_t, uv_t, vv_t)
    {

        x_deriv_vals.clear();
        x_deriv_vals.push_back(&dudx);
        x_deriv_vals.push_back(&dvdx);

        y_deriv_vals.clear();
        y_deriv_vals.push_back(&dudy);
        y_deriv_vals.push_back(&dvdy);

        z_deriv_vals.clear();
        z_deriv_vals.push_back(NULL);
        z_deriv_vals.push_back(NULL);

        #pragma omp for collapse(1) schedule(static)
        for (index = 0; index < num_pts; index++) {

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

            u_t  = u_tilde.at( index);
            v_t  = v_tilde.at( index);
            uu_t = uu_tilde.at(index);
            uv_t = uv_tilde.at(index);
            vv_t = vv_tilde.at(index);

            // rho * h_bar * u_tilde_i,j * tau_tilde(u_i,u_j)
            out_array.at(index) = 
                rho.at(Idepth) * h_bar.at(index) * (
                      ( dudx        ) * ( uu_t - u_t * u_t )
                    + ( dudy + dvdx ) * ( uv_t - u_t * v_t )
                    + (        dvdy ) * ( vv_t - v_t * v_t )
                    );

        }
    }
}

