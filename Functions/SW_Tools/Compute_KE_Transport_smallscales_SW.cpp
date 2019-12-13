#include <math.h>
#include <algorithm>
#include <vector>
#include <omp.h>
#include "../../differentiation_tools.hpp"
#include "../../functions.hpp"
#include "../../functions_sw.hpp"

void Compute_KE_Transport_smallscales_SW(
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
    std::vector<double> field(num_pts, 0.0);

    int Itime, Idepth, Ilat, Ilon, index, ii, jj;

    //
    std::vector<const std::vector<double>*> deriv_fields;
    deriv_fields.push_back(&field);

    std::vector<double*> x_deriv_vals, y_deriv_vals, z_deriv_vals;
    double dfdx, dfdy, ui, tau_uiuj;
    
    //
    #pragma omp parallel \
    default(none) \
    shared(u_tilde, v_tilde, h_bar, uu_tilde, uv_tilde, rho, \
            vv_tilde, deriv_fields, out_array, field, \
            mask, longitude, latitude, Ntime, Ndepth, Nlat, Nlon ) \
    private(Itime, Idepth, Ilat, Ilon, index, ii, jj, \
            x_deriv_vals, y_deriv_vals, z_deriv_vals, \
            dfdx, dfdy, ui, tau_uiuj)
    {

        x_deriv_vals.clear();
        x_deriv_vals.push_back(&dfdx);

        y_deriv_vals.clear();
        y_deriv_vals.push_back(&dfdy);

        z_deriv_vals.clear();
        z_deriv_vals.push_back(NULL);

        // Zero out
        #pragma omp for collapse(1) schedule(static)
        for (index = 0; index < num_pts; index++) {
            out_array.at(index) = 0.;
        }
        ui = 0.;
        tau_uiuj = 0.;

        for (ii = 0; ii < 2; ii++) {
            for (jj = 0; jj < 2; jj++) {

                // Compute field = rho * h * u_i * tau(u_i, u_j)
                #pragma omp for collapse(1) schedule(static)
                for (index = 0; index < num_pts; index++) {

                    Index1to4(index, Itime, Idepth, Ilat, Ilon,
                                     Ntime, Ndepth, Nlat, Nlon);

                    switch (ii) {
                        case 0:
                            ui = u_tilde.at(index);
                            switch (jj) {
                                case 0:
                                    tau_uiuj = uu_tilde.at(index) - ui * ui;
                                    break;
                                case 1:
                                    tau_uiuj = uv_tilde.at(index) - ui * v_tilde.at(index);
                                    break;
                            }
                            break;
                        case 1:
                            ui = v_tilde.at(index);
                            switch (jj) {
                                case 0:
                                    tau_uiuj = uv_tilde.at(index) - ui * u_tilde.at(index);
                                    break;
                                case 1:
                                    tau_uiuj = vv_tilde.at(index) - ui * ui;
                                    break;
                            }
                            break;
                    }

                    field.at(index) = rho.at(Idepth) * h_bar.at(index) * ui * tau_uiuj;
                }

                // compute field_,j
                #pragma omp for collapse(1) schedule(static)
                for (index = 0; index < num_pts; index++) {

                    Index1to4(index, Itime, Idepth, Ilat, Ilon,
                                     Ntime, Ndepth, Nlat, Nlon);

                    if (jj == 0) {
                        x_deriv_vals.at(0) = &dfdx;
                        y_deriv_vals.at(0) = NULL;
                    } else if (jj == 1) {
                        x_deriv_vals.at(0) = NULL;
                        y_deriv_vals.at(0) = &dfdy;
                    }

                    Cart_derivatives_at_point(
                            x_deriv_vals, y_deriv_vals,
                            z_deriv_vals, deriv_fields,
                            latitude, longitude,
                            Itime, Idepth, Ilat, Ilon,
                            Ntime, Ndepth, Nlat, Nlon,
                            mask);

                    // ( rho * h_bar * u_tilde_i,j * tau_tilde(u_i,u_j) )_,j
                    if (jj == 0) {
                        out_array.at(index) += dfdx;
                    } else if (jj == 1) {
                        out_array.at(index) += dfdy;
                    }
                }
            }
        }
    }
}
