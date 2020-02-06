#include <math.h>
#include <algorithm>
#include <vector>
#include <omp.h>
#include "../../differentiation_tools.hpp"
#include "../../functions.hpp"
#include "../../functions_sw.hpp"

void Compute_Gamma_SW(
        std::vector<double> & out_array, 
        const std::vector<double> & h_bar, 
        const std::vector<double> & u_bar,
        const std::vector<double> & v_bar,
        const std::vector<double> & u_tilde,
        const std::vector<double> & v_tilde,
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

    //   g * rho * h * tau( u_i, h)
    // = g * rho * h * ( h (ui_tilde - ui_bar) )_,i

    const int num_pts = Ntime * Ndepth * Nlat * Nlon;
    std::vector<double> field(num_pts, 0.0);

    int Itime, Idepth, Ilat, Ilon, index, ii;

    //
    std::vector<const std::vector<double>*> deriv_fields;
    deriv_fields.push_back(&field);

    std::vector<double*> x_deriv_vals, y_deriv_vals, z_deriv_vals;
    double dfdx, dfdy, ui_bar, ui_tilde, temp;
    
    //
    #pragma omp parallel \
    default(none) \
    shared(u_tilde, v_tilde, h_bar, u_bar, v_bar, rho, \
            deriv_fields, out_array, mask, longitude, latitude, field, \
            Ntime, Ndepth, Nlat, Nlon ) \
    private(Itime, Idepth, Ilat, Ilon, index, ii, \
            x_deriv_vals, y_deriv_vals, z_deriv_vals, \
            dfdx, dfdy, ui_bar, ui_tilde, temp)
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
        ui_bar   = 0.;
        ui_tilde = 0.;
        temp     = 0.;

        for (ii = 0; ii < 2; ii++) {

            // Compute field = rho * h * u_i * tau(u_i, u_j)
            #pragma omp for collapse(1) schedule(static)
            for (index = 0; index < num_pts; index++) {

                Index1to4(index, Itime, Idepth, Ilat, Ilon,
                                 Ntime, Ndepth, Nlat, Nlon);

                switch (ii) {
                    case 0:
                        ui_bar   = u_bar.at(index);
                        ui_tilde = u_tilde.at(index);
                        break;
                    case 1:
                        ui_bar   = v_bar.at(index);
                        ui_tilde = v_tilde.at(index);
                        break;
                }

                field.at(index) = h_bar.at(index) * ( ui_tilde - ui_bar );
            }

            // compute field_,j
            #pragma omp for collapse(1) schedule(static)
            for (index = 0; index < num_pts; index++) {

                Index1to4(index, Itime, Idepth, Ilat, Ilon,
                                 Ntime, Ndepth, Nlat, Nlon);

                if (ii == 0) {
                    x_deriv_vals.at(0) = &dfdx;
                    y_deriv_vals.at(0) = NULL;
                } else if (ii == 1) {
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

                //   g * rho * h * tau( u_i, h)
                // = g * rho * h * ( h (ui_tilde - ui_bar) )_,i
                if (ii == 0) {
                    temp = constants::g * rho.at(Idepth) * h_bar.at(index) * dfdx;
                } else if (ii == 1) {
                    temp = constants::g * rho.at(Idepth) * h_bar.at(index) * dfdy;
                }
                out_array.at(index) += temp;
            }
        }
    }
}
