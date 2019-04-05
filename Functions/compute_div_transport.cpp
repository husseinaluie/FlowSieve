#include "../constants.hpp"
#include "../functions.hpp"
#include "../differentiation_tools.hpp"
#include <algorithm>
#include <vector>
#include <omp.h>
#include <math.h>

//  
//  J_transport =   0.5 * rho0 * | u_l |^2 * u_l
//                + P_l * u_l
//                - nu * 0.5 * rho * grad( | u_l |^2 )
//                + rho0 * u_l * tau(u_l, u_l)
//
//
//  (Spherical)
//     div(J) = ( 
//                (1 / (r * cos(lat)) ) d/dlon (J_lon),
//                (1 / (r * cos(lat)) ) d/dlat (J_lat * cos(lat)),
//                (1 /  r^2           ) d/dr   (J_r * r^2)
//              )
//
//  (Cartesian)
//     div(J) = ( 
//                d/dx (J_x),
//                d/dy (J_y),
//                d/dz (J_z) 
//              )
//
//
//  Term 1: 0.5 * rho0 * | u_l |^2 * u_l
//      This is advection of large-scale KE by the large-scale 
//      velocity. 
//
//      (index form: 0.5 * rho0 * [ (u_i*u_i) * u_j ]    )
//      (   of grad: 0.5 * rho0 * [ (u_i*u_i) * u_j ],j  )
//               = 0.5 * rho0 * (u_i*u_i),j * u_j 
//               =       rho0 * u_i * u_i,j * u_j
//
//
//
//  Term 2: P_l * u_l
//      Transport caused by pressure 
//
//      (index form:  p * u_j     )
//      (   of grad: (p * u_j),j  )
//
//
//
//  Term 3: - nu * 0.5 * rho * grad( | u_l |^2 )
//      This is diffusion
//      NOT YET IMPLEMENTED
//
//
//
//  Term 4: rho0 * u_l * tau(u, u)
//      This is advection of the large-scale KE
//      by the small-scale flow.
//
//      (index form:  rho0 *   u_i * tau_ij     )
//      (   of grad:  rho0 * [ u_i * tau_ij ],j )
//        = rho0 * ( u_i,j * tau_ij  + u_i * tau_ij,j )
//        = rho0 * 
//            (  
//               u_i,j * ( bar(u_i*u_j)   - bar(u_i  )*bar(u_j) )
//             + u_i   * ( bar(u_i*u_j),j - bar(u_i,j)*bar(u_j) )
//            )
//

void compute_div_transport(
        std::vector<double> & div_J,
        const std::vector<double> & u_x,
        const std::vector<double> & u_y,
        const std::vector<double> & u_z,
        const std::vector<double> & uxux,
        const std::vector<double> & uxuy,
        const std::vector<double> & uxuz,
        const std::vector<double> & uyuy,
        const std::vector<double> & uyuz,
        const std::vector<double> & uzuz,
        const std::vector<double> & coarse_p,
        const std::vector<double> & longitude,
        const std::vector<double> & latitude,
        const int Ntime,
        const int Ndepth,
        const int Nlat,
        const int Nlon,
        const std::vector<double> & mask
        ) {

    double div_J_tmp; 
    
    #if COMP_BC_TRANSFERS
    double dpdx, dpdy, dpdz;
    #endif

    int Ilat, Ilon, index, mask_index;

    double ux, uy, uz;

    double ux_x, uy_x, uz_x;
    double ux_y, uy_y, uz_y;
    double ux_z, uy_z, uz_z;

    double uxux_x, uxuy_y, uxuz_z;
    double uyux_x, uyuy_y, uyuz_z;
    double uzux_x, uzuy_y, uzuz_z;

    double uxux_loc, uxuy_loc, uxuz_loc;
    double uyux_loc, uyuy_loc, uyuz_loc;
    double uzux_loc, uzuy_loc, uzuz_loc;
    
    for (int Itime = 0; Itime < Ntime; Itime++) {
        for (int Idepth = 0; Idepth < Ndepth; Idepth++) {
            #pragma omp parallel \
            default(none) \
            shared( div_J, Idepth, Itime, stdout,\
                    latitude, longitude, mask, \
                    u_x, u_y, u_z, uxux, uxuy, uxuz,\
                    uyuy, uyuz, uzuz)\
            private(Ilat, Ilon, index, mask_index, \
                    ux,   uy,   uz,\
                    ux_x, uy_x, uz_x,\
                    ux_y, uy_y, uz_y,\
                    ux_z, uy_z, uz_z,\
                    uxux_x,   uxuy_y,   uxuz_z,\
                    uyux_x,   uyuy_y,   uyuz_z,\
                    uzux_x,   uzuy_y,   uzuz_z,\
                    uxux_loc, uxuy_loc, uxuz_loc,\
                    uyux_loc, uyuy_loc, uyuz_loc,\
                    uzux_loc, uzuy_loc, uzuz_loc,\
                    div_J_tmp)
            {
                #pragma omp for collapse(2) schedule(guided)
                for (Ilat = 0; Ilat < Nlat; Ilat++) {
                    for (Ilon = 0; Ilon < Nlon; Ilon++) {

                        div_J_tmp = 0.;

                        // Convert our four-index to a one-index
                        index = Index(Itime, Idepth, Ilat, Ilon,
                                Ntime, Ndepth, Nlat, Nlon);
                        mask_index = Index(0,     0,      Ilat, Ilon,
                                Ntime, Ndepth, Nlat, Nlon);

                        if (mask.at(mask_index) == 1) { // Skip land areas

                            // u_i
                            ux = u_x.at(index);
                            uy = u_y.at(index);
                            uz = u_z.at(index);

                            // u_i,j
                            ux_x = Cart_derivative_at_point(u_x, latitude, longitude, "x",
                                    Itime, Idepth, Ilat, Ilon,
                                    Ntime, Ndepth, Nlat, Nlon,
                                    mask);
                            uy_x = Cart_derivative_at_point(u_y, latitude, longitude, "x",
                                    Itime, Idepth, Ilat, Ilon,
                                    Ntime, Ndepth, Nlat, Nlon,
                                    mask);
                            uz_x = Cart_derivative_at_point(u_z, latitude, longitude, "x",
                                    Itime, Idepth, Ilat, Ilon,
                                    Ntime, Ndepth, Nlat, Nlon,
                                    mask);
                            ux_y = Cart_derivative_at_point(u_x, latitude, longitude, "y",
                                    Itime, Idepth, Ilat, Ilon,
                                    Ntime, Ndepth, Nlat, Nlon,
                                    mask);
                            uy_y = Cart_derivative_at_point(u_y, latitude, longitude, "y",
                                    Itime, Idepth, Ilat, Ilon,
                                    Ntime, Ndepth, Nlat, Nlon,
                                    mask);
                            uz_y = Cart_derivative_at_point(u_z, latitude, longitude, "y",
                                    Itime, Idepth, Ilat, Ilon,
                                    Ntime, Ndepth, Nlat, Nlon,
                                    mask);
                            ux_z = Cart_derivative_at_point(u_x, latitude, longitude, "z",
                                    Itime, Idepth, Ilat, Ilon,
                                    Ntime, Ndepth, Nlat, Nlon,
                                    mask);
                            uy_z = Cart_derivative_at_point(u_y, latitude, longitude, "z",
                                    Itime, Idepth, Ilat, Ilon,
                                    Ntime, Ndepth, Nlat, Nlon,
                                    mask);
                            uz_z = Cart_derivative_at_point(u_z, latitude, longitude, "z",
                                    Itime, Idepth, Ilat, Ilon,
                                    Ntime, Ndepth, Nlat, Nlon,
                                    mask);

                            // Advection by coarse velocity field
                            //    0.5 * rho0 * [ (u_i*u_i) * u_j ],j
                            //  =       rho0 * u_i * u_i,j * u_j
                            div_J_tmp += constants::rho0 *
                                ( // j across, i down
                                    ux*ux_x*ux  +  ux*ux_y*uy + ux*ux_z*uz
                                  + uy*uy_x*ux  +  uy*uy_y*uy + uy*uy_z*uz
                                  + uz*uz_x*ux  +  uz*uz_y*uy + uz*uz_z*uz
                                );

                            // Advection by small scale velocity field
                            // rho0 * [ u_i * tau_ij ],j
                            // = rho0 (   u_i,j * ( bar(u_i*u_j)   - bar(u_i  )*bar(u_j) )
                            //          + u_i   * ( bar(u_i*u_j),j - bar(u_i,j)*bar(u_j) )
                            //        )

                            // u_iu_j
                            uxux_loc = uxux.at(index);
                            uxuy_loc = uxuy.at(index);
                            uxuz_loc = uxuz.at(index);
                            uyux_loc = uxuy.at(index);
                            uyuy_loc = uyuy.at(index);
                            uyuz_loc = uyuz.at(index);
                            uzux_loc = uxuz.at(index);
                            uzuy_loc = uyuz.at(index);
                            uzuz_loc = uzuz.at(index);

                            // tau_ij,j
                            uxux_x = Cart_derivative_at_point(uxux, latitude, longitude, "x",
                                    Itime, Idepth, Ilat, Ilon,
                                    Ntime, Ndepth, Nlat, Nlon,
                                    mask);
                            uyux_x = Cart_derivative_at_point(uxuy, latitude, longitude, "x",
                                    Itime, Idepth, Ilat, Ilon,
                                    Ntime, Ndepth, Nlat, Nlon,
                                    mask);
                            uzux_x = Cart_derivative_at_point(uxuz, latitude, longitude, "x",
                                    Itime, Idepth, Ilat, Ilon,
                                    Ntime, Ndepth, Nlat, Nlon,
                                    mask);
                            uxuy_y = Cart_derivative_at_point(uxuy, latitude, longitude, "y",
                                    Itime, Idepth, Ilat, Ilon,
                                    Ntime, Ndepth, Nlat, Nlon,
                                    mask);
                            uyuy_y = Cart_derivative_at_point(uyuy, latitude, longitude, "y",
                                    Itime, Idepth, Ilat, Ilon,
                                    Ntime, Ndepth, Nlat, Nlon,
                                    mask);
                            uzuy_y = Cart_derivative_at_point(uyuz, latitude, longitude, "y",
                                    Itime, Idepth, Ilat, Ilon,
                                    Ntime, Ndepth, Nlat, Nlon,
                                    mask);
                            uxuz_z = Cart_derivative_at_point(uxuz, latitude, longitude, "z",
                                    Itime, Idepth, Ilat, Ilon,
                                    Ntime, Ndepth, Nlat, Nlon,
                                    mask);
                            uyuz_z = Cart_derivative_at_point(uyuz, latitude, longitude, "z",
                                    Itime, Idepth, Ilat, Ilon,
                                    Ntime, Ndepth, Nlat, Nlon,
                                    mask);
                            uzuz_z = Cart_derivative_at_point(uzuz, latitude, longitude, "z",
                                    Itime, Idepth, Ilat, Ilon,
                                    Ntime, Ndepth, Nlat, Nlon,
                                    mask);

                            // rho0 * ( u_i,j * ( bar(u_i*u_j)   - bar(u_i  )*bar(u_j) ) )
                            div_J_tmp += constants::rho0 *
                                ( // j across, i down
                                   ux_x * ( uxux_loc - ux*ux ) + ux_y * ( uxuy_loc - ux*uy ) + ux_z * ( uxuz_loc - ux*uz )
                                 + uy_x * ( uyux_loc - uy*ux ) + uy_y * ( uyuy_loc - uy*uy ) + uy_z * ( uyuz_loc - uy*uz )
                                 + uz_x * ( uzux_loc - uz*ux ) + uz_y * ( uzuy_loc - uz*uy ) + uz_z * ( uzuz_loc - uz*uz )
                                );
                            // rho0 * ( u_i   * ( bar(u_i*u_j),j - bar(u_i,j)*bar(u_j) ) )
                            div_J_tmp += constants::rho0 *
                                ( // j across, i down
                                   ux * ( (uxux_x - ux_x*ux) + (uxuy_y - ux_y*uy) + ( uxuz_z - ux_z*uz) )
                                 + uy * ( (uyux_x - uy_x*ux) + (uyuy_y - uy_y*uy) + ( uyuz_z - uy_z*uz) )
                                 + uz * ( (uzux_x - uz_x*ux) + (uzuy_y - uz_y*uy) + ( uzuz_z - uz_z*uz) )
                                );

                        }

                        //
                        div_J.at(index) = div_J_tmp;

                    } // end Ilon loop
                } // end Ilat loop
            } // end pragma

            // Pressure term
            #if COMP_BC_TRANSFERS
            #pragma omp parallel \
            default(none) \
            shared( div_J, Idepth, Itime, stdout,\
                    latitude, longitude, mask, coarse_p, \
                    u_x, u_y, u_z)\
            private(Ilat, Ilon, index, mask_index, \
                    dpdx, dpdy, dpdz,\
                    ux,   uy,   uz,\
                    div_J_tmp)
            {
                #pragma omp for collapse(2) schedule(guided)
                for (Ilat = 0; Ilat < Nlat; Ilat++) {
                    for (Ilon = 0; Ilon < Nlon; Ilon++) {

                        div_J_tmp = 0.;

                        // Convert our four-index to a one-index
                        index = Index(Itime, Idepth, Ilat, Ilon,
                                Ntime, Ndepth, Nlat, Nlon);
                        mask_index = Index(0,     0,      Ilat, Ilon,
                                Ntime, Ndepth, Nlat, Nlon);

                        if (mask.at(mask_index) == 1) { // Skip land areas

                            // u_i
                            ux = u_x.at(index);
                            uy = u_y.at(index);
                            uz = u_z.at(index);

                            // p_,j
                            dpdx = Cart_derivative_at_point(
                                    coarse_p,
                                    latitude, longitude, "x",
                                    Itime, Idepth, Ilat, Ilon,
                                    Ntime, Ndepth, Nlat, Nlon,
                                    mask);
                            dpdy = Cart_derivative_at_point(
                                    coarse_p,
                                    latitude, longitude, "y",
                                    Itime, Idepth, Ilat, Ilon,
                                    Ntime, Ndepth, Nlat, Nlon,
                                    mask);
                            dpdz = Cart_derivative_at_point(
                                    coarse_p,
                                    latitude, longitude, "z",
                                    Itime, Idepth, Ilat, Ilon,
                                    Ntime, Ndepth, Nlat, Nlon,
                                    mask);

                            // Pressure term
                            // (p * u_j),j = p_,j * u_j
                            div_J_tmp += ux * dpdx + uy * dpdy + uz * dpdz;

                        }

                        //
                        div_J.at(index) += div_J_tmp;

                    } // end Ilon loop
                } // end Ilat loop
            } // end pragma
            #endif
        } // end Idepth loop
    } // end Itime loop
}
