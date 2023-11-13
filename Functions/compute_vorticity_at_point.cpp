#include <algorithm>
#include <vector>
#include <math.h>
#include "../functions.hpp"
#include "../differentiation_tools.hpp"
#include "../constants.hpp"

/*!
 * \brief Compute the vorticity, divergnce, and OkuboWeiss at a point
 *
 * Since we're computing derivatives anyways, might as well get a few things out of it.
 *
 *
 * At the moment, only computes the vort_r component of vorticity
 *
 * @param[in,out]   vort_r_tmp,vort_lon_tmp,vort_lat_tmp    where to store vorticity components
 * @param[in,out]   div_tmp                                 where to store velocity divergence
 * @param[in,out]   OkuboWeiss_tmp                          where to store OkuboWeiss result
 * @param[in]       source_data                     dataset class instance containing data (Psi, Phi, etc)
 * @param[in]       u_r,u_lon,u_lat                         velocity components
 * @param[in]       Itime,Idepth,Ilat,Ilon                  Current index in time and space
 *
 */
void compute_vorticity_at_point(
        double & vort_r_tmp,
        double & vort_lon_tmp,
        double & vort_lat_tmp,
        double & div_tmp,
        double & OkuboWeiss_tmp,
        double & cyclonic_energy,
        double & anticyclonic_energy,
        double & divergent_strain_energy,
        double & traceless_strain_energy,
        const dataset & source_data,
        const std::vector<double> & u_r,
        const std::vector<double> & u_lon,
        const std::vector<double> & u_lat,
        const int Itime,
        const int Idepth,
        const int Ilat,
        const int Ilon
        ) {

    const int   Ntime   = source_data.Ntime,    // this is the MPI-local Ntime, not the full Ntime
                Ndepth  = source_data.Ndepth,   // this is the MPI-local Ndepth, not the full Ndepth
                Nlat    = source_data.Nlat,
                Nlon    = source_data.Nlon;

    // For the moment, only compute vort_r
    vort_r_tmp   = 0.;
    vort_lon_tmp = 0.;
    vort_lat_tmp = 0.;

    std::vector<const std::vector<double>*> deriv_fields {&u_lon, &u_lat, &u_r};

    if (constants::CARTESIAN) {

        double ux_x, ux_y, ux_z,
               uy_x, uy_y, uy_z,
               uz_x, uz_y, uz_z;

        std::vector<double*>    x_deriv_vals {&ux_x, &uy_x, &uz_x},
                                y_deriv_vals {&ux_y, &uy_y, &uz_y},
                                z_deriv_vals {&ux_z, &uy_z, &uz_z};

        Cart_derivatives_at_point(
           x_deriv_vals, y_deriv_vals, z_deriv_vals, deriv_fields,
           source_data, Itime, Idepth, Ilat, Ilon);

        vort_lon_tmp = uz_y - uy_z;
        vort_lat_tmp = ux_z - uz_x;
        vort_r_tmp   = uy_x - ux_y;

        div_tmp = ux_x + uy_y + uz_z;

        OkuboWeiss_tmp = pow(ux_x - uy_y, 2) + 4 * ux_y * uy_x;
    } else {

        const std::vector<double>   &latitude   = source_data.latitude,
                                    &longitude  = source_data.longitude;

        const std::vector<bool> &mask = source_data.mask;

        size_t index = Index(Itime, Idepth, Ilat, Ilon,
                             Ntime, Ndepth, Nlat, Nlon);
        double ur_r,   ur_lon,   ur_lat, 
               ulon_r, ulon_lon, ulon_lat, 
               ulat_r, ulat_lon, ulat_lat;

        std::vector<double*>    lon_deriv_vals {&ulon_lon, &ulat_lon, &ur_lon},
                                lat_deriv_vals {&ulon_lat, &ulat_lat, &ur_lat},
                                r_deriv_vals   {&ulon_r,   &ulat_r,   &ur_r  };

        spher_derivative_at_point( lat_deriv_vals, deriv_fields, latitude, "lat",
                Itime, Idepth, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon, mask);

        spher_derivative_at_point( lon_deriv_vals, deriv_fields, longitude, "lon",
                Itime, Idepth, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon, mask);

        //spher_derivative_at_point( r_deriv_vals, deriv_fields, depth, "depth",
        //        Itime, Idepth, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon, mask);
        // Don't compute depth derivatives. It requires merging across MPI ranks.
        ur_r   = 0.;
        ulon_r = 0.;
        ulat_r = 0.;

        const double    lat       = latitude.at(Ilat),
                        cos_lat   = cos(lat),
                        tan_lat   = tan(lat),
                        u_r_loc   = u_r.at(index),
                        u_lon_loc = u_lon.at(index),
                        u_lat_loc = u_lat.at(index);

        //
        //// First, do vorticity
        //
        vort_r_tmp   = ( ulat_lon / cos_lat - ulon_lat + tan_lat * u_lon_loc ) / ( constants::R_earth );
        vort_lon_tmp = ( ur_lat - u_lat_loc ) / ( constants::R_earth ) - ulat_r;
        vort_lat_tmp = ( u_lon_loc - ur_lon / cos_lat ) / ( constants::R_earth ) + ulon_r;

        //
        //// Now the divergence
        //
        div_tmp =   ( 2. * u_r_loc / constants::R_earth )
                  + ( ur_r )
                  + ( ulon_lon / ( constants::R_earth * cos_lat ) )
                  + ( ulat_lat / constants::R_earth )
                  - ( u_lat_loc * tan_lat / constants::R_earth );

        //
        //// Now the Okubo-Weiss parameter
        //
        const double    sec_lat = 1. / cos_lat,
                        S_11 = (sec_lat * ulon_lon - tan_lat * u_lat_loc) / constants::R_earth,
                        S_22 = ulat_lat / constants::R_earth,
                        S_12 = 0.5 * ( sec_lat * ulat_lon + tan_lat * u_lon_loc + ulon_lat ) / constants::R_earth;

        const double    s_n = ( S_11 - S_22 ),
                        s_s = 2 * S_12;

        //const double    s_n = ( cos_lat * ulon_lon - ulat_lat ) / constants::R_earth,
        //                s_s = ( cos_lat * ulat_lon + ulon_lat ) / constants::R_earth;
        OkuboWeiss_tmp = pow(s_n, 2) + pow(s_s, 2) - pow(vort_r_tmp, 2);


        // Johnson-decomposed Energy
        if ( vort_r_tmp * lat >= 0 ) {
            cyclonic_energy = 2 * pow(0.5 * vort_r_tmp, 2);
            anticyclonic_energy = 0.;
        } else {
            cyclonic_energy = 0.;
            anticyclonic_energy = 2 * pow(0.5 * vort_r_tmp, 2);
        }
        divergent_strain_energy = pow( S_11, 2) + pow( S_22, 2 );
        traceless_strain_energy = 2 * pow( S_12, 2);

    }
}
