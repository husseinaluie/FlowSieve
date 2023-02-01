#include <vector>
#include <assert.h>
#include "../../constants.hpp"
#include "../../functions.hpp"
#include "../../differentiation_tools.hpp"
#include <math.h>
#include <string>

void Cart_derivatives_at_point(
        const std::vector<double*> & x_deriv_vals,
        const std::vector<double*> & y_deriv_vals,
        const std::vector<double*> & z_deriv_vals,
        const std::vector<const std::vector<double>*> & fields,
        const dataset & source_data,
        const int Itime,
        const int Idepth,
        const int Ilat,
        const int Ilon,
        const int order_of_deriv,
        const int diff_ord,
        const bool include_depth_derivs
        ) {

    const int   Ntime   = source_data.Ntime,    // this is the MPI-local Ntime, not the full Ntime
                Ndepth  = include_depth_derivs ? source_data.full_Ndepth : source_data.Ndepth,   
                Nlat    = source_data.Nlat,
                Nlon    = source_data.Nlon;

    const std::vector<double>   &depth      = source_data.depth,
                                &latitude   = source_data.latitude,
                                &longitude  = source_data.longitude;

    const bool use_depth_mask  = (source_data.use_depth_derivatives         // if the user turned on depth derivatives
                                  and include_depth_derivs                  // and the call to Cart_derivs wants depth derivatives
                                  and ( source_data.Nprocs_in_depth > 1 )   // and we have more than one MPI rank in depth
                                  );
    const std::vector<bool> &mask = use_depth_mask ? source_data.mask_DEPTH : source_data.mask;
    const int Idepth_DEPTH = use_depth_mask ? (Idepth + source_data.myStarts[1]) : Idepth;

    // Confirm that input sizes match
    assert(x_deriv_vals.size() == fields.size());
    assert(y_deriv_vals.size() == fields.size());
    assert(z_deriv_vals.size() == fields.size());
    const int num_deriv = x_deriv_vals.size();

    // Create arrays to store the derivatives.
    //  the "_p" arrays are pointers to the values
    //  in the other arrays
    std::vector<double*> dfields_dlon_p(num_deriv);
    std::vector<double>  dfields_dlon(  num_deriv, 0.);

    std::vector<double*> dfields_dlat_p(num_deriv);
    std::vector<double>  dfields_dlat(  num_deriv, 0.);

    std::vector<double*> dfields_dr_p(num_deriv);
    std::vector<double>  dfields_dr(  num_deriv, 0.);

    for (int ii = 0; ii < num_deriv; ii++) {
        dfields_dlon_p.at(ii) = &(dfields_dlon.at(ii));
        dfields_dlat_p.at(ii) = &(dfields_dlat.at(ii));
        dfields_dr_p.at(  ii) = &(dfields_dr.at(  ii));
    }

    // Compute spherical derivatives
    spher_derivative_at_point(
        dfields_dlon_p, fields, longitude, "lon",
        Itime, Idepth_DEPTH, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon,
        mask, order_of_deriv, diff_ord);

    spher_derivative_at_point(
        dfields_dlat_p, fields, latitude, "lat",
        Itime, Idepth_DEPTH, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon,
        mask, order_of_deriv, diff_ord);

    if ( include_depth_derivs ) {
        // At the moment, can only use second order for
        // depth derivatives since it's non-uniform.
        spher_derivative_at_point(
            dfields_dr_p, fields, depth, "depth",
            Itime, Idepth_DEPTH, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon,
            mask, order_of_deriv, 2);
        if (not(source_data.depth_is_elevation)) {
            // If we have an actual depth grid, multiply by
            // (-1) to account for the grid increasing down, not up
            for (int ii = 0; ii < num_deriv; ii++) {
                dfields_dr.at(ii) = - dfields_dr.at(ii);
            }
        }
    }

    // Define conversion coefficients
    double cx_lon, cx_lat, cx_r,
           cy_lon, cy_lat, cy_r,
           cz_lon, cz_lat, cz_r;
    if (constants::CARTESIAN) {
        cx_r   = 0.;
        cx_lon = 1.;
        cx_lat = 0.;

        cy_r   = 0.;
        cy_lon = 0.;
        cy_lat = 1.;

        cz_r   = 1.;
        cz_lon = 0.;
        cz_lat = 0.;
    } else {
        double lon = longitude.at(Ilon);
        double lat = latitude.at(Ilat);
        double r = constants::R_earth;

        double cos_lat = cos(lat);
        double cos_lon = cos(lon);
        double sin_lat = sin(lat);
        double sin_lon = sin(lon);

        // ddx =   ( cos(lon) * cos(lat)                ) * ddr
        //       - ( sin(lon)            / (r cos(lat)) ) * ddlon   
        //       - ( cos(lon) * sin(lat) /  r           ) * ddlat
        cx_r   =     cos_lon  * cos_lat;
        cx_lon = -   sin_lon             / (r * cos_lat );
        cx_lat = -   cos_lon  * sin_lat  /  r;

        // ddy =   ( sin(lon) * cos(lat)                ) * ddr
        //         ( cos(lon)            / (r cos(lat)) ) * ddlon   
        //       - ( sin(lon) * sin(lat) /  r           ) * ddlat
        cy_r   =     sin_lon  * cos_lat;
        cy_lon =     cos_lon             / (r * cos_lat );
        cy_lat = -   sin_lon  * sin_lat  /  r;

        // ddz =   (            sin(lat)                ) * ddr
        //         (            0.                      ) * ddlon
        //         (            cos(lat) /  r           ) * ddlat
        cz_r   =                sin_lat;
        cz_lon =                0.;
        cz_lat =                cos_lat  /  r;
    }

    // Now do the actual conversion
    for (int ii = 0; ii < num_deriv; ii++) {

        // Compute x derivatives
        if (x_deriv_vals.at(ii) != NULL) {
            *(x_deriv_vals.at(ii)) =      cx_r   * dfields_dr.at(  ii)
                                        + cx_lon * dfields_dlon.at(ii) 
                                        + cx_lat * dfields_dlat.at(ii);
        }

        // Compute y derivatives
        if (y_deriv_vals.at(ii) != NULL) {
            *(y_deriv_vals.at(ii)) =      cy_r   * dfields_dr.at(  ii)
                                        + cy_lon * dfields_dlon.at(ii) 
                                        + cy_lat * dfields_dlat.at(ii);
        }

        // Compute z derivatives
        if (z_deriv_vals.at(ii) != NULL) {
            *(z_deriv_vals.at(ii)) =      cz_r   * dfields_dr.at(  ii)
                                        + cz_lon * dfields_dlon.at(ii) 
                                        + cz_lat * dfields_dlat.at(ii);
        }

    }
}
