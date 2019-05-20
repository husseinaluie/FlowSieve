#include <vector>
#include <assert.h>
#include "../../differentiation_tools.hpp"
#include "../../constants.hpp"
#include "../../functions.hpp"
#include <math.h>
#include <string>

void Cart_derivatives_at_point(
        const std::vector<double*> & x_deriv_vals,
        const std::vector<double*> & y_deriv_vals,
        const std::vector<double*> & z_deriv_vals,
        const std::vector<const std::vector<double>*> & fields,
        const std::vector<double> & latitude,
        const std::vector<double> & longitude,
        const int Itime,
        const int Idepth,
        const int Ilat,
        const int Ilon,
        const int Ntime,
        const int Ndepth,
        const int Nlat,
        const int Nlon,
        const std::vector<double> & mask
        ) {

    // Confirm that input sizes match
    assert(x_deriv_vals.size() == fields.size());
    assert(y_deriv_vals.size() == fields.size());
    assert(z_deriv_vals.size() == fields.size());
    const int num_deriv = x_deriv_vals.size();

    // Currently assuming ddr = 0
    // ddx = - ( sin(lon)            / (r cos(lat)) ) * ddlon   
    //       - ( cos(lon) * sin(lat) /  r           ) * ddlat
    // ddy =   ( cos(lon)            / (r cos(lat)) ) * ddlon   
    //       - ( sin(lon) * sin(lat) /  r           ) * ddlat
    // ddz =                                          
    //         (            cos(lat) /  r           ) * ddlat
    std::vector<double*> dfields_dlon_p(num_deriv);
    std::vector<double>  dfields_dlon(  num_deriv);

    std::vector<double*> dfields_dlat_p(num_deriv);
    std::vector<double>  dfields_dlat(  num_deriv);

    for (int ii = 0; ii < num_deriv; ii++) {
        dfields_dlon_p.at(ii) = &(dfields_dlon.at(ii));
        dfields_dlat_p.at(ii) = &(dfields_dlat.at(ii));
    }

    // Compute spherical derivatives
    spher_derivative_at_point(
        dfields_dlon_p, fields, longitude, "lon",
        Itime, Idepth, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon,
        mask);

    spher_derivative_at_point(
        dfields_dlat_p, fields, latitude, "lat",
        Itime, Idepth, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon,
        mask);

    // Define conversion coefficients
    double c1x, c2x, c1y, c2y, c1z, c2z;
    if (constants::CARTESIAN) {
        c1x = 1.;
        c2x = 0.;
        c1y = 0.;
        c2y = 1.;
        c1z = 0.;
        c2z = 0.;
    } else {
        double lon = longitude.at(Ilon);
        double lat = latitude.at(Ilat);
        double r = constants::R_earth;

        double cos_lat = cos(lat);
        double cos_lon = cos(lon);
        double sin_lat = sin(lat);
        double sin_lon = sin(lon);

        c1x = - sin_lon           / (r * cos_lat );
        c2x = - cos_lon * sin_lat /  r;
        c1y =   cos_lon           / (r * cos_lat );
        c2y = - sin_lon * sin_lat /  r;
        c1z =             cos_lat /  r;
        c2z =   0.;
    }

    // Now do the actual conversion
    for (int ii = 0; ii < num_deriv; ii++) {
        if (x_deriv_vals.at(ii) != NULL) {
            *(x_deriv_vals.at(ii)) = 
                  c1x * dfields_dlon.at(ii) 
                + c2x * dfields_dlat.at(ii);
        }
        if (y_deriv_vals.at(ii) != NULL) {
            *(y_deriv_vals.at(ii)) = 
                  c1y * dfields_dlon.at(ii) 
                + c2y * dfields_dlat.at(ii);
        }
        if (z_deriv_vals.at(ii) != NULL) {
            *(z_deriv_vals.at(ii)) = 
                  c1z * dfields_dlon.at(ii) 
                + c2z * dfields_dlat.at(ii);
        }
    }
}
