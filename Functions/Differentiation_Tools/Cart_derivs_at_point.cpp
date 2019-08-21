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
    double cx_lon, cx_lat, 
           cy_lon, cy_lat, 
           cz_lon, cz_lat;
    if (constants::CARTESIAN) {
        cx_lon = 1.;
        cx_lat = 0.;

        cy_lon = 0.;
        cy_lat = 1.;

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

        // Currently assuming ddr = 0 (i.e. on a shell)

        // ddx =   ( cos(lon) * cos(lat)                ) * ddr
        //       - ( sin(lon)            / (r cos(lat)) ) * ddlon   
        //       - ( cos(lon) * sin(lat) /  r           ) * ddlat
        //cx_r   =     cos_lon  * cos_lat;
        cx_lon = -   sin_lon             / (r * cos_lat );
        cx_lat = -   cos_lon  * sin_lat  /  r;

        // ddy =   ( sin(lon) * cos(lat)                ) * ddr
        //         ( cos(lon)            / (r cos(lat)) ) * ddlon   
        //       - ( sin(lon) * sin(lat) /  r           ) * ddlat
        //cy_r   =     sin_lon  * cos_lat;
        cy_lon =     cos_lon             / (r * cos_lat );
        cy_lat = -   sin_lon  * sin_lat  /  r;

        // ddz =   (            sin(lat)                ) * ddr
        //         (            0.                      ) * ddlon
        //         (            cos(lat) /  r           ) * ddlat
        //cz_r   =                sin_lat;
        cz_lon =                0.;
        cz_lat =                cos_lat  /  r;
    }

    // Now do the actual conversion
    for (int ii = 0; ii < num_deriv; ii++) {

        // Compute x derivatives
        if (x_deriv_vals.at(ii) != NULL) {
            *(x_deriv_vals.at(ii)) = 
                  cx_lon * dfields_dlon.at(ii) 
                + cx_lat * dfields_dlat.at(ii);
        }

        // Compute y derivatives
        if (y_deriv_vals.at(ii) != NULL) {
            *(y_deriv_vals.at(ii)) = 
                  cy_lon * dfields_dlon.at(ii) 
                + cy_lat * dfields_dlat.at(ii);
        }

        // Compute z derivatives
        if (z_deriv_vals.at(ii) != NULL) {
            *(z_deriv_vals.at(ii)) = 
                  cz_lon * dfields_dlon.at(ii) 
                + cz_lat * dfields_dlat.at(ii);
        }

    }
}
