#include <math.h>
#include <vector>
#include <omp.h>
#include "../functions.hpp"
#include "../constants.hpp"
#include "../differentiation_tools.hpp"

void  compute_baroclinic_transfer(
    std::vector<double> & baroclinic_transfer,
    const std::vector<double> & coarse_vort_r,
    const std::vector<double> & coarse_vort_lon,
    const std::vector<double> & coarse_vort_lat,
    const std::vector<double> & coarse_rho,
    const std::vector<double> & coarse_p,
    const int Ntime,
    const int Ndepth,
    const int Nlat,
    const int Nlon,
    const std::vector<double> & longitude,
    const std::vector<double> & latitude,
    const std::vector<double> & mask
    ) {

    // For the moment, only use vort_r
    double drhodlat, drhodlon, dpdlat, dpdlon, cos_lat;
    const double R2 = pow(constants::R_earth, 2.);
    int index, mask_index, Itime, Idepth, Ilat, Ilon;

    std::vector<double*> lat_deriv_vals, lon_deriv_vals;
    std::vector<const std::vector<double>*> deriv_fields;

    deriv_fields.push_back(&coarse_rho);
    deriv_fields.push_back(&coarse_p);

    #pragma omp parallel \
    default(none) \
    shared(mask, latitude, longitude,\
            coarse_rho, baroclinic_transfer, coarse_vort_r,\
            deriv_fields)\
    private(Itime, Idepth, Ilat, Ilon, index, mask_index,\
            drhodlat, dpdlat, drhodlon, dpdlon, cos_lat,\
            lat_deriv_vals, lon_deriv_vals)
    {
        lat_deriv_vals.push_back(&drhodlat);
        lat_deriv_vals.push_back(&dpdlat);

        lon_deriv_vals.push_back(&drhodlon);
        lon_deriv_vals.push_back(&dpdlon);
        #pragma omp for collapse(2) schedule(guided)
        for (Ilat = 0; Ilat < Nlat; Ilat++) {
            for (Ilon = 0; Ilon < Nlon; Ilon++) {

                mask_index = Index(0,     0,      Ilat, Ilon,
                                   Ntime, Ndepth, Nlat, Nlon);
                cos_lat = cos(latitude.at(Ilat));

                for (Itime = 0; Itime < Ntime; Itime++) {
                    for (Idepth = 0; Idepth < Ndepth; Idepth++) {

                        index = Index(Itime, Idepth, Ilat, Ilon,
                                      Ntime, Ndepth, Nlat, Nlon);

                        if (mask.at(mask_index) == 1) { // Skip land areas

                            // We need a few derivatives
                            spher_derivative_at_point(
                                    lat_deriv_vals, deriv_fields,
                                    latitude, "lat",
                                    Itime, Idepth, Ilat, Ilon,
                                    Ntime, Ndepth, Nlat, Nlon,
                                    mask);

                            spher_derivative_at_point(
                                    lon_deriv_vals, deriv_fields,
                                    longitude, "lon",
                                    Itime, Idepth, Ilat, Ilon,
                                    Ntime, Ndepth, Nlat, Nlon,
                                    mask);

                            baroclinic_transfer.at(index) = 
                                coarse_vort_r.at(index)
                                    * ( drhodlon * dpdlat  -  drhodlat * dpdlon ) 
                                    / ( coarse_rho.at(index) * R2 * cos_lat );

                        } // end if(water) block
                        else { // if(land)
                            baroclinic_transfer.at(index) = constants::fill_value;
                        }  // end if(land) block
                    } // end depth loop
                } // end time loop
            } // end lon loop
        } // end lat loop
    } // end pragma block
} // end function
