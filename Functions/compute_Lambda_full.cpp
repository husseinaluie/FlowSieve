#include <math.h>
#include <vector>
#include <omp.h>
#include "../functions.hpp"
#include "../constants.hpp"
#include "../differentiation_tools.hpp"

/*!
 * \brief Compute the full baroclinic transfer term Lambda (see Lees and Aluie 2019)
 *
 * Specifically, it computes
 * \f[
 *      \Lambda_{\mathrm{rot}} = \frac{1}{\overline{\rho}} \overline{P}_{,j}\overline{\tau}(\rho,u_j) = \overline{P}_{,j}\left(\widetilde{u}-\overline{u}\right)
 * \f]
 *
 * @param[in,out]   Lambda                                      Storage array for computed values
 * @param[in]       coarse_u_r, coarse_u_lon, coarse_u_lat      Components of velocity (bar filtered)
 * @param[in]       tilde_u_r,  tilde_u_lon,  tilde_u_lat       Components of velocity (tilde filtered)
 * @param[in]       coarse_p                                    Coarse pressure 
 * @param[in]       Ntime, Ndepth, Nlat, Nlon                   Size of time, depth, lat, lon dimensions (respectively)
 * @param[in]       longitude, latitude                         Grid vectors
 * @param[in]       mask                                        Mask to distinguish land from water
 */
void  compute_Lambda_full(
    std::vector<double> & Lambda,
    const std::vector<double> & coarse_u_r,
    const std::vector<double> & coarse_u_lon,
    const std::vector<double> & coarse_u_lat,
    const std::vector<double> & tilde_u_r,
    const std::vector<double> & tilde_u_lon,
    const std::vector<double> & tilde_u_lat,
    const std::vector<double> & coarse_p,
    const int Ntime,
    const int Ndepth,
    const int Nlat,
    const int Nlon,
    const std::vector<double> & longitude,
    const std::vector<double> & latitude,
    const std::vector<bool>   & mask
    ) {

    const int OMP_chunksize = get_omp_chunksize(Nlat,Nlon);

    // For the moment, only use vort_r
    double dpdlat, dpdlon, cos_lat; 
    int Itime, Idepth, Ilat, Ilon;
    size_t index;

    const double R_earth = constants::R_earth;

    std::vector<double*> lat_deriv_vals, lon_deriv_vals;
    std::vector<const std::vector<double>*> deriv_fields;

    deriv_fields.push_back(&coarse_p);

    #pragma omp parallel \
    default(none) \
    shared(mask, latitude, longitude,\
            coarse_p, Lambda, \
            coarse_u_r, coarse_u_lon, coarse_u_lat, \
            tilde_u_r, tilde_u_lon, tilde_u_lat, \
            deriv_fields)\
    private(Itime, Idepth, Ilat, Ilon, index,\
            dpdlat, dpdlon, lat_deriv_vals, lon_deriv_vals, cos_lat)
    {
        lat_deriv_vals.push_back(&dpdlat);
        lon_deriv_vals.push_back(&dpdlon);
        #pragma omp for collapse(1) schedule(guided, OMP_chunksize)
        for (index = 0; index < Lambda.size(); ++index) {

            Index1to4(index, Itime, Idepth, Ilat, Ilon,
                             Ntime, Ndepth, Nlat, Nlon);

            if ( mask.at(index) ) { // Skip land areas

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

                cos_lat = cos(latitude.at(Ilat));

                Lambda.at(index) = 
                      ( dpdlon / ( R_earth * cos_lat) ) * ( tilde_u_lon.at(index) - coarse_u_lon.at(index) ) 
                    + ( dpdlat /   R_earth            ) * ( tilde_u_lat.at(index) - coarse_u_lat.at(index) );

            } 
            else { 
                Lambda.at(index) = constants::fill_value;
            }  
        } // end index loop
    } // end pragma block
} // end function

