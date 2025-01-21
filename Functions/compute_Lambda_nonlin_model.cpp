#include <math.h>
#include <vector>
#include <omp.h>
#include "../functions.hpp"
#include "../constants.hpp"
#include "../differentiation_tools.hpp"

/*!
 * \brief Compute the non-linear model of the baroclinic transfer term Lambda (see Lees and Aluie 2019)
 *
 * Specifically, it computes
 * \f[
 *      \Lambda_{\mathrm{rot}} = \frac{1}{2}\alpha_{\mathrm{kernel}}l^2 \frac{1}{\overline{\rho}} \overline{P}_{,j}\overline{\rho}_{,k}\overline{u}_{j,k}
 *                             = \frac{1}{2}\alpha_{\mathrm{kernel}}l^2 \frac{1}{\overline{\rho}} \nabla\overline{P}\cdot \nabla\overline{\vec{u}} \cdot \nabla \overline{\rho}
 * \f]
 * where \f$ \alpha_{\mathrm{kernel}} \f$ is a multiplicative coefficient that depends on the kernel (see kernel_alpha.cpp) and \f$ l\f$ is the filter scale.
 *
 * @param[in,out]   Lambda_rot                                  Storage array for computed values
 * @param[in]       coarse_u_r, coarse_u_lon, coarse_u_lat      Components of vorticity vector
 * @param[in]       coarse_rho, coarse_p                        Coarse density and pressure (respectively)
 * @param[in]       source_data                                 Class object with dataset info
 * @param[in]       scale_factor                                Multiplicative scale factor
 */
void  compute_Lambda_nonlin_model(
    std::vector<double> & Lambda_nonlin,
    const std::vector<double> & coarse_u_r,
    const std::vector<double> & coarse_u_lon,
    const std::vector<double> & coarse_u_lat,
    const std::vector<double> & coarse_rho,
    const std::vector<double> & coarse_p,
    const dataset & source_data,
    const double scale_factor,
    const MPI_Comm comm
    ) {

    const std::vector<short int> &mask = source_data.mask;

    #if DEBUG >= 2
    int wRank, wSize;
    MPI_Comm_rank( comm, &wRank );
    MPI_Comm_size( comm, &wSize );

    if (wRank == 0) { fprintf(stdout, "  Starting Lambda_nonlin_model computation.\n"); }
    #endif

    const int   Ntime   = source_data.Ntime,
                Ndepth  = source_data.Ndepth,
                Nlat    = source_data.Nlat,
                Nlon    = source_data.Nlon;

    const std::vector<double>   &latitude   = source_data.latitude,
                                &longitude  = source_data.longitude;

    // For the moment, only use vort_r
    double drho_dlat, drho_dlon, dp_dlat, dp_dlon, dulon_dlon, dulon_dlat, dulat_dlon, dulat_dlat,
           lon_factor, lat_factor;
    int Itime, Idepth, Ilat, Ilon;
    size_t index;

    std::vector<double*> lat_deriv_vals, lon_deriv_vals;
    std::vector<const std::vector<double>*> deriv_fields;

    deriv_fields.push_back(&coarse_rho);
    deriv_fields.push_back(&coarse_p);
    deriv_fields.push_back(&coarse_u_lon);
    deriv_fields.push_back(&coarse_u_lat);

    #pragma omp parallel \
    default(none) \
    shared(mask, latitude, longitude, source_data, \
            coarse_rho, Lambda_nonlin, coarse_u_lon, coarse_u_lat,\
            deriv_fields)\
    private(Itime, Idepth, Ilat, Ilon, index,\
            lon_factor, lat_factor, \
            drho_dlat, dp_dlat, dulon_dlat, dulat_dlat, \
            drho_dlon, dp_dlon, dulon_dlon, dulat_dlon, \
            lat_deriv_vals, lon_deriv_vals) \
    firstprivate( Ntime, Ndepth, Nlat, Nlon, scale_factor )
    {
        lat_deriv_vals.push_back(&drho_dlat);
        lat_deriv_vals.push_back(&dp_dlat);
        lat_deriv_vals.push_back(&dulon_dlat);
        lat_deriv_vals.push_back(&dulat_dlat);

        lon_deriv_vals.push_back(&drho_dlon);
        lon_deriv_vals.push_back(&dp_dlon);
        lon_deriv_vals.push_back(&dulon_dlon);
        lon_deriv_vals.push_back(&dulat_dlon);
        #pragma omp for collapse(1) schedule(guided)
        for (index = 0; index < Lambda_nonlin.size(); ++index) {

            if ( mask.at(index) ) { // Skip land areas
                Index1to4(index, Itime, Idepth, Ilat, Ilon,
                                 Ntime, Ndepth, Nlat, Nlon);

                // Curvature terms for derivatives
                if (constants::CARTESIAN) {
                    lon_factor = 1.;
                    lat_factor = 1.;
                } else {
                    lon_factor = constants::R_earth * cos(latitude.at(Ilat));
                    lat_factor = constants::R_earth;
                }


                // We need a few derivatives
                spher_derivative_at_point(
                        lat_deriv_vals, deriv_fields,
                        latitude, "lat",
                        source_data, 
                        Itime, Idepth, Ilat, Ilon,
                        mask);

                spher_derivative_at_point(
                        lon_deriv_vals, deriv_fields,
                        longitude, "lon",
                        source_data, 
                        Itime, Idepth, Ilat, Ilon,
                        mask);

                Lambda_nonlin.at(index) = 
                    ( coarse_rho.at(index) == 0 ) 
                    ?
                    0. //constants::fill_value
                    :
                    scale_factor * 1e4  // 1e4 is dbar to Pa conversion
                        * (    dp_dlon * drho_dlon * dulon_dlon / ( lon_factor * lon_factor * lon_factor )
                            +  dp_dlon * drho_dlat * dulon_dlat / ( lon_factor * lat_factor * lat_factor )
                            +  dp_dlat * drho_dlon * dulat_dlon / ( lat_factor * lon_factor * lon_factor )
                            +  dp_dlat * drho_dlat * dulat_dlat / ( lat_factor * lat_factor * lat_factor )
                          )
                        / coarse_rho.at(index);

            } // end if(water) block
            else { // if(land)
                Lambda_nonlin.at(index) = constants::fill_value;
            }  // end if(land) block
        } // end index loop
    } // end pragma block
    #if DEBUG >= 2
    if (wRank == 0) { fprintf(stdout, "     ... done.\n"); }
    #endif
} // end function

