#include <math.h>
#include <vector>
#include <omp.h>
#include "../functions.hpp"
#include "../constants.hpp"
#include "../differentiation_tools.hpp"

/*!
 * \brief Compute the rotational component of the non-linear model of the baroclinic transfer term Lambda (see Lees and Aluie 2019)
 *
 * Specifically, it computes
 * \f[
 *      \Lambda_{\mathrm{rot}} = \frac{1}{2}\alpha_{\mathrm{kernel}}l^2 \frac{1}{\overline{\rho}}\left[ \frac{1}{2}\overline{\omega} \cdot \left( \nabla\overline{\rho}\times\nabla\overline{P} \right)  \right] 
 * \f]
 * where \f$ \alpha_{\mathrm{kernel}} \f$ is a multiplicative coefficient that depends on the kernel (see kernel_alpha.cpp) and \f$ l\f$ is the filter scale.
 *
 * @param[in,out]   Lambda_rot                                          Storage array for computed values
 * @param[in]       coarse_vort_r, coarse_vort_lon, coarse_vort_lat     Components of vorticity vector
 * @param[in]       coarse_rho, coarse_p                                Coarse density and pressure (respectively)
 * @param[in]       source_data                                 Class object with dataset info
 * @param[in]       scale_factor                                        Multiplicative scale factor
 */
void  compute_Lambda_rotational(
    std::vector<double> & Lambda_rot,
    const std::vector<double> & coarse_vort_r,
    const std::vector<double> & coarse_vort_lon,
    const std::vector<double> & coarse_vort_lat,
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
    double drhodlat, drhodlon, dpdlat, dpdlon, cos_lat;
    const double R2 = pow(constants::R_earth, 2.);
    int index, Itime, Idepth, Ilat, Ilon;

    std::vector<double*> lat_deriv_vals, lon_deriv_vals;
    std::vector<const std::vector<double>*> deriv_fields;

    deriv_fields.push_back(&coarse_rho);
    deriv_fields.push_back(&coarse_p);

    #pragma omp parallel \
    default(none) \
    shared(mask, latitude, longitude, source_data, \
            coarse_rho, Lambda_rot, coarse_vort_r,\
            deriv_fields)\
    private(Itime, Idepth, Ilat, Ilon, index,\
            drhodlat, dpdlat, drhodlon, dpdlon, cos_lat,\
            lat_deriv_vals, lon_deriv_vals) \
    firstprivate( Ntime, Ndepth, Nlat, Nlon, scale_factor )
    {
        lat_deriv_vals.push_back(&drhodlat);
        lat_deriv_vals.push_back(&dpdlat);

        lon_deriv_vals.push_back(&drhodlon);
        lon_deriv_vals.push_back(&dpdlon);
        #pragma omp for collapse(2) schedule(guided)
        for (Ilat = 0; Ilat < Nlat; Ilat++) {
            for (Ilon = 0; Ilon < Nlon; Ilon++) {

                cos_lat = cos(latitude.at(Ilat));

                for (Itime = 0; Itime < Ntime; Itime++) {
                    for (Idepth = 0; Idepth < Ndepth; Idepth++) {

                        index = Index(Itime, Idepth, Ilat, Ilon,
                                      Ntime, Ndepth, Nlat, Nlon);

                        if ( mask.at(index) ) { // Skip land areas

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

                            Lambda_rot.at(index) = 
                                ( coarse_rho.at(index) == 0 ) 
                                ?
                                0. //constants::fill_value
                                :
                                scale_factor * 1e4 // 1e4 is dbar to Pa
                                    * coarse_vort_r.at(index)
                                    * ( drhodlon * dpdlat  -  drhodlat * dpdlon ) 
                                    / ( 2 * coarse_rho.at(index) * R2 * cos_lat );

                        } // end if(water) block
                        else { // if(land)
                            Lambda_rot.at(index) = constants::fill_value;
                        }  // end if(land) block
                    } // end depth loop
                } // end time loop
            } // end lon loop
        } // end lat loop
    } // end pragma block
    #if DEBUG >= 2
    if (wRank == 0) { fprintf(stdout, "     ... done.\n"); }
    #endif
} // end function

