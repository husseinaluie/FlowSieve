#include "../constants.hpp"
#include "../functions.hpp"
#include "../differentiation_tools.hpp"
#include <algorithm>
#include <vector>
#include <omp.h>
#include <math.h>


/*!
 * \brief Compute KE transport caused by div(J)
 *
 *
 * Currently implements:
 *    - advection by coarse-scale velocity
 *    - pressure-induced transport
 *    - advection by fine-scale velocity
 *
 * NOT implemented:
 *    - Diffusion
 *
 * \verbatim
 *  J_transport =   0.5 * rho0 * | u_l |^2 * u_l
 *                + P_l * u_l
 *                - nu * 0.5 * rho * grad( | u_l |^2 )
 *                + rho0 * u_l * tau(u_l, u_l)
 *
 *
 *  (Spherical)
 *     div(J) = ( 
 *                (1 / (r * cos(lat)) ) d/dlon (J_lon),
 *                (1 / (r * cos(lat)) ) d/dlat (J_lat * cos(lat)),
 *                (1 /  r^2           ) d/dr   (J_r * r^2)
 *              )
 *
 *  (Cartesian)
 *     div(J) = ( 
 *                d/dx (J_x),
 *                d/dy (J_y),
 *                d/dz (J_z) 
 *              )
 *
 *
 *  Term 1: 0.5 * rho0 * | u_l |^2 * u_l
 *      This is advection of large-scale KE by the large-scale 
 *      velocity. 
 *
 *      (index form: 0.5 * rho0 * [ (u_i*u_i) * u_j ]    )
 *      (   of grad: 0.5 * rho0 * [ (u_i*u_i) * u_j ],j  )
 *               = 0.5 * rho0 * (u_i*u_i),j * u_j 
 *               =       rho0 * u_i * u_i,j * u_j
 *
 *
 *
 *  Term 2: P_l * u_l
 *      Transport caused by pressure 
 *
 *      (index form:  p * u_j     )
 *      (   of grad: (p * u_j),j  )
 *
 *
 *
 *  Term 3: - nu * 0.5 * rho * grad( | u_l |^2 )
 *      This is diffusion
 *      NOT YET IMPLEMENTED
 *
 *
 *
 *  Term 4: rho0 * u_l * tau(u, u)
 *      This is advection of the large-scale KE
 *      by the small-scale flow.
 *
 *      (index form:  rho0 *   u_i * tau_ij     )
 *      (   of grad:  rho0 * [ u_i * tau_ij ],j )
 *        = rho0 * ( u_i,j * tau_ij  + u_i * tau_ij,j )
 *        = rho0 * 
 *            (  
 *               u_i,j * ( bar(u_i*u_j)   - bar(u_i  )*bar(u_j) )
 *             + u_i   * ( bar(u_i*u_j),j - bar(u_i,j)*bar(u_j) )
 *            )
 *  \endverbatim
 *
 *
 * @param[in,out]   div_J                           where to store the computed values (array)
 * @param[in]       u_x,u_y,u_z                     coarse Cartesian velocity components
 * @param[in]       uxux,uxuy,uxuz,uyuy,uyuz,uzuz   coarse velocity products (e.g. bar(u*v) )  
 * @param[in]       coarse_p                        coarse pressure
 * @param[in]       longitude,latitude              1D dimension vectors
 * @param[in]       Ntime,Ndepth,Nlat,Nlon          Size of dimensions (MPI-local sizes)
 * @param[in]       mask                            2D array to distinguish land from water
 *
 */

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
        const std::vector<bool> & mask
        ) {

    const int OMP_chunksize = get_omp_chunksize(Nlat,Nlon);

    double div_J_tmp; 
    
    double dpdx, dpdy, dpdz;

    int Itime, Idepth, Ilat, Ilon;
    size_t index;
    const size_t Npts = u_x.size();

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

    // Set up the derivatives to pass through the differentiation functions
    std::vector<double*> x_deriv_vals, y_deriv_vals, z_deriv_vals;
    std::vector<const std::vector<double>*> deriv_fields;

    deriv_fields.push_back(&u_x);
    deriv_fields.push_back(&u_y);
    deriv_fields.push_back(&u_z);
    deriv_fields.push_back(&uxux);
    deriv_fields.push_back(&uxuy);
    deriv_fields.push_back(&uxuz);
    deriv_fields.push_back(&uyuy);
    deriv_fields.push_back(&uyuz);
    deriv_fields.push_back(&uzuz);
    if (constants::COMP_BC_TRANSFERS) {
        deriv_fields.push_back(&coarse_p);
    }
    
    #pragma omp parallel \
    default(none) \
    shared( div_J, stdout,\
            latitude, longitude, mask, \
            u_x, u_y, u_z, uxux, uxuy, uxuz,\
            uyuy, uyuz, uzuz,\
            deriv_fields)\
    private(Itime, Idepth, Ilat, Ilon, index, \
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
            dpdx, dpdy, dpdz,\
            x_deriv_vals, y_deriv_vals, z_deriv_vals,\
            div_J_tmp)
    {
        x_deriv_vals.push_back(&ux_x);
        x_deriv_vals.push_back(&uy_x);
        x_deriv_vals.push_back(&uz_x);
        x_deriv_vals.push_back(&uxux_x);
        x_deriv_vals.push_back(&uyux_x);
        x_deriv_vals.push_back(&uzux_x);
        x_deriv_vals.push_back(NULL);
        x_deriv_vals.push_back(NULL);
        x_deriv_vals.push_back(NULL);
        if (constants::COMP_BC_TRANSFERS) {
            x_deriv_vals.push_back(&dpdx);
        }

        y_deriv_vals.push_back(&ux_y);
        y_deriv_vals.push_back(&uy_y);
        y_deriv_vals.push_back(&uz_y);
        y_deriv_vals.push_back(NULL);
        y_deriv_vals.push_back(&uxuy_y);
        y_deriv_vals.push_back(NULL);
        y_deriv_vals.push_back(&uyuy_y);
        y_deriv_vals.push_back(&uzuy_y);
        y_deriv_vals.push_back(NULL);
        if (constants::COMP_BC_TRANSFERS) {
            y_deriv_vals.push_back(&dpdy);
        }

        z_deriv_vals.push_back(&ux_z);
        z_deriv_vals.push_back(&uy_z);
        z_deriv_vals.push_back(&uz_z);
        z_deriv_vals.push_back(NULL);
        z_deriv_vals.push_back(NULL);
        z_deriv_vals.push_back(&uxuz_z);
        z_deriv_vals.push_back(NULL);
        z_deriv_vals.push_back(&uyuz_z);
        z_deriv_vals.push_back(&uzuz_z);
        if (constants::COMP_BC_TRANSFERS) {
            z_deriv_vals.push_back(&dpdz);
        }

        #pragma omp for collapse(1) schedule(dynamic, OMP_chunksize)
        for (index = 0; index < Npts; index++) {

            div_J_tmp = constants::fill_value;

            if ( mask.at(index) ) { // Skip land areas

                Index1to4(index, Itime, Idepth, Ilat, Ilon,
                                 Ntime, Ndepth, Nlat, Nlon);

                div_J_tmp = 0.;

                // Compute the desired derivatives
                Cart_derivatives_at_point(
                        x_deriv_vals, y_deriv_vals,
                        z_deriv_vals, deriv_fields,
                        latitude, longitude,
                        Itime, Idepth, Ilat, Ilon,
                        Ntime, Ndepth, Nlat, Nlon,
                        mask);

                // u_i
                ux = u_x.at(index);
                uy = u_y.at(index);
                uz = u_z.at(index);

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

                // rho0 * ( u_i,j * ( bar(u_i*u_j) - bar(u_i)*bar(u_j) ) )
                div_J_tmp += constants::rho0 *
                    ( // j across, i down
                        ux_x * ( uxux_loc - ux*ux ) + ux_y * ( uxuy_loc - ux*uy ) + ux_z * ( uxuz_loc - ux*uz )
                      + uy_x * ( uyux_loc - uy*ux ) + uy_y * ( uyuy_loc - uy*uy ) + uy_z * ( uyuz_loc - uy*uz )
                      + uz_x * ( uzux_loc - uz*ux ) + uz_y * ( uzuy_loc - uz*uy ) + uz_z * ( uzuz_loc - uz*uz )
                    );

                // rho0 * ( u_i * ( bar(u_i*u_j),j - bar(u_i,j)*bar(u_j) ) )
                div_J_tmp += constants::rho0 *
                    ( // j across, i down
                        ux * ( (uxux_x - ux_x*ux) + (uxuy_y - ux_y*uy) + ( uxuz_z - ux_z*uz) )
                      + uy * ( (uyux_x - uy_x*ux) + (uyuy_y - uy_y*uy) + ( uyuz_z - uy_z*uz) )
                      + uz * ( (uzux_x - uz_x*ux) + (uzuy_y - uz_y*uy) + ( uzuz_z - uz_z*uz) )
                    );


                // Pressure term
                // (p * u_j),j = u_j * p_,j
                if (constants::COMP_BC_TRANSFERS) {
                    div_J_tmp += ux * dpdx + uy * dpdy + uz * dpdz;
                }

            } // end if(water) block

            div_J.at(index) = div_J_tmp;

        } // end index loop
    } // end pragma
}
