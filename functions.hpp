#ifndef FUNCTIONS_HPP
#define FUNCTIONS_HPP 1

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <mpi.h>
#include "constants.hpp"

/*!
 * \file
 * \brief Collection of all computation-related functions.
 */

/*! 
 * \brief Compute the area of each computational cell
 *
 * Currently assumes spherical coordinates.
 */
void compute_areas(
        std::vector<double> & areas, 
        const std::vector<double> & longitude, 
        const std::vector<double> & latitude);

/*!
 * \brief Convenience tool to convert physical index (time, depth, lat, lon) to a logical index.
 *
 * Index is a function to convert a four-point (physical) index
 *   (Itime, Idepth, Ilat, Ilon) into a one-point (logical) index
 *   to access the double arrays.
 *
 * Assumes standard CF ordering: time-depth-lat-lon
 */
int Index( const int Itime, const int Idepth, const int Ilat, const int Ilon,
           const int Ntime, const int Ndepth, const int Nlat, const int Nlon  );

/*!
 * \brief Compute the distance (in metres) between two points on a sphere.
 *
 * This function computes the distance between two points on
 *   a spherical shell along a great circle.
 *   (see https://en.wikipedia.org/wiki/Great-circle_distance)
 * It should avoid floating point issues for the grid scales that
 *   we're considering. Worst case, increase to quads for better
 *   accuracy.
 *
 * The last two arguments (Llon and Llat) give the physical 
 *   length of the two dimensions. This is used in the case
 *   of periodic Cartesian grids. They are otherwise unused.
 *
 */
double distance(const double lon1,     const double lat1, 
                const double lon2,     const double lat2,
                const double Llon = 0, const double Llat = 0);

/*!
 * \brief Compute an array of the distances (in metres) from a 
 * given reference point to every other point in the domain
 *
 * (ref_ilat, ref_ilon) is the reference point from which 
 *   distances are computed.
 */
void compute_distances(
        std::vector<double> & distances,
        const std::vector<double> & longitude,
        const std::vector<double> & latitude,
        const int ref_ilat, const int ref_ilon,
        const int Ntime, const int Ndepth,
        const int Nlat,  const int Nlon);

/*!
 * \brief Convert single spherical velocity to Cartesian velocity
 *
 * Convert Spherical velocities to Cartesian
 *   velocities.
 * (u_r, u_lon, u_lat) -> (u_x, u_y, u_z)
 *
 * Note: we are using linear velocities (m/s),
 *   not angular (rad/s), so 
 *   \f{eqnarray*}{
 *   u_\lambda = r\cos(\phi)\cdot\hat{u}_\lambda \\
 *   u_\phi    = r\cdot\hat{u}_\phi
 *   \f}
 *
 * Note: we are still using a Spherical
 *   coordinate system, we are only converting
 *   the velocity fields.
 */
void vel_Spher_to_Cart(
            double & u_x, double & u_y, double & u_z,
            const double u_r, const double u_lon, const double u_lat,
            const double lon, const double lat );

/*!
 * \brief Convert single Cartesian velocity to spherical velocity
 *
 * Convert Cartesian velocities to spherical
 *   velocities.
 *   (u_x, u_y, u_z) -> (u_r, u_lon, u_lat)
 *
 * Note: we are using linear velocities (m/s),
 *   not angular (rad/s), so 
 *   \f{eqnarray*}{
 *   u_\lambda = r\cos(\phi)\cdot\hat{u}_\lambda \\
 *   u_\phi    = r\cdot\hat{u}_\phi
 *   \f}
 *
 * Note: we are still using a Spherical
 *   coordinate system, we are only converting
 *   the velocity fields.
 */
void vel_Cart_to_Spher(
            double & u_r, double & u_lon, double & u_lat,
            const double u_x, const double u_y, const double u_z,
            const double lon, const double lat );

/*!
 * \brief Main filtering driver
 *
 * This function is the main filtering driver. It sets up the appropriate
 * loop sequences, calls the other funcations (velocity conversions), and
 * calls the IO functionality.
 */
void filtering(const std::vector<double> & u_r, 
               const std::vector<double> & u_lon, 
               const std::vector<double> & u_lat,
               const std::vector<double> & rho,
               const std::vector<double> & p,
               const std::vector<double> & scales, 
               const std::vector<double> & dAreas, 
               const std::vector<double> & time, 
               const std::vector<double> & depth,
               const std::vector<double> & longitude, 
               const std::vector<double> & latitude,
               const std::vector<double> & mask,
               const std::vector<int>    & myCounts,
               const std::vector<int>    & myStarts,
               const MPI_Comm comm = MPI_COMM_WORLD);

/*!
 * \brief Compute filtered field at a single point
 *
 * Computes the integral of the provided field with the
 * kernel().
 *
 * dArea for integration computed in compute_areas() 
 */
void apply_filter_at_point(
        double & coarse_val,   
        const std::vector<double> & field, 
        const int Ntime,  const int Ndepth, const int Nlat, const int Nlon,
        const int Itime,  const int Idepth, const int Ilat, const int Ilon,
        const std::vector<double> & longitude, 
        const std::vector<double> & latitude,
        const std::vector<double> & dAreas, 
        const double scale,
        const std::vector<double> & mask,
        const bool use_mask,
        const std::vector<double> * distances);

/*!
 * \brief Primary kernel function coarse-graining procedure (G in publications)
 */
double kernel(const double distance, const double scale);

/*!
 * \brief Compute the vorticity at a given point.
 *
 * Assumes that the lon-lat grid is uniform.
 *
 * Currently only computes the vort_r component.
 */
void compute_vorticity_at_point(
        double & vort_r_tmp, double & vort_lon_tmp, double & vort_lat_tmp,
        const std::vector<double> & u_r, 
        const std::vector<double> & u_lon, 
        const std::vector<double> & u_lat,
        const int Ntime,  const int Ndepth, const int Nlat, const int Nlon,
        const int Itime,  const int Idepth, const int Ilat, const int Ilon,
        const std::vector<double> & longitude, 
        const std::vector<double> & latitude,
        const std::vector<double> & mask);

/*!
 * \brief Wrapper for computing vorticity
 *
 * This wrapper applies compute_vorticity_at_point() at each point in the
 * grid. If the compiler flag COMP_VORT is set to false, then this procedure
 * is skipped.
 */
void compute_vorticity(
        std::vector<double> & vort_r,    
        std::vector<double> & vort_lon,    
        std::vector<double> & vort_lat,
        const std::vector<double> & u_r, 
        const std::vector<double> & u_lon, 
        const std::vector<double> & u_lat,
        const int Ntime, const int Ndepth, const int Nlat, const int Nlon,
        const std::vector<double> & longitude, 
        const std::vector<double> & latitude,
        const std::vector<double> & mask);

/*!
 * \brief Compute filtered quadratic velocities at a single point
 *
 * Computes the integral of each quadratic Cartesian velocity
 * with the kernel().
 *
 * In particular, the quadratic terms being filtered are:
 * \$u_xu_x\$, \$u_xu_y\$, \$u_xu_z\$, \$u_yu_y\$, \$u_yu_z\$, \$u_zu_z\$
 *
 * dArea for integration computed in compute_areas() 
 */
void apply_filter_at_point_for_quadratics(
        double & uxux_tmp,   double & uxuy_tmp,   double & uxuz_tmp,
        double & uyuy_tmp,   double & uyuz_tmp,   double & uzuz_tmp,
        const std::vector<double> & u_x, 
        const std::vector<double> & u_y, 
        const std::vector<double> & u_z,
        const int Ntime,  const int Ndepth, const int Nlat, const int Nlon,
        const int Itime,  const int Idepth, const int Ilat, const int Ilon,
        const std::vector<double> & longitude, 
        const std::vector<double> & latitude,
        const std::vector<double> & dAreas, 
        const double scale,
        const std::vector<double> & mask,
        const std::vector<double> * distances);

/*!
 * \brief Compute the energy transfer through the current filter scale
 */
void compute_energy_transfer_through_scale(
        std::vector<double> & energy_transfer,
        const std::vector<double> & ux,   
        const std::vector<double> & uy,   
        const std::vector<double> & uz,
        const std::vector<double> & uxux, 
        const std::vector<double> & uxuy, 
        const std::vector<double> & uxuz,
        const std::vector<double> & uyuy, 
        const std::vector<double> & uyuz, 
        const std::vector<double> & uzuz,
        const int Ntime, const int Ndepth, const int Nlat, const int Nlon,
        const std::vector<double> & longitude, 
        const std::vector<double> & latitude,
        const std::vector<double> & mask);

/*!
 * \brief Compute the large-scale strain tensor S
 *
 * Uses Cartesian derivatives of Cartesian velocities on a spherical coordinate system.
 *
 * \f[ S = \frac{1}{2}\left( \nabla\overline{u}_l + \nabla\overline{u}_l^T \right) \f]
 */
void compute_largescale_strain(
        double & S_xx, double & S_xy, double & S_xz,
        double & S_yy, double & S_yz, double & S_zz,
        const std::vector<double> & u_x, 
        const std::vector<double> & u_y, 
        const std::vector<double> & u_z,
        const int Itime, const int Idepth, const int Ilat, const int Ilon,
        const int Ntime, const int Ndepth, const int Nlat, const int Nlon,
        const std::vector<double> & longitude, 
        const std::vector<double> & latitude, 
        const std::vector<double> & mask);

/*!
 * \brief Compute the baroclinic energy transfer through the current filter scale
 */
void compute_baroclinic_transfer(
    std::vector<double> & baroclinic_transfer,
    const std::vector<double> & full_vort_r,
    const std::vector<double> & full_vort_lon,
    const std::vector<double> & full_vort_lat,
    const std::vector<double> & coarse_rho,
    const std::vector<double> & coarse_p,
    const int Ntime, const int Ndepth, const int Nlat, const int Nlon,
    const std::vector<double> & longitude,
    const std::vector<double> & latitude,
    const std::vector<double> & mask);

/*!
 *  \brief Interpolate the given field over the land cells
 */
void interpolate_over_land(
        std::vector<double> &field,
        const std::vector<double> &time,
        const std::vector<double> &depth,
        const std::vector<double> &latitude,
        const std::vector<double> &longitude,
        const std::vector<double> &mask);

/*!
 *  \brief Convert potential temperature to actual temperature
 */
double depotential_temperature( 
        const double p, 
        const double theta);

/*!
 * \brief UNESCO equation of state
 */
double equation_of_state(
        const double T,
        const double S,
        const double p);

/*!
 * \brief Compute KE transport caused by div(J)
 *
 * Currently implements:
 *    - advection by coarse-scale velocity
 *    - pressure-induced transport
 *    - advection by fine-scale velocity
 *
 * NOT implemented:
 *    - Diffusion
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
        const std::vector<double> & mask);

/*!
 * \brief Given a field, compute the horizontal average.
 *
 * The result, means, is a function of time and depth.
 */
void compute_mean(
        std::vector<double> & means,
        const std::vector<double> & field,
        const std::vector<double> & dArea,
        const int Ntime,
        const int Ndepth,
        const int Nlat,
        const int Nlon,
        const std::vector<double> & mask);

/*!
 * \brief Compute the divergence of the filtered velocity fields.
 */
void compute_div_vel(
        std::vector<double> & div,
        const std::vector<double> & u_x,
        const std::vector<double> & u_y,
        const std::vector<double> & u_z,
        const std::vector<double> & longitude,
        const std::vector<double> & latitude,
        const int Ntime,
        const int Ndepth,
        const int Nlat,
        const int Nlon,
        const std::vector<double> & mask);

/*!
 * \brief Print summary of compile-time variables.
 *
 * This is triggered by passing --version to the executable.
 */
void print_compile_info(
        const std::vector<double> &scales);

#endif
