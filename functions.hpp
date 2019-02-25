#ifndef FUNCTIONS_HPP
#define FUNCTIONS_HPP 1

#include <stdio.h>
#include <stdlib.h>

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
        double * areas, 
        const double * longitude, 
        const double * latitude, 
        const int Nlon, 
        const int Nlat);

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
 */
double distance(const double lon1, const double lat1, 
                const double lon2, const double lat2);

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
void filtering(const double * u_r, const double * u_lon, const double * u_lat,
               const double * scales, const int Nscales,
               const double dlon, const double dlat,
               const int Ntime, const int Ndepth, const int Nlon, const int Nlat,
               const double * dAreas, 
               const double * time, const double * depth,
               const double * longitude, const double * latitude,
               const double * mask);

/*!
 * \brief Computed filtered Cartesian velocities at a single point
 *
 * Computes the integral of each Cartesian velocity with the
 * kernel().
 *
 * dArea for integration computed in compute_areas() 
 */
void apply_filter_at_point(
        double & u_x_tmp,   double & u_y_tmp,   double & u_z_tmp,
        const double * u_x, const double * u_y, const double * u_z,
        const int dlon_N, const int dlat_N, 
        const int Ntime,  const int Ndepth, const int Nlat, const int Nlon,
        const int Itime,  const int Idepth, const int Ilat, const int Ilon,
        const double * longitude, const double * latitude,
        const double * dAreas, const double scale,
        const double * mask);

/*!
 * \brief Primary kernel function coarse-graining procedure (G in publications)
 */
double kernel(const double distance, const double scale);

/*!
 * \brief Assign appropriate differentiation vector
 *
 * This function produces the differentiation vector necessary 
 *   for computing the value of a derivative at a single point.
 *
 * Currently assumes fourth-order finite differencing on a five-point uniform stencil.
 *
 * The grid need not be centred (to account for coastlines).
 */
void differentiation_vector(double * diff_array, const double delta, const int index);

/*!
 * \brief Compute the (spherical) vorticity at a given point.
 *
 * Uses fourth-order finite differentiation for derivatives.
 *
 * Assumes that the lon-lat grid is uniform.
 *
 * Stencil placement side-steps coast/land.
 *
 * Currently only computes the vort_r component.
 */
void compute_vorticity_at_point(
        double & vort_r_tmp, double & vort_lon_tmp, double & vort_lat_tmp,
        const double * u_r, const double * u_lon, const double * u_lat,
        const int Ntime,  const int Ndepth, const int Nlat, const int Nlon,
        const int Itime,  const int Idepth, const int Ilat, const int Ilon,
        const double * longitude, const double * latitude,
        const double * mask);

#endif
