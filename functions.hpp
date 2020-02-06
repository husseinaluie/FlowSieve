#ifndef FUNCTIONS_HPP
#define FUNCTIONS_HPP 1

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <string>
#include <map>
#include <mpi.h>
#include "constants.hpp"

/*!
 * \file
 * \brief Collection of all computation-related functions.
 */

/*! 
 * \brief Compute the area of each computational cell
 *
 * Works for both spherical and Cartesian coordinated (based on constants.hpp)
 * as well as uniform and non-uniform grids (also based on constants.hpp).
 *
 * @param[in,out]   areas       array in which areas will be stored
 * @param[in]       longitude   longitude vector (1D)
 * @param[in]       latitude    latitude vector (1D)
 *
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
 *
 * @param[in] Itime,Idepth,Ilat,Ilon    4-indices to be converted to a 1-index
 * @param[in] Ntime,Ndepth,Nlat,Nlon    dimension sizes
 *
 * @returns The effective 1-index that corresponds to the 4-index tuplet
 *
 */
int Index( const int Itime, const int Idepth, const int Ilat, const int Ilon,
           const int Ntime, const int Ndepth, const int Nlat, const int Nlon  );

/*!
 * \brief Convenience tool to convert logical index to physical index (time, depth, lat, lon).
 *
 * Index is a function to convert a one-point (logical) index
 *   into a four-point (physical) index (Itime, Idepth, Ilat, Ilon)
 *   to access the double arrays.
 *
 * Assumes standard CF ordering: time-depth-lat-lon
 *
 * @param[in]       index                   logical index to be converted to a 4-index
 * @param[in,out]   Itime,Idepth,Ilat,Ilon  4-indices to be returned
 * @param[in]       Ntime,Ndepth,Nlat,Nlon  dimension sizes
 *
 */
void Index1to4( const size_t index, 
        int & Itime, int & Idepth, int & Ilat, int & Ilon,
        const int & Ntime, const int & Ndepth, const int & Nlat, const int & Nlon  );

/*!
 * \brief Compute the distance (in metres) between two points in the domain.
 *
 * In spherical coordinates,computes the distance between two points on
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
 * @param[in]   lon1,lat1   coordinates for the first position
 * @param[in]   lon2,lat2   coordinates for the second position
 * @param[in]   Llon,Llat   physical length of the dimensions
 *
 * @returns returns the distance (in metres) between two points.
 *
 */
double distance(const double lon1,     const double lat1, 
                const double lon2,     const double lat2,
                const double Llon = 0, const double Llat = 0);

/*!
 * \brief Compute an array of the kernel values from a 
 * given reference point to every other point in the domain
 *
 * (ref_ilat, ref_ilon) is the reference point from which 
 *   the kernel values are computed.
 *
 * LAT_lb and LAT_ub are the (pre-computed) latitudinal bounds for the kernel.
 *
 * @param[in,out]   local_kernel            where to store the local kernel
 * @param[in]       scale                   Filtering scale
 * @param[in]       longitude,latitude      grid vectors (1D)
 * @param[in]       ref_ilat,ref_ilon       reference coordinate (kernel centre)
 * @param[in]       Ntime,Ndepth,Nlat,Nlon  dimension sizes
 * @param[in]       LAT_lb,LAT_ub           upper and lower latitudinal bounds for kernel
 *
 */
void compute_local_kernel(
        std::vector<double> & local_kernel,
        const double scale,
        const std::vector<double> & longitude,
        const std::vector<double> & latitude,
        const int ref_ilat, const int ref_ilon,
        const int Ntime,    const int Ndepth,
        const int Nlat,     const int Nlon,
        const int LAT_lb,   const int LAT_ub);


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
 *
 * @param[in,out]   u_x,u_y,u_z         Computed Cartesian velocities
 * @param[in]       u_r,u_lon,u_lat     Spherical velocities to be converted
 * @param[in]       lon,lat             coordinates of the location of conversion
 *
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
 *
 * @param[in,out]   u_r,u_lon,u_lat     Computed Spherical velocities
 * @param[in]       u_x,u_y,u_z         Cartesian velocities to be converted
 * @param[in]       lon,lat             coordinates of the location of conversion
 *
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
 *
 * @param[in]   u_r,u_lon,u_lat                     full velocity components
 * @param[in]   rho                                 full density
 * @param[in]   p                                   full pressure
 * @param[in]   scales                              scales at which to filter the data
 * @param[in]   dAreas                              cell areas (2D)
 * @param[in]   time, depth, longitude, latitude    dimension vectors (1D)
 * @param[in]   mask                                vector to distinguish land/water
 * @param[in]   myCounts                            Local (to MPI process) dimension sizes
 * @param[in]   myStarts                            Vector indicating where the local (to MPI process) region fits in the whole
 * @param[in]   comm                                MPI communicator (default MPI_COMM_WORLD)
 *
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
 * \brief Alternate filtering driver
 *
 * This function applies straight filtering to the input fields.
 * No secondary fields (Pi, etc.) are computed. Does NOT handle velocites
 * (i.e. does not convert to Cartesian velocities), just handles scalars.
 */
void filter_fields(
        const std::vector<const std::vector<double>*> & fields,
        const std::vector<std::string> var_names,
        const std::vector<double> & scales,
        const std::vector<double> & dAreas,
        const std::vector<double> & time,
        const std::vector<double> & depth,
        const std::vector<double> & longitude,
        const std::vector<double> & latitude,
        const std::vector<double> & mask,
        const std::vector<int>    & myCounts,
        const std::vector<int>    & myStarts,
        const MPI_Comm comm = MPI_COMM_WORLD
        );

/*!
 * \brief Alternate filtering driver for subsetting
 *
 */
void filtering_subsets(
        const std::vector<double> & u_r, 
        const std::vector<double> & u_lon, 
        const std::vector<double> & u_lat,
        const std::vector<double> & rho,
        const std::vector<double> & p,
        const std::vector<double> & scales, 
        const std::vector<double> & windows,
        const int & Nsamples,
        const std::vector<double> & dAreas, 
        const std::vector<double> & time, 
        const std::vector<double> & depth,
        const std::vector<double> & longitude, 
        const std::vector<double> & latitude,
              std::vector<double> & mask,
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
 *
 * @param[in,out]   coarse_val              where to store filtered value
 * @param[in]       fields                  fields to filter
 * @param[in]       Ntime,Ndepth,Nlat,Nlon  length of time dimension
 * @param[in]       Itime,Idepth,Ilat,Ilon  current position in time dimension
 * @param[in]       longitude,latitude      grid vectors (lon,lat)
 * @param[in]       LAT_lb,LAT_ub           lower/upper boundd on latitude for kernel
 * @param[in]       dAreas                  array of cell areas (2D - lat,lon)
 * @param[in]       scale                   filtering scale
 * @param[in]       mask                    array to distinguish land from water
 * @param[in]       use_mask                array of booleans indicating whether or not to use mask (i.e. zero out land) or to use the array value
 * @param[in]       local_kernel            pointer to pre-computed kernel (NULL indicates not provided)
 *
 */
void apply_filter_at_point(
        std::vector<double*> & coarse_val,   
        const std::vector<const std::vector<double>*> & fields,
        const int Ntime,  const int Ndepth, const int Nlat, const int Nlon,
        const int Itime,  const int Idepth, const int Ilat, const int Ilon,
        const std::vector<double> & longitude, 
        const std::vector<double> & latitude,
        const int LAT_lb,
        const int LAT_ub,
        const std::vector<double> & dAreas, 
        const double scale,
        const std::vector<double> & mask,
        const std::vector<bool> & use_mask,
        const std::vector<double> * local_kernel);

/*!
 * \brief Primary kernel function coarse-graining procedure (G in publications)
 *
 * @param[in]   distance    distance for evaluating the kernel
 * @param[in]   scale       filter scale (in metres)
 * 
 * @returns The kernel value for a given distance and filter scale
 *
 */
double kernel(const double distance, const double scale);


/*!
 * \brief Compute alpha value for kernel (for baroclinic transfer)
 *
 * @returns Returns a 'size' coefficient for the kernel. 
 *
 */
double kernel_alpha(void);

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
        const std::vector<double> & mask,
        const MPI_Comm comm = MPI_COMM_WORLD);

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
        const int LAT_lb,
        const int LAT_ub,
        const std::vector<double> & dAreas, 
        const double scale,
        const std::vector<double> & mask,
        const std::vector<double> * local_kernel);

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

void compute_Pi(
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
 *
 * \param Where to store the computed values
 * \param Full vorticity (r   component)
 * \param Full vorticity (lon component)
 * \param Full vorticity (lat component)
 * \param Coarse density field
 * \param Coarse pressure field
 * \param Length of time dimension
 * \param Length of depth dimension
 * \param Length of latitude dimension
 * \param Length of longitude dimension
 * \param Longitude dimension (1D)
 * \param Latitude dimension (1D)
 * \param Mask array (2D) to distinguish land from water
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
 * \brief Get latitude integratation bounds
 */
void get_lat_bounds(
        int & LAT_lb,
        int & LAT_ub,
        const std::vector<double> & latitude,
        const int Ilat,
        const double scale);

/*!
 * \brief Get longitude integratation bounds
 */
void get_lon_bounds(
        int & LON_lb,
        int & LON_ub,
        const std::vector<double> & longitude,
        const int Ilon,
        const double centre_lat,
        const double curr_lat,
        const double scale);

/*!
 * \brief Print summary of compile-time variables.
 *
 * This is triggered by passing --version to the executable.
 */
void print_compile_info(
        const std::vector<double> &scales);


/*!
 * \brief Class for storing internal timings.
 *
 * This class is used to wrap the interal timings into
 *   clean expressions. 
 */
class Timing_Records {
    /*! Main dictionary for storing timings
     * 
     * Keys are strings, which are human-readable labels.
     * Values are doubles, which are the amount of time that falls under the label
     */
    std::map< std::string, double  > time_records;

    public:
        //! Constructor. Simply initializes entries to time_records as zero
        Timing_Records();

        //! Zero out each value in time_records.
        void reset();

        /*! 
         * \brief Add delta to the record given by record_name: time_records[record_name] += delta
         * @param delta a double indicating the amount of time to add to the record
         * @param record_name a string indicating which record should be updating
         */
        void add_to_record(double delta, std::string record_name);

        //! Print the timing information in a human-readable format.
        void print();
};

#endif
