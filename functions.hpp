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
 * \brief Class to store main variables.
 *
 * This class is used to hold the various source data variables
 *    such as velocity components, mask, region definitions, density, etc
 *
 */
class dataset {

    public:

        // Storage for processor assignments
        int Nprocs_in_time, Nprocs_in_depth;

        // Vectors to store the dimension variables
        std::vector<double> time, depth, latitude, longitude;
        int Ntime = -1, Ndepth = -1, Nlat = -1, Nlon = -1;
        int full_Ntime = -1, full_Ndepth = -1;

        // Store cell areas
        std::vector<double> areas;

        // Dictionary for the variables (velocity components, density, etc)
        std::map< std::string , std::vector<double> > variables;

        // Vectors and dictionary for storing the regions.
        //    these are what will be used for the post-processing
        std::vector< std::string > region_names;
        std::map< std::string, std::vector<bool> > regions;
        std::vector<double> region_areas;

        // Store mask data (i.e. land vs water)
        std::vector<bool> mask;

        // Store data-chunking info. These keep track of the MPI divisions to ensure 
        // that the output is in the same order as the input.
        std::vector<int> myCounts, myStarts;

        // Store variable descriptions to add to netcdf files
        std::map< std::string, std::string > variable_descriptions;

        // Constructor
        dataset();

        // Dimension loaders
        void load_time(      const std::string dim_name, const std::string filename );
        void load_depth(     const std::string dim_name, const std::string filename );
        void load_latitude(  const std::string dim_name, const std::string filename );
        void load_longitude( const std::string dim_name, const std::string filename );

        // Compute areas
        void compute_cell_areas();

        // Load in variable and store in dictionary
        void load_variable( const std::string var_name, 
                            const std::string var_name_in_file, 
                            const std::string filename,
                            const bool read_mask = true,
                            const bool load_counts = true );

        void build_variable_descriptions();

        // Load in region definitions
        void load_region_definitions(   const std::string filename, 
                                        const std::string dim_name, 
                                        const std::string var_name, 
                                        const MPI_Comm = MPI_COMM_WORLD );
        void compute_region_areas();

        // Check the processors divions between dimensions
        void check_processor_divisions( const int Nprocs_in_time_input, const int Nprocs_in_depth_input, const MPI_Comm = MPI_COMM_WORLD );

};

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
size_t Index( const int Itime, const int Idepth, const int Ilat, const int Ilon,
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
 * \brief Compute (local) KE density from the provided velocities.
 *
 * KE is simply computed pointwise as 0.5 * rho0 * (u1^2 + u2^2 + u3^2)
 *
 * @param[in,out]   KE          Computed KE
 * @param[in]       u1,u2,u3    Velocities
 * @param[in]       mask        differentiate land from water 
 * @param[in]       rho0        constant density
 *
 */
void KE_from_vels(
            std::vector<double> & KE,
            std::vector<double> * u1,
            std::vector<double> * u2,
            std::vector<double> * u3,
            const std::vector<bool> & mask,
            const double rho0 = constants::rho0
        );

/*!
 * \brief Wrapper that applies vel_Spher_to_Cart_at_point to every point in the domain.
 *
 * @param[in,out]   u_x,u_y,u_z                         Computed Cartesian velocities
 * @param[in]       u_r,u_lon,u_lat                     Spherical velocities to convert
 * @param[in]       mask                                differentiate land from water 
 * @param[in]       time, depth, latitude, longitude    grid vectors (1D)
 *
 */
void vel_Spher_to_Cart(
            std::vector<double> & u_x,
            std::vector<double> & u_y,
            std::vector<double> & u_z,
            const std::vector<double> & u_r,
            const std::vector<double> & u_lon,
            const std::vector<double> & u_lat,
            const std::vector<bool> & mask,
            const std::vector<double> & time,
            const std::vector<double> & depth,
            const std::vector<double> & latitude,
            const std::vector<double> & longitude
            );

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
void vel_Spher_to_Cart_at_point(
            double & u_x,
            double & u_y,
            double & u_z,
            const double u_r,
            const double u_lon,
            const double u_lat,
            const double lon,
            const double lat
        );


/*!
 * \brief Wrapper that applies vel_Spher_to_Cart_at_point to every point in the domain.
 *
 * @param[in,out]   u_r,u_lon,u_lat                     Computed Spherical velocities
 * @param[in]       u_x,u_y,u_z                         Cartesian velocities to convert
 * @param[in]       mask                                differentiate land from water 
 * @param[in]       time, depth, latitude, longitude    grid vectors (1D)
 *
 */
void vel_Cart_to_Spher(
            std::vector<double> & u_r,
            std::vector<double> & u_lon,
            std::vector<double> & u_lat,
            const std::vector<double> & u_x,
            const std::vector<double> & u_y,
            const std::vector<double> & u_z,
            const std::vector<bool> & mask,
            const std::vector<double> & time,
            const std::vector<double> & depth,
            const std::vector<double> & latitude,
            const std::vector<double> & longitude
            );


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
void vel_Cart_to_Spher_at_point(
            double & u_r, 
            double & u_lon, 
            double & u_lat,
            const double u_x, 
            const double u_y, 
            const double u_z,
            const double lon, 
            const double lat 
            );

/*!
 * \brief Main filtering driver
 *
 * This function is the main filtering driver. It sets up the appropriate
 * loop sequences, calls the other funcations (velocity conversions), and
 * calls the IO functionality.
 *
 * @param[in]   source_data     dataset class instance containing data (velocities, etc)
 * @param[in]   scales          scales at which to filter the data
 * @param[in]   comm            MPI communicator (default MPI_COMM_WORLD)
 *
 */
void filtering(const dataset & source_data,
               const std::vector<double> & scales, 
               const MPI_Comm comm = MPI_COMM_WORLD);



/*!
 * \brief Main filtering driver for Helmholtz decomposed data
 *
 * This function is the main filtering driver. It sets up the appropriate
 * loop sequences, calls the other funcations (velocity conversions), and
 * calls the IO functionality.
 *
 * @param[in]   source_data     dataset class instance containing data (Psi, Phi, etc)
 * @param[in]   scales          scales at which to filter the data
 * @param[in]   comm            MPI communicator (default MPI_COMM_WORLD)
 *
 */
void filtering_helmholtz(
        const dataset & source_data,
        const std::vector<double> & scales,
        const MPI_Comm comm = MPI_COMM_WORLD
        );


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
        const std::vector<bool> & mask,
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
              std::vector<bool> & mask,
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
 * @param[in]       weight                  pointer to spatial weight (i.e. rho) (NULL indicates not provided)
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
        const std::vector<bool> & mask,
        const std::vector<bool> & use_mask,
        const std::vector<double> * local_kernel,
        const std::vector<double> * weight = NULL
        );

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
 * Also computes the divergence and OkuboWeiss parameter (since we're at it anyways)
 */
void compute_vorticity_at_point(
        double & vort_r_tmp, 
        double & vort_lon_tmp, 
        double & vort_lat_tmp,
        double & div_tmp,
        double & OkuboWeiss_tmp,
        const std::vector<double> & u_r, 
        const std::vector<double> & u_lon, 
        const std::vector<double> & u_lat,
        const int Ntime,  const int Ndepth, const int Nlat, const int Nlon,
        const int Itime,  const int Idepth, const int Ilat, const int Ilon,
        const std::vector<double> & longitude, 
        const std::vector<double> & latitude,
        const std::vector<bool> & mask);

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
        std::vector<double> & vel_div,
        std::vector<double> & OkuboWeiss,
        const std::vector<double> & u_r, 
        const std::vector<double> & u_lon, 
        const std::vector<double> & u_lat,
        const int Ntime, const int Ndepth, const int Nlat, const int Nlon,
        const std::vector<double> & longitude, 
        const std::vector<double> & latitude,
        const std::vector<bool> & mask,
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
        const std::vector<bool> & mask,
        const std::vector<double> * local_kernel);

/*!
 * \brief Compute the energy transfer through the current filter scale
 */
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
        const std::vector<bool> & mask);

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
 * @param[in]       Ntime, Ndepth, Nlat, Nlon                           Size of time, depth, lat, lon dimensions (respectively)
 * @param[in]       longitude, latitude                                 Grid vectors
 * @param[in]       mask                                                Mask to distinguish land from water
 * @param[in]       scale_factor                                        Multiplicative scale factor
 */
void compute_Lambda_rotational(
    std::vector<double> & Lambda_rot,
    const std::vector<double> & coarse_vort_r,
    const std::vector<double> & coarse_vort_lon,
    const std::vector<double> & coarse_vort_lat,
    const std::vector<double> & coarse_rho,
    const std::vector<double> & coarse_p,
    const int Ntime, const int Ndepth, const int Nlat, const int Nlon,
    const std::vector<double> & longitude,
    const std::vector<double> & latitude,
    const std::vector<bool> & mask,
    const double scale_factor
    );


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
    const std::vector<bool> & mask
    );

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
 * @param[in]       Ntime, Ndepth, Nlat, Nlon                   Size of time, depth, lat, lon dimensions (respectively)
 * @param[in]       longitude, latitude                         Grid vectors
 * @param[in]       mask                                        Mask to distinguish land from water
 * @param[in]       scale_factor                                Multiplicative scale factor
 */
void compute_Lambda_nonlin_model(
    std::vector<double> & Lambda_nonlin,
    const std::vector<double> & coarse_u_r,
    const std::vector<double> & coarse_u_lon,
    const std::vector<double> & coarse_u_lat,
    const std::vector<double> & coarse_rho,
    const std::vector<double> & coarse_p,
    const int Ntime, const int Ndepth, const int Nlat, const int Nlon,
    const std::vector<double> & longitude,
    const std::vector<double> & latitude,
    const std::vector<bool> & mask,
    const double scale_factor
    );


/*!
 *
 * \param Where to store the computed values
 * \param Coarse vorticity (r   component)
 * \param Coarse vorticity (lon component)
 * \param Coarse vorticity (lat component)
 * \param Coarse velocity (r   component)
 * \param Coarse velocity (lon component)
 * \param Coarse velocity (lat component)
 * \param Length of time dimension
 * \param Length of depth dimension
 * \param Length of latitude dimension
 * \param Length of longitude dimension
 * \param Longitude dimension (1D)
 * \param Latitude dimension (1D)
 * \param Mask to distinguish land from water
 */
void compute_vort_stretch(
    std::vector<double> & baroclinic_transfer,
    const std::vector<double> & full_vort_r,
    const std::vector<double> & full_vort_lon,
    const std::vector<double> & full_vort_lat,
    const std::vector<double> & coarse_u_r,
    const std::vector<double> & coarse_u_lon,
    const std::vector<double> & coarse_u_lat,
    const int Ntime, const int Ndepth, const int Nlat, const int Nlon,
    const std::vector<double> & longitude,
    const std::vector<double> & latitude,
    const std::vector<bool> & mask);

/*!
 *  \brief Interpolate the given field over the land cells
 */
void interpolate_over_land(
        std::vector<double> &field,
        const std::vector<double> &time,
        const std::vector<double> &depth,
        const std::vector<double> &latitude,
        const std::vector<double> &longitude,
        const std::vector<bool> &mask);

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
        const std::vector<bool> & mask);

/*!
 * \brief Given a field, compute the horizontal average.
 *
 * The result, means, is a function of time and depth.
 */

void compute_spatial_average(
        std::vector<double> & means,
        const std::vector<double> & field,
        const std::vector<double> & areas,
        const int Ntime,
        const int Ndepth,
        const int Nlat,
        const int Nlon,
        const std::vector<bool> & mask);

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
 *
 * @param[in]   scales      Scales used for filtering
 */
void print_compile_info(
        const std::vector<double> * scales = NULL);

/*! 
 * \brief Compute openmp chunk size for OMP for loops
 *
 * The point of using a function for this is that it standardizes
 *    the chunk size across the different functions.
 *
 * This assumes that the loop being split is a Lat-Lon loop.
 *
 */
int get_omp_chunksize(const int Nlat, const int Nlon);


/*!
 * \brief Apply degree->radian conversion, if necessary.
 *
 * This is just to provide a clean wrapper to improve some case file readability.
 * It applies a Cartesian check (via constants.hpp), so it is safe to include by default.
 *
 * @param[in]   longitude, latitude     Coordinate grids to be converted, if appropriate
 *
 */
void convert_coordinates(
        std::vector<double> & longitude,
        std::vector<double> & latitude
        );

/*!
 *
 * \brief If the grid includes the poles (lat = 90 degrees), then mask it out
 *
 * @param[in]       latitude                    Latitude grid
 * @param[in,out]   mask                        Mask array to differentiate land/water
 * @param[in]       Ntime,Ndepth,Nlat,Nlon      Dimension sizes
 *
 */
void mask_out_pole(
        const std::vector<double> & latitude,
        std::vector<bool> & mask,
        const int Ntime,
        const int Ndepth,
        const int Nlat,
        const int Nlon
        );


/*!
 *
 * \brief Roll field along dimension
 *
 * Currently hard-coded to roll along lon dimension only.
 *
 * @param[in,out]   field_to_roll
 * @param[in]       dimension
 * @param[in]       roll_count
 * @param[in]       Ntime,Ndepth,Nlat,Nlon
 *
 */
void roll_field(
        std::vector<double> & field_to_roll,
        const std::string dimension,
        const int roll_count,
        const int Ntime,
        const int Ndepth,
        const int Nlat,
        const int Nlon
        );

/*!
 * \brief Class for storing internal timings.
 *
 * This class is used to wrap the interal timings into
 *   clean expressions. 
 */
class Timing_Records {

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

    private:
        /*! Main dictionary for storing timings
         * 
         * Keys are strings, which are human-readable labels.
         * Values are doubles, which are the amount of time that falls under the label
         */
        std::map< std::string, double  > time_records;
};


/*!
 * \brief Class to process command-line arguments
 *
 * This parser tool was provided by StackOverflow user 'iain' at the following post.
 *    mild modifications have been made to adapt the code to our purposes.
 *
 * https://stackoverflow.com/questions/865668/how-to-parse-command-line-arguments-in-c
 *
 */
class InputParser {

    public:
        InputParser (int &argc, char **argv);

        const std::string getCmdOption(
                const std::string &option,
                const std::string &default_value = ""
                ) const;

        bool cmdOptionExists(const std::string &option) const;

        void getFilterScales( std::vector<double> &filter_scales, const std::string &argname ) const;

    private:
        std::vector <std::string> tokens;

};

/*!
 * \brief Convert string to boolean
 */
bool string_to_bool( std::string str );

void print_header_info( );

#endif
