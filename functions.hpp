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

        //
        //// Variables
        //

        // Storage for processor assignments
        int Nprocs_in_time, Nprocs_in_depth;

        // Vectors to store the dimension variables
        std::vector<double> time, depth, latitude, longitude;
        int Ntime = -1, Ndepth = -1, Nlat = -1, Nlon = -1;
        int full_Ntime = -1, full_Ndepth = -1;

        // MPI Communicator Objects
        MPI_Comm MPI_Comm_Global = MPI_COMM_WORLD;
        MPI_Comm MPI_subcomm_sametimes;
        MPI_Comm MPI_subcomm_samedepths;

        // Store cell areas
        std::vector<double> areas;

        // Boolean to indicate if we're computing u_r from incompressibility
        bool compute_radial_vel    = false,
             use_depth_derivatives = false,
             depth_is_elevation    = false,
             depth_is_increasing   = true;

        // Dictionary for the variables (velocity components, density, etc)
        std::map< std::string , std::vector<double> > variables;

        // Vectors and dictionary for storing the regions.
        //    these are what will be used for the post-processing
        std::vector< std::string > region_names;
        std::map< std::string, std::vector<bool> > regions;
        std::vector<double> region_areas, region_areas_water_only;

        // The vectors for the coarse lat/lon maps
        std::vector<double> coarse_map_lat, coarse_map_lon, coarse_map_areas;

        // Store mask data (i.e. land vs water)
        std::vector<bool> mask, reference_mask, mask_DEPTH;

        // Store data-chunking info. These keep track of the MPI divisions to ensure 
        // that the output is in the same order as the input.
        std::vector<int> myCounts, myStarts;

        //
        //// Functions
        //

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
                            const bool load_counts = true,
                            const bool do_splits = true );

        // Load in region definitions
        void load_region_definitions(   const std::string filename, 
                                        const std::string dim_name, 
                                        const std::string var_name, 
                                        const MPI_Comm = MPI_COMM_WORLD );
        void compute_region_areas();

        // Prepare for grid-coarsening outputs
        void prepare_for_coarsened_grids(   const std::string filename,
                                            const MPI_Comm = MPI_COMM_WORLD );

        // Check the processors divions between dimensions
        void check_processor_divisions( const int Nprocs_in_time_input, 
                                        const int Nprocs_in_depth_input, 
                                        const MPI_Comm = MPI_COMM_WORLD );

        // Function to gather a variable across all depths (i.e. reconstruct depth profile)
        //  this is necessary for things like depth derivatives
        void gather_variable_across_depth( const std::vector<double> & var,
                                            std::vector<double> & gathered_var ) const ;
        void gather_mask_across_depth( const std::vector<bool> & var,
                                                std::vector<bool> & gathered_var
                                              ) const ;

        // Indexing functions
        size_t local_index( const int Itime, 
                            const int Idepth, 
                            const int Ilat, 
                            const int Ilon 
                          ) const;
        size_t global_index(    const int Itime, 
                                const int Idepth, 
                                const int Ilat, 
                                const int Ilon, 
                                const std::string merge_kind 
                           ) const;
        size_t index_local_to_global( const size_t index, const std::string merge_kind ) const;
        size_t index_global_to_local( const size_t index, const std::string merge_kind ) const;
        void index1to4_local(   const size_t index, 
                                int & Itime, int & Idepth, int & Ilat, int & Ilon
                            ) const;
        void index1to4_global(  const size_t index, 
                                int & Itime, int & Idepth, int & Ilat, int & Ilon,
                                const std::string merge_kind
                             ) const;


};

void compute_areas(
        std::vector<double> & areas, 
        const std::vector<double> & longitude, 
        const std::vector<double> & latitude);

size_t Index( const int Itime, const int Idepth, const int Ilat, const int Ilon,
              const int Ntime, const int Ndepth, const int Nlat, const int Nlon  );

void Index1to4( const size_t index, 
        int & Itime, int & Idepth, int & Ilat, int & Ilon,
        const int & Ntime, const int & Ndepth, const int & Nlat, const int & Nlon  );

double distance(const double lon1,     const double lat1, 
                const double lon2,     const double lat2,
                const double Llon = 0, const double Llat = 0);

void compute_local_kernel(
        std::vector<double> & local_kernel,
        std::vector<double> & local_dl_kernel,
        std::vector<double> & local_dl2_kernel,
        const double scale,
        const dataset & source_data,
        const int Ilat,     const int Ilon,
        const int LAT_lb,   const int LAT_ub);

void KE_from_vels(
            std::vector<double> & KE,
            std::vector<double> * u1,
            std::vector<double> * u2,
            std::vector<double> * u3,
            const std::vector<bool> & mask,
            const double rho0 = constants::rho0
        );

void vel_Spher_to_Cart(
            std::vector<double> & u_x,
            std::vector<double> & u_y,
            std::vector<double> & u_z,
            const std::vector<double> & u_r,
            const std::vector<double> & u_lon,
            const std::vector<double> & u_lat,
            const dataset & source_data
            );

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


void vel_Cart_to_Spher(
            std::vector<double> & u_r,
            std::vector<double> & u_lon,
            std::vector<double> & u_lat,
            const std::vector<double> & u_x,
            const std::vector<double> & u_y,
            const std::vector<double> & u_z,
            const dataset & source_data
            );


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

void filtering(const dataset & source_data,
               const std::vector<double> & scales, 
               const MPI_Comm comm = MPI_COMM_WORLD);

void filtering_helmholtz(
        const dataset & source_data,
        const std::vector<double> & scales,
        const MPI_Comm comm = MPI_COMM_WORLD
        );

void apply_filter_at_point(
        std::vector<double*> & coarse_vals,   
        std::vector<double*> & dl_coarse_vals,
        std::vector<double*> & dll_coarse_vals,
        double & dl_kernel_val,
        double & dll_kernel_val,
        const std::vector<const std::vector<double>*> & fields,
        const dataset & source_data,
        const int Itime,  const int Idepth, const int Ilat, const int Ilon,
        const int LAT_lb,
        const int LAT_ub,
        const double scale,
        const std::vector<bool> & use_mask,
        const std::vector<double> & local_kernel,
        const std::vector<double> & local_dl_kernel,
        const std::vector<double> & local_dll_kernel,
        const std::vector<double> * weight = NULL
        );

double kernel(const double distance, const double scale, const int deriv_order = 0);

double kernel_alpha(void);

void compute_vorticity_at_point(
        double & vort_r_tmp, 
        double & vort_lon_tmp, 
        double & vort_lat_tmp,
        double & div_tmp,
        double & OkuboWeiss_tmp,
        double & cyclonic_energy,
        double & anticyclonic_energy,
        double & strain_energy,
        const dataset & source_data,
        const std::vector<double> & u_r, 
        const std::vector<double> & u_lon, 
        const std::vector<double> & u_lat,
        const int Itime,  const int Idepth, const int Ilat, const int Ilon);

void compute_vorticity(
        std::vector<double> & vort_r,    
        std::vector<double> & vort_lon,    
        std::vector<double> & vort_lat,
        std::vector<double> & vel_div,
        std::vector<double> & OkuboWeiss,
        std::vector<double> & cyclonic_energy,
        std::vector<double> & anticyclonic_energy,
        std::vector<double> & strain_energy,
        const dataset & source_data,
        const std::vector<double> & u_r, 
        const std::vector<double> & u_lon, 
        const std::vector<double> & u_lat,
        const MPI_Comm comm = MPI_COMM_WORLD);

void apply_filter_at_point_for_quadratics(
        double & uxux_tmp,    double & uxuy_tmp,    double & uxuz_tmp,
        double & uyuy_tmp,    double & uyuz_tmp,    double & uzuz_tmp,
        double & vort_ux_tmp, double & vort_uy_tmp, double & vort_uz_tmp,
        const std::vector<double> & u_x, 
        const std::vector<double> & u_y, 
        const std::vector<double> & u_z,
        const std::vector<double> & vort_r,
        const dataset & source_data,
        const int Itime,  const int Idepth, const int Ilat, const int Ilon,
        const int LAT_lb, const int LAT_ub,
        const double scale,
        const std::vector<double> & local_kernel);

void compute_Pi(
        std::vector<double> & energy_transfer,
        const dataset & source_data,
        const std::vector<double> & ux,   
        const std::vector<double> & uy,   
        const std::vector<double> & uz,
        const std::vector<double> & uxux, 
        const std::vector<double> & uxuy, 
        const std::vector<double> & uxuz,
        const std::vector<double> & uyuy, 
        const std::vector<double> & uyuz, 
        const std::vector<double> & uzuz,
        const MPI_Comm comm = MPI_COMM_WORLD);

void compute_Pi_shift_deriv(
        std::vector<double> & energy_transfer,
        const dataset & source_data,
        const std::vector<double> & ux,   
        const std::vector<double> & uy,   
        const std::vector<double> & uz,
        const std::vector<double> & uxux, 
        const std::vector<double> & uxuy, 
        const std::vector<double> & uxuz,
        const std::vector<double> & uyuy, 
        const std::vector<double> & uyuz, 
        const std::vector<double> & uzuz,
        const MPI_Comm comm = MPI_COMM_WORLD);

void compute_Pi_Helmholtz(
        std::vector<double> & energy_transfer,
        const dataset & source_data,
        const std::vector<double> & ulon,
        const std::vector<double> & ulat,
        const std::vector<double> & ulon_ulon,
        const std::vector<double> & ulon_ulat,
        const std::vector<double> & ulat_ulat,
        const MPI_Comm comm = MPI_COMM_WORLD
        );

void compute_Z(
        std::vector<double> & enstrophy_transfer,
        const dataset & source_data,
        const std::vector<double> & ux,
        const std::vector<double> & uy,
        const std::vector<double> & uz,
        const std::vector<double> & coarse_vort_r,
        const std::vector<double> & vort_ux,
        const std::vector<double> & vort_uy,
        const std::vector<double> & vort_uz,
        const MPI_Comm comm = MPI_COMM_WORLD);

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

double depotential_temperature( 
        const double p, 
        const double theta);

double equation_of_state(
        const double T,
        const double S,
        const double p);

void compute_div_transport(
        std::vector<double> & div_J,
        const dataset & source_data,
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
        const MPI_Comm comm = MPI_COMM_WORLD
        );

void compute_spatial_average(
        std::vector<double> & means,
        const std::vector<double> & field,
        const std::vector<double> & areas,
        const int Ntime,
        const int Ndepth,
        const int Nlat,
        const int Nlon,
        const std::vector<bool> & mask);

void compute_KE_spectra_and_slopes( 
        std::vector<double> & u_spectrum_tot, 
        std::vector<double> & u_spectrum_tor, 
        std::vector<double> & u_spectrum_pot,
        std::vector<double> & v_spectrum_tot, 
        std::vector<double> & v_spectrum_tor, 
        std::vector<double> & v_spectrum_pot,
        std::vector<double> & spec_slope_tot, 
        std::vector<double> & spec_slope_tor, 
        std::vector<double> & spec_slope_pot,
        const std::vector<double> & u_lon_tot, 
        const std::vector<double> & u_lon_tor, 
        const std::vector<double> & u_lon_pot,
        const std::vector<double> & u_lat_tot, 
        const std::vector<double> & u_lat_tor, 
        const std::vector<double> & u_lat_pot,
        const std::vector<double> & dl_coarse_Phi, 
        const std::vector<double> & dl_coarse_Psi,
        const std::vector<double> & dll_coarse_Phi, 
        const std::vector<double> & dll_coarse_Psi,
        const dataset & source_data,
        const double filter_scale
        );

void get_lat_bounds(
        int & LAT_lb,
        int & LAT_ub,
        const std::vector<double> & latitude,
        const int Ilat,
        const double scale);

void get_lon_bounds(
        int & LON_lb,
        int & LON_ub,
        const std::vector<double> & longitude,
        const int Ilon,
        const double centre_lat,
        const double curr_lat,
        const double scale);

void print_compile_info(
        const std::vector<double> * scales = NULL);

int get_omp_chunksize(const int Nlat, const int Nlon);


void convert_coordinates(
        std::vector<double> & longitude,
        std::vector<double> & latitude
        );

void mask_out_pole(
        const std::vector<double> & latitude,
        std::vector<bool> & mask,
        const int Ntime,
        const int Ndepth,
        const int Nlat,
        const int Nlon
        );


void roll_field(
        std::vector<double> & field_to_roll,
        const std::string dimension,
        const int roll_count,
        const int Ntime,
        const int Ndepth,
        const int Nlat,
        const int Nlon
        );


// Extending to the poles

void extend_latitude_to_poles(
        const std::vector<double> & original_latitude,
        std::vector<double> & extended_latitude,
        int & orig_Ilat_start,
        const bool IS_DEGREES = false,
        const MPI_Comm comm = MPI_COMM_WORLD
        );

void extend_field_to_poles(
        std::vector<double> & array_to_extend,
        const dataset & source_data,
        const std::vector<double> & extended_latitude,
        const int Ilat_start
        );

void extend_mask_to_poles(
        std::vector<bool> & mask_to_extend,
        const dataset & source_data,
        const std::vector<double> & extended_latitude,
        const int Ilat_start,
        const bool extend_val = constants::FILTER_OVER_LAND
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
        void add_to_record( const double delta, const std::string record_name );

        //! Print the timing information in a human-readable format.
        void print() const;

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

        /*!
         *  If command-line option <option> was passed at run-time, return the value, otherwise return the default
         */
        const std::string getCmdOption(
                const std::string &option,
                const std::string &default_value = "",
                const bool help = false,
                const std::string &description = ""
                ) const;

        /*!
         *  Check if command-line option <option> was passed at run-time
         */
        bool cmdOptionExists(const std::string &option) const;

        /*!
         *  Extract, parse, and format, the string of filter scales that was provided at run-time.
         *  Assumes a string of space-delimited numbers (e.g. "1.3e3 2.4e5 8.9e9")
         */
        void getFilterScales( 
                std::vector<double> &filter_scales, 
                const std::string &argname,
                const bool help = false
                ) const;

        /*!
         *  Extract, parse, and format, a list of strings. This is typically a list of variable names.
         *  Assumes a string of space-delimited strings (e.g. "rho u v")
         */
        void getListofStrings( 
                std::vector<std::string> &list_of_strings, 
                const std::string &argname,
                const bool help = false,
                const std::string &description = ""
                ) const;

    private:
        std::vector <std::string> tokens;

};

/*!
 * \brief Convert string to boolean
 */
bool string_to_bool( std::string str );

void print_header_info( );

#endif
