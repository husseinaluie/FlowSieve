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

        // Load in region definitions
        void load_region_definitions(   const std::string filename, 
                                        const std::string dim_name, 
                                        const std::string var_name, 
                                        const MPI_Comm = MPI_COMM_WORLD );
        void compute_region_areas();

        // Check the processors divions between dimensions
        void check_processor_divisions( const int Nprocs_in_time_input, const int Nprocs_in_depth_input, const MPI_Comm = MPI_COMM_WORLD );

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
        const double scale,
        const std::vector<double> & longitude,
        const std::vector<double> & latitude,
        const int Ilat,     const int Ilon,
        const int Ntime,    const int Ndepth,
        const int Nlat,     const int Nlon,
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
            const std::vector<bool> & mask,
            const std::vector<double> & time,
            const std::vector<double> & depth,
            const std::vector<double> & latitude,
            const std::vector<double> & longitude
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
            const std::vector<bool> & mask,
            const std::vector<double> & time,
            const std::vector<double> & depth,
            const std::vector<double> & latitude,
            const std::vector<double> & longitude
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

double kernel(const double distance, const double scale);

double kernel_alpha(void);

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

void compute_spatial_average(
        std::vector<double> & means,
        const std::vector<double> & field,
        const std::vector<double> & areas,
        const int Ntime,
        const int Ndepth,
        const int Nlat,
        const int Nlon,
        const std::vector<bool> & mask);

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
