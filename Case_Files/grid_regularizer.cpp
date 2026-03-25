#include <fenv.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <algorithm>
#include <math.h>
#include <vector>
#include <mpi.h>
#include <omp.h>
#include <cassert>
#include <locale>
#include <iostream>

#include "../netcdf_io.hpp"
#include "../functions.hpp"
#include "../constants.hpp"
#include "../preprocess.hpp"

template<unsigned int power>
double my_pow(double x) {
    if      constexpr (power == 0) { return 1.0; }
    else if constexpr (power == 1) { return x; }
    else                           { return x * my_pow<power-1>(x); }
}


int main(int argc, char *argv[]) {
    
    // PERIODIC_Y implies UNIFORM_LAT_GRID
    static_assert( (constants::UNIFORM_LAT_GRID) or (not(constants::PERIODIC_Y)),
            "PERIODIC_Y requires UNIFORM_LAT_GRID.\n"
            "Please update constants.hpp accordingly.\n");

    // Currently cannot be Cartesian
    static_assert( not(constants::CARTESIAN),
            "Toroidal projection not set to handle Cartesian coordinates.\n"
            );

    // Enable all floating point exceptions but FE_INEXACT
    //feenableexcept(FE_ALL_EXCEPT & ~FE_INEXACT);

    // Specify the number of OpenMP threads
    //   and initialize the MPI world
    int thread_safety_provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &thread_safety_provided);
    //MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI::ERRORS_THROW_EXCEPTIONS);

    int wRank=-1, wSize=-1;
    MPI_Comm_rank( MPI_COMM_WORLD, &wRank );
    MPI_Comm_size( MPI_COMM_WORLD, &wSize );

    // Set the locale to the system/user specs
    //  [ this modifies how numbers are printed ]
    // on startup, the global locale is the "C" locale. we'll replace the C++ global 
    // locale and the "C" locale with the user-preferred locale for future wide character output
    //std::cout << L"User-preferred locale setting is " << std::locale("").name().c_str() << L'\n';
    std::locale::global(std::locale(""));
    std::cout.imbue(std::locale());

    //
    //// Parse command-line arguments
    //
    InputParser input(argc, argv);
    if(input.cmdOptionExists("--version")){
        if (wRank == 0) { print_compile_info(NULL); } 
        return 0;
    }
    const bool asked_help = input.cmdOptionExists("--help");
    if (asked_help) {
        fprintf( stdout, "\033[1;4mThe command-line input arguments [and default values] are:\033[0m\n" );
    }

    // first argument is the flag, second argument is default value (for when flag is not present)
    const std::string   &input_fname      = input.getCmdOption("--input_file",      
                                                               "input.nc",                 
                                                               asked_help,
                                                               "netCDF file containing vector field to Helmholtz decompose"),
                        &adjacency_fname  = input.getCmdOption("--adjacency_file",     
                                                               "adjacency.nc",  
                                                               asked_help,
                                                               "Filename for the adjacency data (from the LLC_build_adjacency routine)"),
                        &output_fname     = input.getCmdOption("--output_file",     
                                                               "regularized_grid.nc",  
                                                               asked_help,
                                                               "Filename for the output");

    const std::string   &time_dim_name      = input.getCmdOption("--time",        
                                                                 "time",       
                                                                 asked_help,
                                                                 "Name of 'time' dimension in netCDF input file."),
                        &depth_dim_name     = input.getCmdOption("--depth",       
                                                                 "depth",      
                                                                 asked_help,
                                                                 "Name of 'depth' dimension in netCDF input file."),
                        &latitude_dim_name  = input.getCmdOption("--latitude",    
                                                                 "latitude",   
                                                                 asked_help,
                                                                 "Name of 'latitude' dimension in netCDF input file."),
                        &longitude_dim_name = input.getCmdOption("--longitude",   
                                                                 "longitude",  
                                                                 asked_help,
                                                                 "Name of 'longitude' dimension in netCDF input file.");

    const std::string &latlon_in_degrees  = input.getCmdOption("--is_degrees",   
                                                               "true", 
                                                               asked_help,
                                                               "Boolean (true/false) indicating if the grid is in degrees (true) or radians (false).");

    const std::string   &Nprocs_in_time_string  = input.getCmdOption("--Nprocs_in_time",  
                                                                     "1", 
                                                                     asked_help,
                                                                     "The number of MPI divisions in time. Optimally divides Ntime evenly.\nIf Ndepth = 1, Nprocs_in_time is automatically determined."),
                        &Nprocs_in_depth_string = input.getCmdOption("--Nprocs_in_depth", 
                                                                     "1", 
                                                                     asked_help,
                                                                     "The number of MPI divisions in depth. Optimally divides Ndepth evenly.\nIf Ntime = 1, Nprocs_in_depth is automatically determined.");
    const int   Nprocs_in_time_input  = stoi(Nprocs_in_time_string),
                Nprocs_in_depth_input = stoi(Nprocs_in_depth_string);

    const std::string   &zonal_vel_name    = input.getCmdOption("--zonal_vel",   
                                                                "uo",          
                                                                asked_help,
                                                                "Name of zonal (eastward) velocity in netCDF input file"),
                        &merid_vel_name    = input.getCmdOption("--merid_vel",   
                                                                "vo",          
                                                                asked_help,
                                                                "Name of meridional (northward) velocity in netCDF input file"),
                        &dArea_field_var_name   = input.getCmdOption("--dArea_field",   "dA",    asked_help, "Name of cell areas in input file."),
                        &tor_seed_name     = input.getCmdOption("--tor_seed",    
                                                                "Psi_seed",    
                                                                asked_help,
                                                                "Name of streamfunction (Psi) seed variable in seed file (if applicable)"),
                        &pot_seed_name     = input.getCmdOption("--pot_seed",    
                                                                "Phi_seed",    
                                                                asked_help,
                                                                "Name of potential function (Phi) seed variable in seed file (if applicable)");

    const std::string &iteration_string = input.getCmdOption("--max_iterations", 
                                                             "1e5", 
                                                             asked_help,
                                                             "Maximum number of iterations for the solver before terminating. Can use exponential notation (e.g. '5e3')");
    const size_t max_iterations = (size_t) stod(iteration_string);  

    const std::string &polar_fix_string = input.getCmdOption("--south_polar_fix", 
                                                             "0", 
                                                             asked_help,
                                                             "Number of grid points to hold fixed for south pole");
    const size_t polar_fix_count = (size_t) stod(polar_fix_string);  

    // repulsion parameter
    const std::string &repulsion_beta_string = input.getCmdOption("--repulsion", 
                                                       "0", 
                                                       asked_help,
                                                       "Repulsion coefficient [0,1].");
    const double repulsion_beta = stod(repulsion_beta_string);  


    if (asked_help) { return 0; }

    // Print processor assignments
    const int max_threads = omp_get_max_threads();
    omp_set_num_threads( max_threads );

    // Print some header info, depending on debug level
    print_header_info();

    // Initialize dataset class instance
    dataset source_data;

    // Read in source data / get size information
    #if DEBUG >= 1
    if (wRank == 0) { fprintf(stdout, "Reading in source data.\n\n"); }
    #endif

    // Read in the grid coordinates
    source_data.load_time(      time_dim_name,      input_fname );
    source_data.load_depth(     depth_dim_name,     input_fname );
    read_LLC_latlon_from_file( source_data.latitude,  latitude_dim_name,  input_fname );
    read_LLC_latlon_from_file( source_data.longitude, longitude_dim_name, input_fname );
     
    // Convert to radians, if appropriate
    if ( (latlon_in_degrees == "true") and (not(constants::CARTESIAN)) ) {
        convert_coordinates( source_data.longitude, source_data.latitude );
    }

    // Apply some cleaning to the processor allotments if necessary. 
    source_data.Nlat = 1;
    source_data.Nlon = 1;
    source_data.check_processor_divisions( Nprocs_in_time_input, Nprocs_in_depth_input );

    // Read in the velocity fields
    source_data.load_variable( "u_lon", zonal_vel_name, input_fname, true, true );
    source_data.load_variable( "u_lat", merid_vel_name, input_fname, true, true );

    // If we're using FILTER_OVER_LAND, then the mask has been wiped out. Load in a mask that still includes land references
    //      so that we have both. Will be used to get 'water-only' region areas.
    if (constants::FILTER_OVER_LAND) { 
        read_mask_from_file( source_data.reference_mask, zonal_vel_name, input_fname,
               source_data.Nprocs_in_time, source_data.Nprocs_in_depth );
    }

    // Get the MPI-local dimension sizes
    source_data.Ntime  = source_data.myCounts[0];
    source_data.Ndepth = source_data.myCounts[1];

    // Compute the area of each 'cell' which will be necessary for integration
    //source_data.compute_cell_areas();
    // Cell areas are trickier, so they will be passed in as an input.
    read_LLC_latlon_from_file( source_data.areas, dArea_field_var_name, adjacency_fname );

    // Build the adjacency matrix and other adjacency-adjacent arrays
    //  down the road, just load in a pre-built one, but for right now
    //  this is easier.
    // Need to have area loaded already!
    source_data.load_adjacency( adjacency_fname );

    size_t Ipt;
    const size_t Npts = source_data.longitude.size();

    // Flag points that are sufficiently far from coast
    std::vector<short int> ill_conditioned( Npts, false );
    size_t num_ill_conditioned = 0;
    #pragma omp parallel default(none) private( Ipt ) \
    shared( source_data, ill_conditioned ) \
    firstprivate(Npts) reduction( +:num_ill_conditioned )
    {
        #pragma omp for collapse(1) schedule(guided) 
        for ( Ipt = 0; Ipt < Npts; Ipt++ ) {

            if ( source_data.mask[Ipt] ) { ill_conditioned[Ipt] = false; continue; }

            bool all_land_neighbours = true;
            for ( size_t In = 0; In < source_data.num_neighbours; In++ ) {
                size_t Ineigh = source_data.adjacency_indices[Ipt][In];

                if ( source_data.mask[Ineigh] ) { all_land_neighbours = false; break; }

                for ( size_t Jn = 0; Jn < source_data.num_neighbours; Jn++ ) {
                    size_t Jneigh = source_data.adjacency_indices[Ineigh][Jn];

                    if ( source_data.mask[Jneigh] ) { all_land_neighbours = false; break; }
                }
            }
            ill_conditioned[Ipt] = all_land_neighbours;
            if ( ill_conditioned[Ipt] ) { num_ill_conditioned++; }
        }
    }
    /*
    for ( Ipt = 0; Ipt < polar_fix_count; Ipt++ ) {
        ill_conditioned[Ipt] = false; // fix the south polar ring if specified
        num_ill_conditioned--;
    }
    */
    fprintf( stdout, "%'zu points flagged as ill conditioned\n", num_ill_conditioned );
    if ( num_ill_conditioned == 0 ) {
        throw std::logic_error("No ill conditioned points detected! Halting now.\n\n");
    }

    // Get Cartesian coordiantes
    std::vector<double> X(Npts,0.), Y(Npts,0.), Z(Npts,0.), X_tmp(Npts,0.), Y_tmp(Npts,0.), Z_tmp(Npts,0.);
    #pragma omp parallel default(none) private( Ipt ) \
    shared( source_data, ill_conditioned, X, Y, Z, X_tmp, Y_tmp, Z_tmp ) \
    firstprivate(Npts)
    {
        #pragma omp for collapse(1) schedule(guided) 
        for ( Ipt = 0; Ipt < Npts; Ipt++ ) {
            if ( source_data.mask[Ipt] ) { continue; }
            // Do for all land points

            double loc_lat = source_data.latitude[Ipt];
            double loc_lon = source_data.longitude[Ipt];
            X[Ipt] = cos( loc_lat ) * cos( loc_lon );
            Y[Ipt] = cos( loc_lat ) * sin( loc_lon );
            Z[Ipt] = sin( loc_lat );

            X_tmp[Ipt] = X[Ipt];
            Y_tmp[Ipt] = Y[Ipt];
            Z_tmp[Ipt] = Z[Ipt];
        }
    }

    // For desired number of iterations, run regularizer
    int counter = 1000 - 1;
    const double d_epsilon = 1. / 6371e3;
    for ( size_t iter = 0; iter < max_iterations; iter++ ) {
        counter++;
        long double mean_displacement = 0.;

        #pragma omp parallel default(none) private( Ipt ) \
        shared( source_data, ill_conditioned, X, Y, Z, X_tmp, Y_tmp, Z_tmp ) \
        firstprivate(Npts, counter, num_ill_conditioned, repulsion_beta, max_iterations ) \
        reduction( +:mean_displacement )
        {

            #pragma omp for collapse(1) schedule(static,1)
            for ( Ipt = 0; Ipt < Npts; Ipt++ ) {

                bool do_repulse = true;
                //bool do_repulse = (source_data.latitude[Ipt]) > (-60. * M_PI / 180.);

                // Do for all land points with land-only neighbours
                if ( ill_conditioned[Ipt] ) {
                    double x_bar = 0., y_bar = 0., z_bar = 0.;
                    double x_disp = 0., y_disp = 0., z_disp = 0.;
                    double mean_delta = 0.;

                    const double X0 = X_tmp[Ipt];
                    const double Y0 = Y_tmp[Ipt];
                    const double Z0 = Z_tmp[Ipt];

                    const size_t UB = source_data.num_neighbours + 1;
                    const double UB_inv = 1. / (double) UB;
                    //const size_t UB = source_data.num_neighbours + ( do_repulse ? 0 : 1 );
                        //  when repulsion is on, just looks at neghbouring points
                        //  when repulsion is off, includes self, for a stabilized laplace
                    for ( size_t In = 0; In < UB; In++ ) {

                        size_t Ineigh = source_data.adjacency_indices[Ipt][In];

                        // Build the centroid of the neihbours
                        x_bar += X[Ineigh] * UB_inv; //source_data.num_neighbours;
                        y_bar += Y[Ineigh] * UB_inv; //source_data.num_neighbours;
                        z_bar += Z[Ineigh] * UB_inv; //source_data.num_neighbours;

                        if ( do_repulse ) {
                            // Also build a 'repulsive' force
                            double x_delta = X0 - X[Ineigh];
                            double y_delta = Y0 - Y[Ineigh];
                            double z_delta = Z0 - Z[Ineigh];

                            double d_delta = sqrt( my_pow<2>(x_delta) + my_pow<2>(y_delta) + my_pow<2>(z_delta) );
                            mean_delta += d_delta / source_data.num_neighbours;

                            //const double d_cube = pow( d_delta + d_epsilon, 3. );
                            /*
                            const double d_cube =   ( d_delta + d_epsilon ) 
                                                  * ( d_delta + d_epsilon ) 
                                                  * ( d_delta + d_epsilon );
                            */
                            const double d_cube_inv = 1. / my_pow<3>( d_delta + d_epsilon );
                            x_disp += ( x_delta * d_cube_inv ) / source_data.num_neighbours;
                            y_disp += ( y_delta * d_cube_inv ) / source_data.num_neighbours;
                            z_disp += ( z_delta * d_cube_inv ) / source_data.num_neighbours;
                        }
                    }

                    // Remove the radial laplace term
                    const double radial_laplace = x_bar * X0 + y_bar * Y0 + z_bar * Z0;
                    double x_Lap = x_bar - radial_laplace * X0;
                    double y_Lap = y_bar - radial_laplace * Y0;
                    double z_Lap = z_bar - radial_laplace * Z0;

                    /*
                    // Normalize the mean onto the sphere
                    const double r = sqrt( x_bar*x_bar + y_bar*y_bar + z_bar*z_bar );

                    x_bar /= r;
                    y_bar /= r;
                    z_bar /= r;
                    */

                    double x_repulse = 0, y_repulse = 0, z_repulse = 0;
                    if ( do_repulse ) {
                        // Remove the radial repulsion piece
                        double radial_repulsion = x_disp * X0 + y_disp * Y0 + z_disp * Z0;
                        x_repulse = x_disp - radial_repulsion * X0;
                        y_repulse = y_disp - radial_repulsion * Y0;
                        z_repulse = z_disp - radial_repulsion * Z0;

                        // Limit the repulsion to based on mean neighbour displacement
                        double repulse_mag = sqrt( x_repulse*x_repulse + y_repulse*y_repulse + z_repulse*z_repulse );
                        double reuplse_coeff = repulsion_beta * mean_delta / repulse_mag;
                        if ( repulse_coeff < 1 ) {
                            x_repulse *= repulse_coeff;
                            y_repulse *= repulse_coeff;
                            z_repulse *= repulse_coeff;
                        }
                    }


                    // Iterate to new point
                    const double delta_t = 0.1;
                    double x_new = X0 + delta_t * ( x_Lap + x_repulse );
                    double y_new = Y0 + delta_t * ( y_Lap + y_repulse );
                    double z_new = Z0 + delta_t * ( z_Lap + z_repulse );

                    // Normalize the mean onto the sphere
                    const double r_new = sqrt( x_new*x_new + y_new*y_new + z_new*z_new );

                    x_new /= r_new;
                    y_new /= r_new;
                    z_new /= r_new;

                    // Track displacements
                    if (counter == 1000) {
                        x_disp = X0 - x_new,
                        y_disp = Y0 - y_new,
                        z_disp = Z0 - z_new;
                        double d_disp = sqrt( x_disp*x_disp + y_disp*y_disp + z_disp*z_disp );
                        mean_displacement += d_disp / num_ill_conditioned;
                    }

                    // Record new points
                    X_tmp[Ipt] = x_new;
                    Y_tmp[Ipt] = y_new;
                    Z_tmp[Ipt] = z_new;

                }
            }
        }

        if (counter == 1000) {
            fprintf( stdout, "Iter %'zu: mean displacement -> %.6Le\n", iter, mean_displacement );
            counter = 0;
        }

        X.swap(X_tmp);
        Y.swap(Y_tmp);
        Z.swap(Z_tmp);
    }

    // Now reconstruct the lat/lon coordinates
    #pragma omp parallel default(none) private( Ipt ) \
    shared( source_data, ill_conditioned, X, Y, Z ) \
    firstprivate(Npts)
    {
        #pragma omp for collapse(1) schedule(guided) 
        for ( Ipt = 0; Ipt < Npts; Ipt++ ) {
            // Do for all land points with land-only neighbours
            if ( ill_conditioned[Ipt] ) {

                source_data.longitude[Ipt] = atan2( Y[Ipt], X[Ipt] );
                source_data.latitude[ Ipt] = asin( Z[Ipt] );
            }
        }
    }

    // Compute cell areas
    source_data.areas.resize(Npts);
    #pragma omp parallel default(none) private( Ipt ) \
    shared( source_data, ill_conditioned ) \
    firstprivate(Npts)
    {
        #pragma omp for collapse(1) schedule(static) 
        for ( Ipt = 0; Ipt < Npts; Ipt++ ) {

            double min_lon = 1e3, max_lon = -1e3, min_lat = 1e3, max_lat = -1e3;
            for ( size_t In = 0; In < source_data.num_neighbours; In++ ) { 
                size_t Ineigh = source_data.adjacency_indices[Ipt][In];
                double lat = source_data.latitude[Ineigh];

                min_lat = std::fmin( min_lat, lat );
                max_lat = std::fmax( max_lat, lat );

                double lon = source_data.longitude[Ineigh] - source_data.longitude[Ipt];
                if ( lon >  M_PI ) { lon -= 2*M_PI; }
                if ( lon < -M_PI ) { lon += 2*M_PI; }

                min_lon = std::fmin( min_lon, lon );
                max_lon = std::fmax( max_lon, lon );
            }

            double dlon = std::fmax( 1./constants::R_earth, max_lon - min_lon );
            double dlat = std::fmax( 1./constants::R_earth, max_lat - min_lat );

            double lat = source_data.latitude[Ipt];
            double area = 0.25 * pow(constants::R_earth, 2.) * dlon * dlat * cos( lat );

            source_data.areas[Ipt] = area;

        }
    }

    // And write to a file
    source_data.write_adjacency( output_fname.c_str() );


    // Done!
    #if DEBUG >= 0
    if (wRank == 0) {
        fprintf(stdout, "\n\n");
        fprintf(stdout, "Process completed.\n");
        fprintf(stdout, "\n");
    }
    #endif

    #if DEBUG >= 1
    fprintf(stdout, "Processor %d / %d waiting to finalize.\n", wRank + 1, wSize);
    #endif 
    MPI_Finalize();
    return 0;
}
