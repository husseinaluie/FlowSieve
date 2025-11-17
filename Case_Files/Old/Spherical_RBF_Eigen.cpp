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

#include "../netcdf_io.hpp"
#include "../functions.hpp"
#include "../constants.hpp"

#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/IterativeSolvers>

int main(int argc, char *argv[]) {
    
    // PERIODIC_Y implies UNIFORM_LAT_GRID
    static_assert( (constants::UNIFORM_LAT_GRID) or (not(constants::PERIODIC_Y)),
            "PERIODIC_Y requires UNIFORM_LAT_GRID.\n"
            "Please update constants.hpp accordingly.\n");

    // Cannot extend to poles AND be Cartesian
    static_assert( not( (constants::EXTEND_DOMAIN_TO_POLES) and (constants::CARTESIAN) ),
            "Cartesian implies that there are no poles, so cannot extend to poles."
            "Please update constants.hpp accordingly.");

    // We need to keep the mask, so cannot have constants::FILTER_OVER_LAND turned on
    static_assert( not( constants::FILTER_OVER_LAND ),
            "Spherical RBF requires the land mask, so FILTER_OVER_LAND must be turned OFF" );

    // Enable all floating point exceptions but FE_INEXACT
    //feenableexcept( FE_ALL_EXCEPT & ~FE_INEXACT & ~FE_UNDERFLOW );
    //fprintf( stdout, " %d : %d \n", FE_ALL_EXCEPT, FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW | FE_INEXACT | FE_UNDERFLOW );
    feenableexcept( FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW );

    // Specify the number of OpenMP threads
    //   and initialize the MPI world
    int thread_safety_provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &thread_safety_provided);
    //MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI::ERRORS_THROW_EXCEPTIONS);

    int wRank=-1, wSize=-1;
    MPI_Comm_rank( MPI_COMM_WORLD, &wRank );
    MPI_Comm_size( MPI_COMM_WORLD, &wSize );

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
        fprintf( stdout, "The command-line input arguments [and default values] are:\n" );
    }

    // first argument is the flag, second argument is default value (for when flag is not present)
    const std::string   &input_fname   = input.getCmdOption("--input_file",     "input.nc", asked_help),
                        &output_fname  = input.getCmdOption("--output_file",    "output.nc", asked_help);

    const std::string   &time_dim_name      = input.getCmdOption("--time",        "time",      asked_help),
                        &depth_dim_name     = input.getCmdOption("--depth",       "depth",     asked_help),
                        &latitude_dim_name  = input.getCmdOption("--latitude",    "latitude",  asked_help),
                        &longitude_dim_name = input.getCmdOption("--longitude",   "longitude", asked_help);

    const std::string &latlon_in_degrees  = input.getCmdOption("--is_degrees",   "true", asked_help);

    const std::string   &Nprocs_in_time_string      = input.getCmdOption("--Nprocs_in_time",  "1", asked_help),
                        &Nprocs_in_depth_string     = input.getCmdOption("--Nprocs_in_depth", "1", asked_help),
                        &num_interp_layers_string   = input.getCmdOption("--Num_interp_layers", "1", asked_help),
                        &max_interp_scale_string    = input.getCmdOption("--max_interp_scale", "2645e3", asked_help),
                        &skip_factor_string         = input.getCmdOption("--skip_factor", "10", asked_help),
                        &max_iterations_string      = input.getCmdOption("--max_iterations", "1e3", asked_help),
                        &tolerance_string           = input.getCmdOption("--tolerance", "1e-4", asked_help);
    const int   Nprocs_in_time_input  = stoi(Nprocs_in_time_string),
                Nprocs_in_depth_input = stoi(Nprocs_in_depth_string),
                num_interp_layers     = stoi(num_interp_layers_string),
                base_skip_factor      = stoi(skip_factor_string),
                max_iterations        = stod(max_iterations_string);

    const double max_interp_scale = stod(max_interp_scale_string),
                 tolerance        = stod(tolerance_string);

    // Also read in the fields to be filtered from commandline
    //   e.g. --variables "rho salinity p" (names must match with input netcdf file)
    std::vector< std::string > vars_to_interpolate;
    input.getListofStrings( vars_to_interpolate, "--variables", asked_help );
    const size_t Nvars = vars_to_interpolate.size();

    if (asked_help) { return 0; }

    // Print processor assignments
    const int max_threads = omp_get_max_threads();
    omp_set_num_threads( max_threads );

    // Print some header info, depending on debug level
    print_header_info();

    Timing_Records timing_records;

    // Initialize dataset class instance
    dataset source_data;

    int nEigen = Eigen::nbThreads( );
    fprintf( stdout, "Eigen will use %d threads.\n", nEigen );

    // Read in source data / get size information
    #if DEBUG >= 1
    if (wRank == 0) { fprintf(stdout, "Reading in source data.\n\n"); }
    #endif

    // Read in the grid coordinates
    //   implicitely assume coordinates are the same between input files
    source_data.load_time(      time_dim_name,      input_fname );
    source_data.load_depth(     depth_dim_name,     input_fname );
    source_data.load_latitude(  latitude_dim_name,  input_fname );
    source_data.load_longitude( longitude_dim_name, input_fname );

    // Apply some cleaning to the processor allotments if necessary. 
    source_data.check_processor_divisions( Nprocs_in_time_input, Nprocs_in_depth_input );

    // Convert to radians, if appropriate
    if ( latlon_in_degrees == "true" ) { convert_coordinates( source_data.longitude, source_data.latitude ); }

    // Compute the area of each 'cell' which will be necessary for integration
    #if DEBUG >= 2
    if (wRank == 0) { fprintf( stdout, "Computing cell areas.\n" ); }
    #endif
    source_data.compute_cell_areas();

    // Read in the scalar fields to filter
    #if DEBUG >= 1
    if (wRank == 0) { fprintf( stdout, "Reading in original fields.\n" ); }
    #endif
    for (size_t field_ind = 0; field_ind < vars_to_interpolate.size(); field_ind++) {
        source_data.load_variable( vars_to_interpolate.at(field_ind), vars_to_interpolate.at(field_ind), input_fname, true, true );
    }

    // Get the MPI-local dimension sizes
    source_data.Ntime  = source_data.myCounts[0];
    source_data.Ndepth = source_data.myCounts[1];
    const size_t Npts = source_data.Ntime * source_data.Ndepth * source_data.Nlat * source_data.Nlon;

    const std::vector<int>  &myStarts = source_data.myStarts;

    //
    const int   Ntime   = source_data.Ntime,
                Ndepth  = source_data.Ndepth,
                Nlat    = source_data.Nlat,
                Nlon    = source_data.Nlon;

    const std::vector<double>   &latitude   = source_data.latitude,
                                &longitude  = source_data.longitude;

    const std::vector<bool> &mask = source_data.mask;

    // Initialize interpolation fields
    #if DEBUG >= 1
    if (wRank == 0) { fprintf( stdout, "Setting up %'zu interp fields.\n", Nvars ); fflush(stdout); }
    #endif
    std::vector< std::vector<double> > interpolated_fields(Nvars);
    for (size_t field_ind = 0; field_ind < Nvars; field_ind++) {
        interpolated_fields.at(field_ind).resize( Npts, 0. );
    }


    // Coastline bool arrays [actually a bit-array under the hood]
    std::vector<bool> is_coastal_land(Npts), is_coastal_water(Npts), is_coastal_water_ALL(Npts);

    // various indices etc
    int Itime, Idepth, Ilat, Ilon, LAT, LON, curr_lat, curr_lon;
    double lat_at_curr, lat_at_ilat;
    size_t index, inner_index, i, j, i_land, i_water, Ivar;
    bool is_coast;

    for ( int interp_layer = 0; interp_layer < num_interp_layers; interp_layer++ ) {

        // Get the length-scale for this layer of interpolation
        double interp_scale = max_interp_scale / pow( 2, interp_layer );
        //double interp_scale = max_interp_scale / pow( 2, num_interp_layers - interp_layer - 1 );
        const double grid_res = constants::R_earth * (latitude[1] - latitude[0]); // right now assume uniform and equal spacing
        const int skip_factor = std::max( 1, (int) floor(interp_scale / ( base_skip_factor * grid_res )) );

        fprintf( stdout, "\n\nApplying RBF layer %d - %gkm, with skip-factor %d\n", interp_layer+1, interp_scale/1e3, skip_factor );

        //
        //// Get the 'coastline' mask for this layer
        //
        fprintf(stdout, "Mapping the coastline.\n");
        std::fill( is_coastal_land.begin(), is_coastal_land.end(), 0);
        std::fill( is_coastal_water.begin(), is_coastal_water.end(), 0);
        std::fill( is_coastal_water_ALL.begin(), is_coastal_water_ALL.end(), 0); // this ignores skip_factor

        size_t num_coastal_land = 0,
               num_coastal_water = 0,
               num_coastal_water_ALL = 0;
        int LAT_lb, LAT_ub, LON_ub, LON_lb;
        #pragma omp parallel \
        default(none) \
        shared( source_data, mask, latitude, longitude, is_coastal_land, is_coastal_water, is_coastal_water_ALL ) \
        private( Itime, Idepth, Ilat, Ilon, is_coast, index, inner_index, \
                 LAT_lb, LAT_ub, LON_lb, LON_ub, \
                 LAT, LON, lat_at_ilat, lat_at_curr, curr_lat, curr_lon) \
        firstprivate( Ntime, Ndepth, Nlat, Nlon, interp_scale, skip_factor ) \
        reduction( +:num_coastal_land,num_coastal_water,num_coastal_water_ALL )
        {
            #pragma omp for collapse(2) schedule(guided)
            for ( Ilat = 0; Ilat < Nlat; Ilat++ ) {
                for ( Ilon = 0; Ilon < Nlon; Ilon++ ) {
                    Itime = 0;  // just hard-code for now
                    Idepth = 0;

                    is_coast = false;

                    index = Index(Itime, Idepth, Ilat, Ilon,
                                  Ntime, Ndepth, Nlat, Nlon);
                    get_lat_bounds(LAT_lb, LAT_ub, latitude,  Ilat, 
                                    interp_scale );

                    lat_at_ilat = latitude[Ilat];

                    for ( LAT = LAT_lb; LAT < LAT_ub; LAT++) {

                        // Handle periodicity if necessary
                        if (constants::PERIODIC_Y) { curr_lat = ( LAT % Nlat + Nlat ) % Nlat; }
                        else                       { curr_lat = LAT; }
                        lat_at_curr = latitude.at(curr_lat);

                        get_lon_bounds(LON_lb, LON_ub, longitude, Ilon, lat_at_ilat, lat_at_curr, 
                                        interp_scale );
                        for ( LON = LON_lb; LON < LON_ub; LON++ ) {

                            // Handle periodicity if necessary
                            if (constants::PERIODIC_X) { curr_lon = ( LON % Nlon + Nlon ) % Nlon; }
                            else                       { curr_lon = LON; }

                            inner_index = Index(Itime, Idepth, curr_lat, curr_lon, Ntime, Ndepth, Nlat, Nlon);

                            // If mask is different between the two points, then it's coastal!
                            if ( mask.at(index) xor mask.at(inner_index) ) {
                                is_coast = true;
                                break; // break out of LON loop
                            }
                        }
                        if ( is_coast ) { break; } // break out of LAT loop
                    }
                    if (is_coast) { 
                        if ( mask.at(index) ) {
                            if ( (Ilat % skip_factor == 0) and (Ilon % skip_factor == 0) ) {
                                // Only keep a subset for the solver
                                is_coastal_water.at(index) = true;
                                num_coastal_water++;
                            }
                            is_coastal_water_ALL.at(index) = true;
                            num_coastal_water_ALL++;
                        } else {
                            is_coastal_land.at(index) = true;
                            num_coastal_land++;
                        }
                    }
                }
            }
        }
        // end of computing coastal mask
        fprintf(stdout, "Mapped the coastline. %'zu land points, %'zu water points.\n",
                num_coastal_land, num_coastal_water);

        fprintf(stdout, "Building utility index arrays.\n");
        // for utility, build an array of coastal water and land indices
        std::vector< size_t > coastal_water_indices( num_coastal_water, 0 );
        index = 0;
        for ( size_t II = 0; II < Npts; II++ ) {
            if ( is_coastal_water[II] ) { 
                coastal_water_indices[index] = II; 
                index++;
            }
        }

        std::vector< size_t > coastal_water_indices_ALL( num_coastal_water_ALL, 0 );
        index = 0;
        for ( size_t II = 0; II < Npts; II++ ) {
            if ( is_coastal_water_ALL[II] ) { 
                coastal_water_indices_ALL[index] = II; 
                index++;
            }
        }

        std::vector< size_t > coastal_land_indices( num_coastal_land, 0 );
        index = 0;
        for ( size_t II = 0; II < Npts; II++ ) {
            if ( is_coastal_land[II] ) { 
                coastal_land_indices[index] = II; 
                index++;
            }
        }


        // Count number of non-zero elements going into the Aij matrix
        double lat_i, lat_j, lon_i, lon_j, kern_val;
        size_t num_nonzero = 0;
        std::vector<unsigned int> nnz( num_coastal_water, 0 );
        #pragma omp parallel default(none) \
        shared( is_coastal_water, coastal_water_indices, latitude, longitude, nnz ) \
        private( i, j, index, inner_index, Itime, Idepth, Ilat, Ilon, curr_lat, curr_lon, \
                 LAT, LON, LAT_lb, LAT_ub, LON_ub, LON_lb, lon_i, lat_i, lon_j, lat_j, \
                 kern_val ) \
        firstprivate( Ntime, Ndepth, Nlat, Nlon, num_coastal_water, \
                      interp_scale ) \
        reduction( +:num_nonzero )
        {
            #pragma omp for collapse(1) schedule(static)
            for ( i = 0; i < num_coastal_water; i++ ) {
                Itime = 0;
                Idepth = 0;
                index = coastal_water_indices[i];
                Index1to4( index, Itime, Idepth, Ilat, Ilon,
                        Ntime, Ndepth, Nlat, Nlon );
                get_lat_bounds(LAT_lb, LAT_ub, latitude, Ilat, interp_scale );

                lon_i = longitude[Ilon];
                lat_i = latitude[Ilat];

                for ( LAT = LAT_lb; LAT < LAT_ub; LAT++ ) {

                    // Handle periodicity if necessary
                    if (constants::PERIODIC_Y) { curr_lat = ( LAT % Nlat + Nlat ) % Nlat; }
                    else                       { curr_lat = LAT; }
                    lat_j = latitude.at(curr_lat);

                    get_lon_bounds(LON_lb, LON_ub, longitude, Ilon, lat_i, lat_j, interp_scale );
                    for ( LON = LON_lb; LON < LON_ub; LON++ ) {

                        // Handle periodicity if necessary
                        if (constants::PERIODIC_X) { curr_lon = ( LON % Nlon + Nlon ) % Nlon; }
                        else                       { curr_lon = LON; }
                        lon_j = longitude[curr_lon];

                        inner_index = Index(Itime, Idepth, curr_lat, curr_lon, Ntime, Ndepth, Nlat, Nlon);

                        if ( not( is_coastal_water[inner_index] ) ) { continue; }

                        //if (inner_index < index) { continue; } // only store upper triangular part

                        //kern_val = (distance( lon_i, lat_i, lon_j, lat_j ) < interp_scale / 2 ) ? 1. : 0.;
                        kern_val = kernel( distance( lon_i, lat_i, lon_j, lat_j ), interp_scale );
                        if ( kern_val > 1e-3 ) {
                            num_nonzero++;
                            nnz[i]++;
                        }
                    }
                }
            }
        }
        fprintf( stdout, "Total non-zeros: %'zu\n", num_nonzero );

        std::vector<size_t> prev_nnz(num_coastal_water, 0);
        for ( j = 1; j < num_coastal_water; j++ ) {
            prev_nnz[j] = prev_nnz[j-1] + nnz[j-1];
        }
        fprintf( stdout, "Built prev_nnz\n" );

        // We can't do parallel writing to Aij directly,
        //  so build list of the entries in parallel, and then
        //  just copy them over in serial
        size_t Innz = 0;
        typedef Eigen::Triplet<double> T;
        std::vector<T> Aij_triplets(num_nonzero);
        #pragma omp parallel default(none) \
        shared( is_coastal_water, coastal_water_indices, latitude, longitude, \
                nnz, prev_nnz, Aij_triplets, stderr ) \
        private( i, j, i_water, index, inner_index, Innz, Itime, Idepth, Ilat, Ilon, curr_lat, curr_lon, \
                 LAT, LON, LAT_lb, LAT_ub, LON_ub, LON_lb, lon_i, lat_i, lon_j, lat_j, kern_val ) \
        firstprivate( Ntime, Ndepth, Nlat, Nlon, num_coastal_water, \
                      interp_scale )
        {
            #pragma omp for collapse(1) schedule(static)
            for ( i = 0; i < num_coastal_water; i++ ) {
                Itime = 0;
                Idepth = 0;
                index = coastal_water_indices[i];
                Index1to4( index, Itime, Idepth, Ilat, Ilon,
                        Ntime, Ndepth, Nlat, Nlon );
                get_lat_bounds(LAT_lb, LAT_ub, latitude, Ilat, interp_scale );

                lon_i = longitude[Ilon];
                lat_i = latitude[Ilat];

                Innz = 0;
                for ( LAT = LAT_lb; LAT < LAT_ub; LAT++ ) {

                    // Handle periodicity if necessary
                    if (constants::PERIODIC_Y) { curr_lat = ( LAT % Nlat + Nlat ) % Nlat; }
                    else                       { curr_lat = LAT; }
                    lat_j = latitude.at(curr_lat);

                    get_lon_bounds(LON_lb, LON_ub, longitude, Ilon, lat_i, lat_j, interp_scale );
                    for ( LON = LON_lb; LON < LON_ub; LON++ ) {

                        // Handle periodicity if necessary
                        if (constants::PERIODIC_X) { curr_lon = ( LON % Nlon + Nlon ) % Nlon; }
                        else                       { curr_lon = LON; }
                        lon_j = longitude[curr_lon];

                        inner_index = Index(Itime, Idepth, curr_lat, curr_lon, Ntime, Ndepth, Nlat, Nlon);

                        if ( not( is_coastal_water[inner_index] ) ) { continue; }

                        //if (inner_index < index) { continue; } // only store upper triangular part

                        //kern_val = (distance( lon_i, lat_i, lon_j, lat_j ) < interp_scale / 2 ) ? 1. : 0.;
                        kern_val = kernel( distance( lon_i, lat_i, lon_j, lat_j ), interp_scale );
                        if ( kern_val > 1e-3 ) {

                            i_water = std::lower_bound( coastal_water_indices.begin(), 
                                                        coastal_water_indices.end(), 
                                                        inner_index ) 
                                - coastal_water_indices.begin();

                            Aij_triplets[prev_nnz[i] + Innz] = T( i, i_water, kern_val );
                            Innz++;
                        }
                    }
                }
                if (Innz != nnz[i] ) {
                    fprintf( stderr, "(i,nnz[i],Innz) = (%'zu, %'u, %'zu)\n", i, nnz[i], Innz );
                }
                assert(Innz == nnz[i]);
            }
        }

        //
        //// Build the linear system
        //
        std::vector<double> RHS_vector( num_coastal_water, 0. );
        Eigen::SparseMatrix<double> Aij( num_coastal_water, num_coastal_water );

        fprintf(stdout, "Building the system.\n");
        Aij.setFromTriplets( Aij_triplets.begin(), Aij_triplets.end() );
        //Aij_triplets.clear(); // free up the memory

        fprintf(stdout, "Converting the system to crs.\n");
        Aij.makeCompressed();

        fprintf(stdout, "  Solving the system!\n");
        //Eigen::ConjugateGradient< Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper > solver;
        Eigen::LeastSquaresConjugateGradient< Eigen::SparseMatrix<double> > solver;
        //Eigen::BiCGSTAB< Eigen::SparseMatrix<double> > solver;
        //Eigen::GMRES< Eigen::SparseMatrix<double> > solver;
        fprintf( stdout, "  Setting solver tolerance and max iterations to %g and %d.\n", tolerance, max_iterations );
        solver.setMaxIterations(max_iterations);
        solver.setTolerance(tolerance);
        solver.compute( Aij );
        if ( solver.info() != Eigen::Success ) {
            // decomposition failed
            fprintf( stderr, "Eigen decomposition failed.\n" );
            return -1;
        }

        fprintf( stdout, "len(coastal_land_indices) = %zu\n", coastal_land_indices.size() );
        fprintf( stdout, "len(coastal_water_indices) = %zu\n", coastal_water_indices.size() );

        for ( Ivar = 0; Ivar < Nvars; Ivar++ ) {
            fprintf(stdout, "Applying to variable %zu of %zu.\n", Ivar+1, Nvars);

            // Do a bit of extra processing if this is the first round
            if ( interp_layer == 0 ) {
                // Get the global mean value
                double val_sum = 0, area_sum = 0, dA;
                #pragma omp parallel default(none) \
                shared( mask, source_data, vars_to_interpolate, dA ) \
                private( i ) \
                firstprivate( Ivar, Npts ) \
                reduction( +:val_sum,area_sum )
                {
                    #pragma omp for collapse(1) schedule(static)
                    for ( i = 0; i < Npts; i++ ) {
                        if ( mask[i] ) {
                            dA = source_data.areas.at(i);
                            val_sum += dA * source_data.variables.at(vars_to_interpolate[Ivar])[i];
                            area_sum += dA;
                        }
                    }
                }
                const double global_mean_val = val_sum / area_sum;

                // For anything outside of the largest coastline, just copy over the original
                // Otherwise, initialize with global mean
                #pragma omp parallel default(none) \
                shared( mask, source_data, vars_to_interpolate, \
                        interpolated_fields, is_coastal_water_ALL ) \
                private( i ) \
                firstprivate( Ivar, Npts, global_mean_val )
                {
                    #pragma omp for collapse(1) schedule(static)
                    for ( i = 0; i < Npts; i++ ) {
                        interpolated_fields[Ivar][i] = global_mean_val;
                        /*
                        if ( mask[i] and not( is_coastal_water_ALL[i] ) ) {
                            // if water, but not coastal water
                            interpolated_fields[Ivar][i] = source_data.variables.at(vars_to_interpolate[Ivar])[i];
                        } else {
                            interpolated_fields[Ivar][i] = global_mean_val;
                        }
                        */
                    }
                }

                // And finally, subtract the global mean off from the original field
                #pragma omp parallel default(none) \
                shared( mask, source_data, vars_to_interpolate ) \
                private( i ) \
                firstprivate( Ivar, Npts, global_mean_val )
                {
                    #pragma omp for collapse(1) schedule(static)
                    for ( i = 0; i < Npts; i++ ) {
                        if ( mask[i] ) {
                            source_data.variables.at(vars_to_interpolate[Ivar])[i] -= global_mean_val;
                        }
                    }
                }
            }

            // Set up RHS vector
            #pragma omp parallel default(none) \
            shared(RHS_vector, source_data, vars_to_interpolate, coastal_water_indices) \
            private( i ) \
            firstprivate( Ivar, num_coastal_water )
            {
                #pragma omp for collapse(1) schedule(static)
                for ( i = 0; i < num_coastal_water; i++ ) {
                    RHS_vector[i] = source_data.variables.at(vars_to_interpolate[Ivar]
                                                            )[coastal_water_indices[i]];
                }
            }

            Eigen::Map<Eigen::VectorXd> RHS( &RHS_vector[0], num_coastal_water );

            // Apply the solved linear system
            fprintf(stdout, "  Applying the solved system.\n");
            Eigen::VectorXd F_Eigen = solver.solve( RHS );
            /*
            if ( solver.info() != Eigen::Success ) {
                // solving failed
                fprintf( stderr, "Eigen solving failed.\n" );
                return -1;
            } else {
                fprintf( stdout, "Solver converged after %ld iterations to error %g.\n", 
                         solver.iterations(), solver.error() );
            }
            */
            fprintf( stdout, "Solver converged after %ld iterations to error %g.\n", 
                    solver.iterations(), solver.error() );
            std::vector<double> wj( F_Eigen.data(), F_Eigen.data() + num_coastal_water);
            fprintf( stdout, " size of wj is %zu\n", wj.size() );
    
            // And now apply the interpolation weights to fill in the water parts
            fprintf(stdout, "  Applying interpolation weights.\n");
            double lat_land, lon_land, lat_water, lon_water;
            #pragma omp parallel \
            default(none) \
            shared( latitude, longitude, wj, interpolated_fields, \
                    coastal_water_indices, coastal_land_indices, is_coastal_water ) \
            private( i_land, i_water, lat_land, lon_land, lat_water, lon_water, kern_val, \
                     Itime, Idepth, Ilat, Ilon, LAT, LON, curr_lat, curr_lon, \
                     LAT_ub, LAT_lb, LON_ub, LON_lb, index, inner_index ) \
            firstprivate( num_coastal_water, num_coastal_land, interp_scale, Ivar,\
                          Ntime, Ndepth, Nlat, Nlon )
            {
                #pragma omp for collapse(1) schedule(static)
                for ( i_land = 0; i_land < num_coastal_land; i_land++ ) {
                    index = coastal_land_indices[i_land];
                    Index1to4( index, Itime, Idepth, Ilat, Ilon,
                                      Ntime, Ndepth, Nlat, Nlon );
                    get_lat_bounds(LAT_lb, LAT_ub, latitude,  Ilat, interp_scale );

                    lon_land = longitude[Ilon];
                    lat_land = latitude[Ilat];

                    for ( LAT = LAT_lb; LAT < LAT_ub; LAT++) {

                        // Handle periodicity if necessary
                        if (constants::PERIODIC_Y) { curr_lat = ( LAT % Nlat + Nlat ) % Nlat; }
                        else                       { curr_lat = LAT; }
                        lat_water = latitude.at(curr_lat);

                        get_lon_bounds(LON_lb, LON_ub, longitude, Ilon, lat_land, lat_water, interp_scale );
                        for ( LON = LON_lb; LON < LON_ub; LON++ ) {

                            // Handle periodicity if necessary
                            if (constants::PERIODIC_X) { curr_lon = ( LON % Nlon + Nlon ) % Nlon; }
                            else                       { curr_lon = LON; }
                            lon_water = longitude[curr_lon];

                            inner_index = Index(Itime, Idepth, curr_lat, curr_lon, Ntime, Ndepth, Nlat, Nlon);

                            // If inner_index isn't on the coast, skip it
                            if (not(is_coastal_water[inner_index])) { continue; }

                            i_water = std::lower_bound( coastal_water_indices.begin(), 
                                                        coastal_water_indices.end(), 
                                                        inner_index ) 
                                - coastal_water_indices.begin();
    
                            kern_val = kernel( distance( lon_land,  lat_land, 
                                                         lon_water, lat_water ),
                                               interp_scale );
                            /*
                            kern_val = ( distance( lon_land,  lat_land, lon_water, lat_water ) < interp_scale / 2 )
                                        ? 1.
                                        : 0.;
                             */
    
                            interpolated_fields[Ivar][coastal_land_indices[i_land]] 
                                += wj[i_water] * kern_val;

                        }
                    }
                }
            }

            if ( interp_layer < num_interp_layers - 1 ) {
                // If this is the last interpolation layer, then no need to spend
                // time computing the residual, since no-one is going to use it
                fprintf(stdout, "  Computing residual.\n");
                double lat_water_ALL, lon_water_ALL;
                size_t i_water_ALL, j_water;
                // Get the residual of the fields for the next interpolation layer
                //  and apply the interpolator to coastal water
                #pragma omp parallel \
                default(none) \
                shared( latitude, longitude, wj, interpolated_fields, source_data, vars_to_interpolate, \
                        Aij, coastal_water_indices, coastal_water_indices_ALL, is_coastal_water ) \
                private( i_water_ALL, i_water, j_water, lat_water_ALL, lon_water_ALL, lat_water, lon_water, kern_val, \
                         Itime, Idepth, Ilat, Ilon, LAT, LON, curr_lat, curr_lon, \
                         LAT_ub, LAT_lb, LON_ub, LON_lb, index, inner_index ) \
                firstprivate( num_coastal_water, num_coastal_water_ALL, interp_scale, Ivar,\
                              Ntime, Ndepth, Nlat, Nlon )
                {
                    #pragma omp for collapse(1) schedule(static)
                    for ( i_water_ALL = 0; i_water_ALL < num_coastal_water_ALL; i_water_ALL++ ) {
                        index = coastal_water_indices_ALL[i_water_ALL];
                        Index1to4( index, Itime, Idepth, Ilat, Ilon,
                                          Ntime, Ndepth, Nlat, Nlon );
                        get_lat_bounds(LAT_lb, LAT_ub, latitude,  Ilat, interp_scale );
    
                        lon_water_ALL = longitude[Ilon];
                        lat_water_ALL = latitude[Ilat];
    
                        for ( LAT = LAT_lb; LAT < LAT_ub; LAT++) {
    
                            // Handle periodicity if necessary
                            if (constants::PERIODIC_Y) { curr_lat = ( LAT % Nlat + Nlat ) % Nlat; }
                            else                       { curr_lat = LAT; }
                            lat_water = latitude.at(curr_lat);
    
                            get_lon_bounds(LON_lb, LON_ub, longitude, Ilon, lat_water_ALL, lat_water, interp_scale );
                            for ( LON = LON_lb; LON < LON_ub; LON++ ) {
    
                                // Handle periodicity if necessary
                                if (constants::PERIODIC_X) { curr_lon = ( LON % Nlon + Nlon ) % Nlon; }
                                else                       { curr_lon = LON; }
                                lon_water = longitude[curr_lon];
    
                                inner_index = Index(Itime, Idepth, curr_lat, curr_lon, Ntime, Ndepth, Nlat, Nlon);
    
                                // If inner_index isn't on the coast, skip it
                                if (not(is_coastal_water[inner_index])) { continue; }
    
                                i_water = std::lower_bound( coastal_water_indices.begin(), 
                                                            coastal_water_indices.end(), 
                                                            inner_index ) 
                                    - coastal_water_indices.begin();
        
                                if ( is_coastal_water[index] ) {
                                    j_water = std::lower_bound( coastal_water_indices.begin(), 
                                                                coastal_water_indices.end(), 
                                                                index ) 
                                        - coastal_water_indices.begin();
                                    kern_val = Aij.coeff( i_water, j_water );
                                } else {
                                    kern_val = kernel( distance( lon_water_ALL,  lat_water_ALL, 
                                                                 lon_water,      lat_water      ),
                                                       interp_scale );
                                }
                                /*
                                kern_val = ( distance( lon_water_ALL, lat_water_ALL, lon_water, lat_water ) < interp_scale / 2 )
                                        ? 1.
                                            : 0.;
                                 */
        
                                interpolated_fields[Ivar][coastal_water_indices_ALL[i_water_ALL]] 
                                    += wj[i_water] * kern_val;
    
                                // Subtract off of the original data to leave residual
                                source_data.variables.at(vars_to_interpolate[Ivar]
                                        )[coastal_water_indices_ALL[i_water_ALL]]
                                    -= wj[i_water] * kern_val;

                            }
                        }
                    }
                }
            }
        } // end of Ivar loop

    } // End of interp layers loop


    //
    //// Create output file
    //
    if (not(constants::NO_FULL_OUTPUTS)) {
        initialize_output_file( source_data, vars_to_interpolate, output_fname, -1 );

        add_attr_to_file( "max_interp_scale",   max_interp_scale,   output_fname );
        add_attr_to_file( "num_interp_layers",  num_interp_layers,  output_fname );
        add_attr_to_file( "max_iterations",     max_iterations,     output_fname );
        add_attr_to_file( "tolerance",          tolerance,          output_fname );
        add_attr_to_file( "base_skip_factor",   base_skip_factor,   output_fname );

        const int ndims = 4;
        size_t starts[ndims];
        starts[0] = size_t(myStarts.at(0));
        starts[1] = size_t(myStarts.at(1));
        starts[2] = size_t(myStarts.at(2));
        starts[3] = size_t(myStarts.at(3));
        size_t counts[ndims] = { size_t(Ntime), size_t(Ndepth), size_t(Nlat), size_t(Nlon) };

        for ( Ivar = 0; Ivar < Nvars; Ivar++ ) {
            write_field_to_output( interpolated_fields.at(Ivar), 
                                   vars_to_interpolate.at(Ivar), 
                                   starts, counts, output_fname, NULL );
        }
    }


    // DONE!
    #if DEBUG >= 1
    fprintf(stdout, "Processor %d / %d waiting to finalize.\n", wRank + 1, wSize);
    #endif
    MPI_Finalize();

    fprintf( stdout, "\n\nProcess Complete\n\n" );
    return 0;

}
