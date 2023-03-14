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

//#include "../ALGLIB/stdafx.h"
//#include "../ALGLIB/interpolation.h"

#include "../netcdf_io.hpp"
#include "../functions.hpp"
#include "../constants.hpp"
#include "../preprocess.hpp"
#include "../differentiation_tools.hpp"

// Pragma-magic to allow reduction over vector operator
//      thanks to: https://stackoverflow.com/questions/43168661/openmp-and-reduction-on-stdvector
#pragma omp declare reduction(vec_double_plus : std::vector<double> : \
                              std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<double>())) \
                    initializer(omp_priv = decltype(omp_orig)(omp_orig.size()))


int main(int argc, char *argv[]) {
    
    // PERIODIC_Y implies UNIFORM_LAT_GRID
    static_assert( (constants::UNIFORM_LAT_GRID) or (not(constants::PERIODIC_Y)),
            "PERIODIC_Y requires UNIFORM_LAT_GRID.\n"
            "Please update constants.hpp accordingly.\n");
    static_assert( not(constants::CARTESIAN),
            "Toroidal projection now set to handle Cartesian coordinates.\n"
            );

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

    // first argument is the flag, second argument is default value (for when flag is not present)
    const std::string   &input_fname    = input.getCmdOption("--input_file",    "input.nc"),
                        &output_fname   = input.getCmdOption("--output_file",   "output.nc");

    const std::string   &time_dim_name      = input.getCmdOption("--time",        "time"),
                        &depth_dim_name     = input.getCmdOption("--depth",       "depth"),
                        &latitude_dim_name  = input.getCmdOption("--latitude",    "latitude"),
                        &longitude_dim_name = input.getCmdOption("--longitude",   "longitude");

    const std::string &latlon_in_degrees  = input.getCmdOption("--is_degrees",   "true");

    const std::string   &Nprocs_in_time_string  = input.getCmdOption("--Nprocs_in_time",  "1"),
                        &Nprocs_in_depth_string = input.getCmdOption("--Nprocs_in_depth", "1");
    const int   Nprocs_in_time_input  = stoi(Nprocs_in_time_string),
                Nprocs_in_depth_input = stoi(Nprocs_in_depth_string);

    std::vector< std::string > vars_to_interp, vars_in_output;
    input.getListofStrings( vars_to_interp, "--input_variables" );
    input.getListofStrings( vars_in_output, "--output_variables" );
    const int Nvars = vars_to_interp.size();

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

    // Apply some cleaning to the processor allotments if necessary. 
    source_data.Nlat = 1;
    source_data.Nlon = 1;
    source_data.check_processor_divisions( Nprocs_in_time_input, Nprocs_in_depth_input );
     
    // Convert to radians, if appropriate
    if ( latlon_in_degrees == "true" ) {
        convert_coordinates( source_data.longitude, source_data.latitude );
    }

    // Build the adjacency matrix and other adjacency-adjacent arrays
    source_data.build_adjacency();

    // Write the adjacency
    source_data.write_adjacency( "adjacency.nc" );

    // Read in field to get size info
    source_data.load_variable( "to_interp", vars_to_interp.at(0), input_fname, true, true );

    // Initialize file and write out coarsened fields
    if (wRank == 0) { fprintf( stdout, "Preparing output file\n" ); }
    initialize_output_file( source_data, vars_in_output, output_fname.c_str() );

    // Get the relevant sizes
    const size_t    Npts  = source_data.latitude.size();
    size_t starts[4] = { source_data.myStarts.at(0), source_data.myStarts.at(1), 0    };
    size_t counts[4] = { source_data.myCounts.at(0), source_data.myCounts.at(1), Npts };

    const size_t num_neighbours = source_data.num_neighbours;

    std::vector<double> lat_deriv(Npts, 0.), lon_deriv(Npts, 0.);

    // Next, the coarse velocities
    for ( int Ivar = 0; Ivar < Nvars; Ivar++ ) {

        // Load the data
        source_data.load_variable( "to_interp", vars_to_interp.at(Ivar), input_fname, true, true );
        std::vector<const std::vector<double>*> deriv_fields {&source_data.variables["to_interp"]};

        double lon_deriv_val, lat_deriv_val;

        std::vector<double*>    lon_deriv_vals {&lon_deriv_val},
                                lat_deriv_vals {&lat_deriv_val};

        size_t index, II;
        double weight, value;
        #pragma omp parallel default(none)\
        private( index, value, weight, II, \
                lon_deriv_val, lat_deriv_val, lon_deriv_vals, lat_deriv_vals ) \
        shared( source_data, deriv_fields ) \
        firstprivate( Npts, num_neighbours, stdout )\
        reduction( vec_double_plus : lat_deriv, lon_deriv )
        {
            lon_deriv_vals.clear();
            lon_deriv_vals.push_back(&lon_deriv_val);

            lat_deriv_vals.clear();
            lat_deriv_vals.push_back(&lat_deriv_val);

            #pragma omp for collapse(1) schedule(static)
            for ( index = 0; index < Npts; index++ ) {

                //fprintf( stdout, "%'zu - start lon deriv\n", index );
                spher_derivative_at_point(
                        lon_deriv_vals, deriv_fields, source_data.longitude, "lon", source_data,
                        0, 0, index, index, source_data.mask
                        );
                //fprintf( stdout, "%'zu - start lat deriv\n", index );
                spher_derivative_at_point(
                        lat_deriv_vals, deriv_fields, source_data.latitude, "lat", source_data,
                        0, 0, index, index, source_data.mask
                        );
                //fprintf( stdout, "%'zu - done derivs\n", index );
                lon_deriv.at(index) = lon_deriv_val;
                lat_deriv.at(index) = lat_deriv_val;

                /*
                for ( II = 0; II < num_neighbours+1; II++ ) {
                    if (II < num_neighbours) {
                        value = source_data.variables["to_interp"].at(source_data.adjacency_indices.at(index).at(II));
                    } else {
                        value = source_data.variables["to_interp"].at(index);
                    }

                    weight = source_data.adjacency_ddlat_weights.at(index).at(II);
                    lat_deriv.at(index) += value * weight;

                    weight = source_data.adjacency_ddlon_weights.at(index).at(II);
                    lon_deriv.at(index) += value * weight;
                }
                */
            }

        }

        // Write the data
        write_field_to_output( lat_deriv, "ddlat_" + vars_to_interp.at(Ivar), 
                starts, counts, output_fname );
        write_field_to_output( lon_deriv, "ddlon_" + vars_to_interp.at(Ivar), 
                starts, counts, output_fname );
    }

    #if DEBUG >= 1
    fprintf(stdout, "Processor %d / %d waiting to finalize.\n", wRank + 1, wSize);
    #endif
    MPI_Finalize();
    return 0;
}
