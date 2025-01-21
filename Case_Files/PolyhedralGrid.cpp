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
#include "../preprocess.hpp"
#include "../differentiation_tools.hpp"

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
    const bool asked_help = input.cmdOptionExists("--help");
    if (asked_help) {
        fprintf( stdout, "\033[1;4mThe command-line input arguments [and default values] are:\033[0m\n" );
    }

    // first argument is the flag, second argument is default value (for when flag is not present)
    const std::string   &output_fname   = input.getCmdOption("--output_file",   "output.nc");
    const std::string   &adj_out_fname  = input.getCmdOption("--adjacency_file", "adjacency.nc");

    const std::string &Npts_string = input.getCmdOption("--target_grid_size", 
                                                        "100000", 
                                                        asked_help,
                                                        "Desired number of cells on the computational grid. NOTE: Output will be close, but may not have the exact size requested.");
    const size_t Npts_requested = (size_t) std::stoll(Npts_string);  
    if (asked_help) { return 0; }

    // Print processor assignments
    const int max_threads = omp_get_max_threads();
    omp_set_num_threads( max_threads );

    // Print some header info, depending on debug level
    print_header_info();

    // Initialize dataset class instance
    dataset polyhedral_grid;
    polyhedral_grid.time.resize(1, 0);
    polyhedral_grid.depth.resize(1, 0);

    // Apply some cleaning to the processor allotments if necessary. 
    polyhedral_grid.Ntime = 1;
    polyhedral_grid.Ndepth = 1;
    polyhedral_grid.Nlat = 1;
    polyhedral_grid.Nlon = 1;

    // Build the grid
    BuildPolyhedralGrid( polyhedral_grid, Npts_requested );
    const size_t Npts = polyhedral_grid.longitude.size();
    polyhedral_grid.mask.resize( Npts, true );
    polyhedral_grid.write_adjacency( adj_out_fname );

    // Initialize file and write out coarsened fields
    if (wRank == 0) { fprintf( stdout, "Preparing output file\n" ); }
    std::vector< std::string > vars_in_output;
    vars_in_output.push_back("cell_area");
    vars_in_output.push_back("scalar");
    vars_in_output.push_back("ddlon_scalar");
    vars_in_output.push_back("ddlat_scalar");
    initialize_output_file( polyhedral_grid, vars_in_output, output_fname );

    // Get the relevant sizes
    size_t starts[3] = { 0, 0, 0    };    
    size_t counts[3] = { 1, 1, Npts };

    std::vector<double> scalar(Npts, 0.);
    for ( size_t Ipt = 0; Ipt < Npts; Ipt++ ) {

        double lon = polyhedral_grid.longitude[Ipt],
               lat = polyhedral_grid.latitude[Ipt];

        scalar[Ipt] = cos( lon ) * exp( -pow( lat / (M_PI/6) , 2.) );
    }

    std::vector<double> ddlon_scalar(Npts, 0.), ddlat_scalar(Npts, 0.);
    std::vector< const std::vector<double>* > deriv_fields;
    deriv_fields.push_back( &scalar );

    double dfdlon, dfdlat;
    std::vector< double* > lon_deriv_vals, lat_deriv_vals;
    lon_deriv_vals.push_back( &dfdlon );
    lat_deriv_vals.push_back( &dfdlat );
    for ( size_t Ipt = 0; Ipt < Npts; Ipt++ ) {
        spher_derivative_at_point(
                lon_deriv_vals, deriv_fields, polyhedral_grid.longitude, "lon",
                polyhedral_grid, 0, 0, Ipt, Ipt, polyhedral_grid.mask);
        ddlon_scalar[Ipt] = dfdlon;

        spher_derivative_at_point(
                lat_deriv_vals, deriv_fields, polyhedral_grid.latitude, "lat",
                polyhedral_grid, 0, 0, Ipt, Ipt, polyhedral_grid.mask);
        ddlat_scalar[Ipt] = dfdlat;
    }

    // Write the data
    write_field_to_output( polyhedral_grid.areas, "cell_area", starts, counts, output_fname );
    write_field_to_output( scalar, "scalar", starts, counts, output_fname );
    write_field_to_output( ddlon_scalar, "ddlon_scalar", starts, counts, output_fname );
    write_field_to_output( ddlat_scalar, "ddlat_scalar", starts, counts, output_fname );

    #if DEBUG >= 1
    fprintf(stdout, "Processor %d / %d waiting to finalize.\n", wRank + 1, wSize);
    #endif
    MPI_Finalize();
    return 0;
}
