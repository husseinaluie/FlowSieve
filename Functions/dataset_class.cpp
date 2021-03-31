#include <math.h>
#include <algorithm>
#include "../constants.hpp"
#include "../functions.hpp"
#include "../netcdf_io.hpp"
#include <cassert>


// Class constructor
dataset::dataset() {

};

void dataset::load_time( const std::string dim_name, const std::string filename ) {
    read_var_from_file(time, dim_name, filename);
    full_Ntime = time.size();
};

void dataset::load_depth( const std::string dim_name, const std::string filename ) {
    read_var_from_file(depth, dim_name, filename);
    full_Ndepth = depth.size();
};

void dataset::load_latitude( const std::string dim_name, const std::string filename ) {
    read_var_from_file(latitude, dim_name, filename);
    Nlat = latitude.size();
};

void dataset::load_longitude( const std::string dim_name, const std::string filename ) {
    read_var_from_file(longitude, dim_name, filename);
    Nlon = longitude.size();
};

void dataset::compute_cell_areas() {

    assert( (Nlat > 0) and (Nlon > 0) );    

    areas.resize( Nlat * Nlon );
    compute_areas( areas, longitude, latitude );
}

void dataset::load_variable( 
        const std::string var_name, 
        const std::string var_name_in_file, 
        const std::string filename,
        const bool read_mask,
        const bool load_counts
        ) {

    // Add a new entry to the variables dictionary with an empty array
    variables.insert( std::pair< std::string, std::vector<double> >( var_name, std::vector<double>() ) );

    // Now read in from the file and store in the variables dictionary
    read_var_from_file( variables.at( var_name ), var_name_in_file, filename, 
                        read_mask ? &mask : NULL, 
                        load_counts ? &myCounts : NULL, 
                        load_counts ? &myStarts : NULL, 
                        Nprocs_in_time, Nprocs_in_depth );
};

void dataset::check_processor_divisions( const int Nprocs_in_time_input, const int Nprocs_in_depth_input, const MPI_Comm ) {

    assert( (full_Ntime > 0) and (full_Ndepth > 0) and (Nlon > 0) and (Nlat > 0) ); // Must read in dimensions before checking processor divisions.

    int wRank=-1, wSize=-1;
    MPI_Comm_rank( MPI_COMM_WORLD, &wRank );
    MPI_Comm_size( MPI_COMM_WORLD, &wSize );

    // Apply some cleaning to the processor allotments if necessary. 
    Nprocs_in_time  = ( Ntime  == 1 ) ? 1 : 
                      ( Ndepth == 1 ) ? wSize : 
                                        Nprocs_in_time_input;
    Nprocs_in_depth = ( Ndepth == 1 ) ? 1 : 
                      ( Ntime  == 1 ) ? wSize : 
                                        Nprocs_in_depth_input;

    #if DEBUG >= 0
    if (Nprocs_in_time != Nprocs_in_time_input) { 
        if (wRank == 0) { fprintf(stdout, " WARNING!! Changing number of processors in time to %'d from %'d\n", Nprocs_in_time, Nprocs_in_time_input); }
    }
    if (Nprocs_in_depth != Nprocs_in_depth_input) { 
        if (wRank == 0) { fprintf(stdout, " WARNING!! Changing number of processors in depth to %'d from %'d\n", Nprocs_in_depth, Nprocs_in_depth_input); }
    }
    if (wRank == 0) { fprintf(stdout, " Nproc(time, depth) = (%'d, %'d)\n\n", Nprocs_in_time, Nprocs_in_depth); }
    #endif

    assert( Nprocs_in_time * Nprocs_in_depth == wSize );
}
