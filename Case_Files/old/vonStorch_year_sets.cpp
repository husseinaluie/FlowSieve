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
#include <bitset>

#include "../netcdf_io.hpp"
#include "../functions.hpp"
#include "../differentiation_tools.hpp"
#include "../constants.hpp"
#include "../postprocess.hpp"

size_t factorial( size_t n ) {
    return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

void merge_time_averages(   
            std::map< std::string, std::vector<double> > & time_means,
            const std::vector< std::string > & filenames,
            const size_t Npts 
        ){

    double Ntime = 0, curr_Ntime;

    std::vector<double> storage;

    std::vector< std::string > list_of_vars { "u", "v", "uu", "uv", "vv" };

    // First, zero out the time_means
    for ( size_t Ivar = 0; Ivar < list_of_vars.size(); Ivar++ ) {
        std::fill( time_means[list_of_vars.at(Ivar)].begin(), 
                   time_means[list_of_vars.at(Ivar)].end(), 
                   0.);
    }

    size_t index;
    for ( size_t Ifile = 0; Ifile < filenames.size(); Ifile++ ) {

        // Get number of time points from current average
        read_attr_from_file( curr_Ntime, "Ntime_used_in_average", filenames.at(Ifile) );
        assert( curr_Ntime >= 1 );
        Ntime += curr_Ntime;

        for ( size_t Ivar = 0; Ivar < list_of_vars.size(); Ivar++ ) {

            // Read in the time average
            read_var_from_file( storage, list_of_vars.at(Ivar), filenames.at(Ifile) );

            // Add on to the stored time_mean
            #pragma omp parallel default(none) shared( time_means, list_of_vars, storage, curr_Ntime, Ivar ) private( index )
            {
                #pragma omp for collapse(1) schedule(static)
                for ( index = 0; index < Npts; index++ ) {
                    time_means[ list_of_vars.at(Ivar) ].at(index) += curr_Ntime * storage.at(index);
                }
            }
        }
    }

    // Finally, divide out by the total accumulated number of days
    for ( size_t Ivar = 0; Ivar < list_of_vars.size(); Ivar++ ) {
        #pragma omp parallel default(none) shared( time_means, list_of_vars, Ntime, Ivar ) private( index )
        {
            #pragma omp for collapse(1) schedule(static)
            for ( index = 0; index < Npts; index++ ) {
                time_means[ list_of_vars.at(Ivar) ].at(index) *= 1. / Ntime;
            }
        }
    }

}


void compute_vonStorch(
        std::vector<double> & vonStorch,
        const std::map< std::string, std::vector<double> > & time_means,
        const dataset & source_data
        ) {

    // Some handy references
    const std::vector<double>   &u  = time_means.at("u"),
                                &v  = time_means.at("v"),
                                &uu = time_means.at("uu"),
                                &uv = time_means.at("uv"),
                                &vu = time_means.at("uv"),
                                &vv = time_means.at("vv");

    const std::vector<double>   &latitude  = source_data.latitude,
                                &longitude = source_data.longitude;

    const int   Ntime   = 1,
                Ndepth  = source_data.Ndepth,
                Nlat    = source_data.Nlat,
                Nlon    = source_data.Nlon;

    const size_t Npts = vonStorch.size();

    // Compute vonStorch point-wise
    std::fill( vonStorch.begin(), vonStorch.end(), 0.);
    std::vector<const std::vector<double>*> deriv_fields { &u, &v };

    double dudlon, dudlat, dvdlon, dvdlat;
    std::vector<double*> lon_deriv_vals, lat_deriv_vals;

    size_t index;
    int Itime, Idepth, Ilat, Ilon;

    const int OMP_chunksize = get_omp_chunksize( Nlat, Nlon );

    double u_loc, v_loc, tau_uu, tau_uv, tau_vu, tau_vv;
    #pragma omp parallel \
    default(none) \
    shared( source_data, latitude, longitude, deriv_fields, u, v, uu, uv, vu, vv, vonStorch )\
    private( Itime, Idepth, Ilat, Ilon, index, dudlon, dudlat, dvdlon, dvdlat, lon_deriv_vals, lat_deriv_vals,\
             u_loc, v_loc, tau_uu, tau_uv, tau_vu, tau_vv )
    {

        lon_deriv_vals.push_back( &dudlon );
        lon_deriv_vals.push_back( &dvdlon );

        lat_deriv_vals.push_back( &dudlat );
        lat_deriv_vals.push_back( &dvdlat );

        #pragma omp for collapse(1) schedule(dynamic, OMP_chunksize)
        for ( index = 0; index < Npts; index++ ) {

            if (source_data.mask.at(index)) {
                Index1to4( index, Itime, Idepth, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon );

                spher_derivative_at_point( lat_deriv_vals, deriv_fields, latitude,  "lat", Itime, Idepth, Ilat, Ilon, 
                                                                                           Ntime, Ndepth, Nlat, Nlon, source_data.mask );
                spher_derivative_at_point( lon_deriv_vals, deriv_fields, longitude, "lon", Itime, Idepth, Ilat, Ilon, 
                                                                                           Ntime, Ndepth, Nlat, Nlon, source_data.mask );

                // rho0 * [ tau( u, vec(u) ) dot grad( u )  +  tau( v, vec(u) ) dot grad ( v ) ]
                //      where tau( a, b ) = bar(ab) - bar(a) bar(b)
                //      and bar(.) is a time mean
                
                u_loc = u.at(index);
                v_loc = v.at(index);
                tau_uu = uu.at(index) - u_loc * u_loc;
                tau_uv = uv.at(index) - u_loc * v_loc;
                tau_vu = vu.at(index) - v_loc * u_loc;
                tau_vv = vv.at(index) - v_loc * v_loc;
                      
                vonStorch.at(index) = ( constants::rho0 / constants::R_earth ) * (
                              tau_uu * dudlon / cos(latitude.at(Ilat))
                            + tau_uv * dudlat 
                            + tau_vu * dvdlon / cos(latitude.at(Ilat))
                            + tau_vv * dvdlat
                        );
            }
        }
    }
}



int main(int argc, char *argv[]) {
    
    static_assert ( not(constants::FILTER_OVER_LAND), "Cannot have FILTER_OVER_LAND on when computing vonStorch" );

    static_assert( not( constants::DO_OKUBOWEISS_ANALYSIS ), "No OkuboWeiss available to pass to post-processing." );

    static_assert( not( constants::POSTPROCESS_DO_TIME_MEANS ), "year_sets does not handle spatial maps output" );

    // Enable all floating point exceptions but FE_INEXACT
    //feenableexcept( FE_ALL_EXCEPT & ~FE_INEXACT & ~FE_UNDERFLOW );
    //fprintf( stdout, " %d : %d \n", FE_ALL_EXCEPT, FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW | FE_INEXACT | FE_UNDERFLOW );
    feenableexcept( FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW );

    // Specify the number of OpenMP threads
    //   and initialize the MPI world
    int thread_safety_provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &thread_safety_provided);
    //MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI::ERRORS_THROW_EXCEPTIONS);
    const double start_time = MPI_Wtime();

    int wRank=-1, wSize=-1;
    MPI_Comm_rank( MPI_COMM_WORLD, &wRank );
    MPI_Comm_size( MPI_COMM_WORLD, &wSize );
    assert(wSize == 1);

    //
    //// Parse command-line arguments
    //
    InputParser input(argc, argv);
    if(input.cmdOptionExists("--version")){
        if (wRank == 0) { print_compile_info(NULL); } 
        return 0;
    }

    // first argument is the flag, second argument is default value (for when flag is not present)
    const std::string   &postproc_fname     = input.getCmdOption("--postproc_file", "postprocess");

    const std::string   &time_dim_name      = input.getCmdOption("--time",        "time"),
                        &depth_dim_name     = input.getCmdOption("--depth",       "depth"),
                        &latitude_dim_name  = input.getCmdOption("--latitude",    "latitude"),
                        &longitude_dim_name = input.getCmdOption("--longitude",   "longitude");

    const std::string &latlon_in_degrees  = input.getCmdOption("--is_degrees",   "true");

    const std::string   &Nprocs_in_time_string  = input.getCmdOption("--Nprocs_in_time",  "1"),
                        &Nprocs_in_depth_string = input.getCmdOption("--Nprocs_in_depth", "1");
    const int   Nprocs_in_time_input  = stoi(Nprocs_in_time_string),
                Nprocs_in_depth_input = stoi(Nprocs_in_depth_string);

    const std::string   &u_name     = input.getCmdOption("--u",     "u"),
                        &v_name     = input.getCmdOption("--v",     "v"),
                        &uu_name    = input.getCmdOption("--uu",    "uu"),
                        &uv_name    = input.getCmdOption("--uv",    "uv"),
                        &vv_name    = input.getCmdOption("--vv",    "vv");

    const std::string   &region_defs_fname    = input.getCmdOption("--region_definitions_file",    "region_definitions.nc"),
                        &region_defs_dim_name = input.getCmdOption("--region_definitions_dim",     "region"),
                        &region_defs_var_name = input.getCmdOption("--region_definitions_var",     "region_definition");

    // Read in the list of year files
    std::vector< std::string > list_of_year_files;
    input.getListofStrings( list_of_year_files, "--data_files" );
    const size_t Nyears = list_of_year_files.size();
    const size_t total_year_subset = pow(2, Nyears);

    // Get list of year-durations to sample
    std::vector<double> sample_durations;
    input.getFilterScales( sample_durations, "--sample_durations" );
    assert( sample_durations.size() > 0 );

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
    //   implicitely assume coordinates are the same between input files
    source_data.load_time(      time_dim_name,      list_of_year_files[0] );
    source_data.load_depth(     depth_dim_name,     list_of_year_files[0] );
    source_data.load_latitude(  latitude_dim_name,  list_of_year_files[0] );
    source_data.load_longitude( longitude_dim_name, list_of_year_files[0] );

    // Apply some cleaning to the processor allotments if necessary. 
    source_data.check_processor_divisions( Nprocs_in_time_input, Nprocs_in_depth_input );
     
    // Convert to radians, if appropriate
    if ( latlon_in_degrees == "true" ) { convert_coordinates( source_data.longitude, source_data.latitude ); }

    // Compute the area of each 'cell' which will be necessary for integration
    source_data.compute_cell_areas();

    // Read in one sample to get the grid and mask, then erase it to save memory
    source_data.load_variable( "temp", u_name, list_of_year_files[0], true, true );
    source_data.variables.at("temp").clear();

    const std::vector<double>   &latitude  = source_data.latitude,
                                &longitude = source_data.longitude;

    // Get the MPI-local dimension sizes
    source_data.Ntime  = source_data.myCounts[0];
    source_data.Ndepth = source_data.myCounts[1];
    const int   Stime   = source_data.myStarts.at(0),
                Sdepth  = source_data.myStarts.at(1);

    //
    const int   Ntime   = source_data.Ntime,
                Ndepth  = source_data.Ndepth,
                Nlat    = source_data.Nlat,
                Nlon    = source_data.Nlon;

    const size_t Npts = Ndepth * Nlat * Nlon;

    // Mask out the pole, if necessary (i.e. set lat = 90 to land)
    mask_out_pole( source_data.latitude, source_data.mask, Ntime, Ndepth, Nlat, Nlon );


    // Storage for the computed time averages
    std::map< std::string, std::vector<double> > subset_time_means = {
        { "u",  std::vector<double>(Npts, 0.) },
        { "v",  std::vector<double>(Npts, 0.) },
        { "uu", std::vector<double>(Npts, 0.) },
        { "uv", std::vector<double>(Npts, 0.) },
        { "vv", std::vector<double>(Npts, 0.) }
    };

    std::vector<double> subset_vonStorch( Npts, 0. );

    // Set up to pass vonStorch to area averaging
    std::vector<const std::vector<double>*> postprocess_fields;
    postprocess_fields.push_back( &subset_vonStorch );

    std::vector<std::string> postprocess_names;
    postprocess_names.push_back( "vonStorch" );

    // Storage for the horizontally-averaged values
    std::vector< std::vector< double > >    field_averages( 1 ), 
                                            field_std_devs( 1 ),
                                            subset_field_averages( 1, std::vector<double>( Ndepth, 0. )),
                                            subset_field_std_devs( 1, std::vector<double>( Ndepth, 0. ));
    std::vector< double > Okubo_placeholder;

    // Array to hold the subset of file names
    std::vector< std::string > file_subset;

    // Otherwise, just make a single region which is the entire domain
    source_data.region_names.push_back("full_domain");
    source_data.regions.insert( std::pair< std::string, std::vector<bool> >( 
                "full_domain", std::vector<bool>( source_data.Nlat * source_data.Nlon, true) ) 
            );
    source_data.compute_region_areas();


    // Loop through all samples sizes ranging from 1 year to Nyears (maximum)
    for ( size_t duration_index = 0; duration_index <= sample_durations.size(); duration_index++ ) {

        size_t years_in_sample = sample_durations.at(duration_index);

        size_t sample_count = 0;

        // Resize the field_averages appropriately
        size_t num_samples_of_size = factorial( Nyears ) / ( factorial( years_in_sample ) * factorial( Nyears - years_in_sample ) );
        field_averages[0].resize( num_samples_of_size * Ndepth );
        field_std_devs[0].resize( num_samples_of_size * Ndepth );

        // Loop through every element in the power set of the available years
        for ( size_t sample_index = 1; sample_index < total_year_subset; sample_index++ ) {

            #if DEBUG >= 1
            fprintf( stdout, "  Testing sample index %'zu of %'zu\n", sample_index, total_year_subset-1 );
            fflush( stdout );
            #endif

            std::bitset<32> sample_expansion( sample_index );

            #if DEBUG >= 1
            std::string mystring = sample_expansion.to_string<char, std::string::traits_type, std::string::allocator_type>();
            fprintf( stdout, "    has bitset %s\n", mystring.c_str() );
            fflush( stdout );
            #endif
            
            // Check if this element in the power set has the right number of years
            if ( sample_expansion.count() == years_in_sample ) {

                #if DEBUG >= 1
                fprintf( stdout, "    Processing bitset...\n" );
                #endif

                // Then make a list of the selected year samples
                file_subset.clear();
                for( size_t Iyear = 0; Iyear < Nyears; Iyear++ ) {
                    if ( sample_expansion.test(Iyear) ) {
                        file_subset.push_back( list_of_year_files.at(Iyear) );
                    }
                }

                // Then, get the time averages corresponding to that year subset
                merge_time_averages( subset_time_means, file_subset, Npts );

                // Next, compute the vonStorch values for that year subset
                compute_vonStorch( subset_vonStorch, subset_time_means, source_data );

                // Then apply the area averaging routines
                compute_region_avg_and_std( subset_field_averages, subset_field_std_devs, source_data, postprocess_fields );

                // And now copy those area average results into the array for the whole set
                for ( int Idepth = 0; Idepth < Ndepth; Idepth++ ) {
                    field_averages[0].at( Idepth + sample_count * Ndepth ) = subset_field_averages[0].at( Idepth );
                }

                // And finally increment our total sample count
                sample_count++;
                #if DEBUG >= 0
                fprintf( stdout, "    Finished sample %'zu of %'zu\n", sample_count, num_samples_of_size );
                #endif

            }
        }

        fprintf( stdout, "Extracted %'zu subsets (expects %'zu) of size %'zu from the provided file list.\n", 
                sample_count, num_samples_of_size, years_in_sample );
        fflush( stdout );

        // And now write those horizontally-averaged values to file
        char postproc_output_filename[50];
        snprintf( postproc_output_filename, 50, (postproc_fname + "_%zu.nc").c_str(), years_in_sample);

        // Update time and region areas before outputting
        source_data.Ntime = num_samples_of_size;
        source_data.time.resize( num_samples_of_size );
        for ( size_t Itime = 0; Itime < num_samples_of_size; Itime++ ) { source_data.time[Itime] = Itime; }

        source_data.region_areas.resize( num_samples_of_size * Ndepth );
        for ( size_t Idepth = 0; Idepth < Ndepth; Idepth++ ) {
            for ( size_t Itime = 1; Itime < num_samples_of_size; Itime++ ) { 
                source_data.region_areas.at( Idepth + Itime * Ndepth ) = source_data.region_areas.at( Idepth );
            }
        }

        // Set up output file
        initialize_postprocess_file( source_data, Okubo_placeholder, postprocess_names, postproc_output_filename, -1, false );

        // Put time back to what it was
        source_data.Ntime = 1;
        source_data.time.resize( 1 );
        source_data.time[0] = 0.;
        source_data.region_areas.resize( Ndepth );

        write_region_avg_and_std( field_averages, field_std_devs, postprocess_names, postproc_output_filename, 
                Stime, Sdepth, num_samples_of_size, Ndepth, 1, 1 );

    } 


    // Done!
    MPI_Finalize();
    return 0;
}
