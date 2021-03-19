#include <math.h>
#include <mpi.h>
#include <omp.h>
#include <vector>
#include <cassert>

#include "../constants.hpp"
#include "../functions.hpp"
#include "../postprocess.hpp"
#include "../netcdf_io.hpp"

void Apply_Postprocess_Routines(
        const std::vector<const std::vector<double>*> & postprocess_fields,
        const std::vector<std::string> & vars_to_process,
        const std::vector<double> & OkuboWeiss,
        const std::vector<double> & time,
        const std::vector<double> & depth,
        const std::vector<double> & latitude,
        const std::vector<double> & longitude,
        const std::vector<bool>   & mask,
        const std::vector<double> & areas,
        const std::vector<int>    & myCounts,
        const std::vector<int>    & myStarts,
        const double filter_scale,
        const std::string filename_base,
        const MPI_Comm comm
        ) {

    //static_assert( not(constants::APPLY_POSTPROCESS) or not(constants::CAST_TO_INT) );

    const int Ntime  = myCounts.at(0);
    const int Ndepth = myCounts.at(1);
    const int Nlat   = myCounts.at(2);
    const int Nlon   = myCounts.at(3);

    // Get full number of time points
    int full_Ntime;
    MPI_Allreduce(&Ntime, &full_Ntime, 1, MPI_INT, MPI_SUM, comm);

    const int Stime  = myStarts.at(0);
    const int Sdepth = myStarts.at(1);

    int wRank=-1, wSize=-1;
    MPI_Comm_rank( MPI_COMM_WORLD, &wRank );
    MPI_Comm_size( MPI_COMM_WORLD, &wSize );

    //double dA, local_area;
    int Ilat, Ilon, Itime, Idepth;
    size_t index, area_index;

    // Extract a common mask that determines what points are always masked.
    //    Also keep a tally of how often a cell is masked
    std::vector<bool> always_masked(   Ndepth * Nlat * Nlon, true );
    std::vector<int>  mask_count(      Ndepth * Nlat * Nlon, 0 ),
                      mask_count_loc(  Ndepth * Nlat * Nlon, 0 );

    #pragma omp parallel default(none) \
    private( index, area_index, Itime, Idepth, Ilat, Ilon ) \
    shared( mask_count_loc, mask )
    { 
        #pragma omp for collapse(1) schedule(static)
        for (index = 0; index < mask.size(); ++index) {
            Index1to4(index, Itime, Idepth, Ilat, Ilon,
                             Ntime, Ndepth, Nlat, Nlon);

            area_index = Index(0, Idepth, Ilat, Ilon,
                               1, Ndepth, Nlat, Nlon);

            // Add up the number of times a cell is water (not masked)
            if ( mask.at(index) ) { mask_count_loc.at(area_index) = mask_count_loc.at(area_index) + 1; }
        }
    }
    MPI_Allreduce( &(mask_count_loc[0]),
                   &(mask_count[0]),
                   Ndepth * Nlat * Nlon, MPI_INT, MPI_SUM, comm);

    #pragma omp parallel default(none) \
    private( index ) shared( mask_count, always_masked )
    { 
        #pragma omp for collapse(1) schedule(static)
        for (index = 0; index < mask_count.size(); ++index) {
            always_masked.at(index) = mask_count.at(index) == 0 ? true : false;
        }
    }

    //
    //// Write a file that defines the regions
    //
    char filename[50];
    snprintf(filename, 50, "postprocess_regions.nc");
    write_regions(filename, latitude, longitude, mask, areas, myCounts, myStarts);
    write_regions_to_post(filename);


    //
    //// Start by initializing the postprocess file
    //
   
    // If we're using OkuboWeiss, create dimension values
    const bool do_OkuboWeiss = constants::DO_OKUBOWEISS_ANALYSIS;
    const int N_Okubo = 601; // Make it odd
    const int N_Ok_by_2 = N_Okubo / 2;
    std::vector<double> OkuboWeiss_dim_vals( do_OkuboWeiss ? N_Okubo : 0 );
    const double Okubo_max = 1.e10 / pow( 24. * 60. * 60., 2.);
    if (do_OkuboWeiss) {

        // This does log spacing (except for 0)
        const double decades = 30,
                     log_UB = log10( Okubo_max ),
                     log_LB = log_UB - decades;
        OkuboWeiss_dim_vals.at( N_Ok_by_2 ) = 0.;
        for (int II = 1; II <= N_Ok_by_2; ++II) {
            OkuboWeiss_dim_vals.at(N_Ok_by_2 + II) =  pow( 10., log_LB + II * ( decades / N_Ok_by_2 ) );
            OkuboWeiss_dim_vals.at(N_Ok_by_2 - II) = -OkuboWeiss_dim_vals.at( N_Ok_by_2 + II );
        }

    }

    // Filename + file
    snprintf(filename, 50, (filename_base + "_%.6gkm.nc").c_str(), filter_scale/1e3);
    initialize_postprocess_file(
            time, depth, latitude, longitude, OkuboWeiss_dim_vals,
            RegionTest::region_names,
            vars_to_process,
            filename,
            filter_scale,
            do_OkuboWeiss
            );

    // Add some attributes to the file
    const double kern_alpha = kernel_alpha();
    add_attr_to_file("kernel_alpha", 
            kern_alpha * pow(filter_scale, 2), 
            filename);

    //
    //// Compute the area of each region
    //
    const int num_fields  = vars_to_process.size();
    const int num_regions = RegionTest::all_regions.size();

    if (wRank == 0) { fprintf(stdout, "  .. computing the area of each region\n"); }
    fflush(stdout);

    std::vector<double> region_areas(num_regions * Ntime * Ndepth, 0.);
    compute_region_areas(region_areas, areas, mask, latitude, longitude,
            num_regions, Ntime, Ndepth, Nlat, Nlon);

    // Write the region areas
    size_t start_r[3], count_r[3];
    start_r[0] = (size_t) Stime;
    count_r[0] = (size_t) Ntime;

    start_r[1] = (size_t) Sdepth;
    count_r[1] = (size_t) Ndepth;

    start_r[2] = 0;
    count_r[2] = (size_t) num_regions;

    write_field_to_output(region_areas, "region_areas", start_r, count_r, filename, NULL);

    if (wRank == 0) { fprintf(stdout, "  .. .. wrote region areas\n"); }
    fflush(stdout);

    // Write region names
    //   this has to be done separately for reasons
    write_regions_to_post(filename);

    if (wRank == 0) { fprintf(stdout, "  .. .. wrote region names\n"); }
    fflush(stdout);


    //
    //// Region averages and standard deviations
    //

    if (wRank == 0) { fprintf(stdout, "  .. computing the average and std. dev. in each region\n"); }
    fflush(stdout);

    std::vector< std::vector< double > >
        field_averages(num_fields, std::vector<double>(Ntime * Ndepth * num_regions, 0.)), 
        field_std_devs(num_fields, std::vector<double>(Ntime * Ndepth * num_regions, 0.));

    compute_region_avg_and_std(
            field_averages, field_std_devs, postprocess_fields, region_areas,
            areas, mask, latitude, longitude,
            num_fields, num_regions, Ntime, Ndepth, Nlat, Nlon
            );

    if (wRank == 0) { fprintf(stdout, "  .. writing region averages and deviations\n"); }
    fflush(stdout);

    write_region_avg_and_std(
            field_averages, field_std_devs, vars_to_process, filename,
            Stime, Sdepth, Ntime, Ndepth, num_regions, num_fields
            );


    //
    //// If we have OkuboWeiss data, then also do processing along OW contours
    //

    if (wRank == 0) { fprintf(stdout, "  .. Applying Okubo-Weiss processing\n"); }
    fflush(stdout);

    std::vector< double > OkuboWeiss_areas;
    if (do_OkuboWeiss) {

        std::vector< std::vector< double > > 
            field_averages_OW(num_fields, std::vector<double>(Ntime * Ndepth * N_Okubo * num_regions, 0.)), 
            field_std_devs_OW(num_fields, std::vector<double>(Ntime * Ndepth * N_Okubo * num_regions, 0.));

        OkuboWeiss_areas.resize(Ntime * Ndepth * N_Okubo * num_regions, 0.);

        compute_region_avg_and_std_OkuboWeiss(
                field_averages_OW, field_std_devs_OW, OkuboWeiss_areas, 
                postprocess_fields, OkuboWeiss, region_areas,
                areas, mask, OkuboWeiss_dim_vals, latitude, longitude,
                num_fields, num_regions, Ntime, Ndepth, Nlat, Nlon, N_Okubo
                );

        if (wRank == 0) { fprintf(stdout, "  .. writing Okubo results\n"); }
        fflush(stdout);

        write_region_avg_and_std_OkuboWeiss(
                field_averages_OW, field_std_devs_OW, OkuboWeiss_areas,
                vars_to_process, filename,
                Stime, Sdepth, Ntime, Ndepth, N_Okubo, num_regions, num_fields
                );
    }

    //
    //// Time averages
    //

    if (wRank == 0) { fprintf(stdout, "  .. computing time-averages of fields\n"); }
    fflush(stdout);

    // Time-average the fields
    std::vector<std::vector<double>> time_average(num_fields), time_std_dev(num_fields);
    int Ifield;
    for (Ifield = 0; Ifield < num_fields; ++Ifield) {
        time_average.at( Ifield ).resize( Ndepth * Nlat * Nlon, 0. );
        time_std_dev.at( Ifield ).resize( Ndepth * Nlat * Nlon, 0. );
    }

    compute_time_avg_std(
            time_average, time_std_dev, postprocess_fields,
            mask, areas, latitude, longitude, mask_count, always_masked,
            full_Ntime, num_fields, Ntime, Ndepth, Nlat, Nlon
            );

    // Write the time averages
    //   dimension order: depth - lat - lon
    const int Slat   = myStarts.at(2);
    const int Slon   = myStarts.at(3);

    size_t start[3], count[3];
    start[0] = Sdepth;
    count[0] = Ndepth;

    start[1] = Slat;
    count[1] = Nlat;

    start[2] = Slon;
    count[2] = Nlon;

    for (int Ifield = 0; Ifield < num_fields; ++Ifield) {
        write_field_to_output( time_average.at(Ifield), vars_to_process.at(Ifield) + "_time_average", 
                start, count, filename, &mask );

        /*
        write_time_average_to_post(time_average.at(Ifield), vars_to_process.at(Ifield), 
                "_time_average", start, count, filename, &mask);
        write_time_average_to_post(time_std_dev.at(Ifield), vars_to_process.at(Ifield), 
                "_time_std_dev", start, count, filename, &mask);
        */
    }
}
