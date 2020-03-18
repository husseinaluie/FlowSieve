#include <math.h>
#include <mpi.h>
#include <omp.h>
#include <vector>

#include "../constants.hpp"
#include "../functions.hpp"
#include "../postprocess.hpp"
#include "../netcdf_io.hpp"

void Apply_Postprocess_Routines(
        const std::vector<const std::vector<double>*> & postprocess_fields,
        const std::vector<std::string> & vars_to_process,
        const std::vector<double> & time,
        const std::vector<double> & depth,
        const std::vector<double> & latitude,
        const std::vector<double> & longitude,
        const std::vector<double> & mask,
        const std::vector<double> & areas,
        const std::vector<int>    & myCounts,
        const std::vector<int>    & myStarts,
        const double filter_scale,
        const MPI_Comm comm
        ) {

    const int Ntime  = myCounts.at(0);
    const int Ndepth = myCounts.at(1);
    const int Nlat   = myCounts.at(2);
    const int Nlon   = myCounts.at(3);

    const int chunk_size = get_omp_chunksize(Nlat, Nlon);

    // Get full number of time points
    int full_Ntime;
    MPI_Allreduce(&Ntime, &full_Ntime, 1, MPI_INT, MPI_SUM, comm);

    const int Stime  = myStarts.at(0);
    const int Sdepth = myStarts.at(1);
    const int Slat   = myStarts.at(2);
    const int Slon   = myStarts.at(3);

    int wRank=-1, wSize=-1;
    MPI_Comm_rank( MPI_COMM_WORLD, &wRank );
    MPI_Comm_size( MPI_COMM_WORLD, &wSize );

    double dA, local_area;
    int Ifield, Ilat, Ilon, Itime, Idepth, Iregion, 
        area_index, int_index, reg_index, space_index;

    size_t index;

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
            mask_count_loc.at(area_index) += (int) mask.at(index);
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
    
    // Filename + file
    snprintf(filename, 50, "postprocess_%.6gkm.nc", filter_scale/1e3);
    initialize_postprocess_file(
            time, depth, latitude, longitude, 
            RegionTest::region_names,
            vars_to_process,
            filename,
            filter_scale
            );

    //
    //// Compute the area of each region
    //
    const int num_fields  = vars_to_process.size();
    const int num_regions = RegionTest::all_regions.size();

    if (wRank == 0) { fprintf(stdout, "  .. computing the area of each region\n"); }
    fflush(stdout);

    std::vector<double> region_areas(num_regions * Ntime * Ndepth, 0.);
    for (Iregion = 0; Iregion < num_regions; ++Iregion) {
        for (Itime = 0; Itime < Ntime; ++Itime) {
            for (Idepth = 0; Idepth < Ndepth; ++Idepth) {

                reg_index = Index(0, Itime, Idepth, Iregion,
                                  1, Ntime, Ndepth, num_regions);

                local_area = 0.;

                #pragma omp parallel default(none)\
                private(Ilat, Ilon, index, dA, area_index )\
                shared(latitude, longitude, areas, mask, Iregion, Itime, Idepth) \
                reduction(+ : local_area)
                { 
                    #pragma omp for collapse(2) schedule(guided, chunk_size)
                    for (Ilat = 0; Ilat < Nlat; ++Ilat) {
                        for (Ilon = 0; Ilon < Nlon; ++Ilon) {

                            index = Index(Itime, Idepth, Ilat, Ilon,
                                          Ntime, Ndepth, Nlat, Nlon);

                            if (mask.at(index) == 1) { // Skip land areas

                                area_index = Index(0, 0, Ilat, Ilon,
                                                   1, 1, Nlat, Nlon);

                                dA = areas.at(area_index);

                                if ( RegionTest::all_regions.at(Iregion)(
                                            latitude.at(Ilat), longitude.at(Ilon)) ) 
                                { local_area += dA; }
                            }
                        }
                    }
                }

                region_areas.at(reg_index) = local_area;
            }
        }
    }


    //
    //// Domain integrals
    //
    std::vector< std::vector< std::vector< double > > > 
        field_averages(num_fields), field_std_devs(num_fields);

    for (int Ifield = 0; Ifield < num_fields; ++Ifield) {
        field_averages.at(Ifield).resize(num_regions,
            std::vector<double>(Ntime * Ndepth, 0.));
        field_std_devs.at(Ifield).resize(num_regions,
            std::vector<double>(Ntime * Ndepth, 0.));
    }

    std::vector<double> field_values(num_fields, 0.);

    // Compute region averages
    if (wRank == 0) { fprintf(stdout, "  .. computing region integrals\n"); }
    fflush(stdout);

    double int_val, reg_area;
    for (Ifield = 0; Ifield < num_fields; ++Ifield) {
        for (Iregion = 0; Iregion < num_regions; ++Iregion) {
            for (Itime = 0; Itime < Ntime; ++Itime) {
                for (Idepth = 0; Idepth < Ndepth; ++Idepth) {
                    int_index = Index(Itime, Idepth, 0, 0,
                                      Ntime, Ndepth, 1, 1);

                    reg_index = Index(0, Itime, Idepth, Iregion,
                                      1, Ntime, Ndepth, num_regions);
                    reg_area = region_areas.at(reg_index);

                    int_val = 0.;

                    #pragma omp parallel default(none)\
                    private(Ilat, Ilon, index, dA, area_index )\
                    shared(Ifield, Iregion, Itime, Idepth, int_index, mask, areas,\
                            latitude, longitude, postprocess_fields) \
                    reduction(+ : int_val)
                    { 
                        #pragma omp for collapse(2) schedule(guided, chunk_size)
                        for (Ilat = 0; Ilat < Nlat; ++Ilat) {
                            for (Ilon = 0; Ilon < Nlon; ++Ilon) {

                                index = Index(Itime, Idepth, Ilat, Ilon,
                                              Ntime, Ndepth, Nlat, Nlon);

                                if (mask.at(index) == 1) { // Skip land areas

                                    area_index = Index(0, 0, Ilat, Ilon,
                                                       1, 1, Nlat, Nlon);

                                    dA = areas.at(area_index);

                                    if ( RegionTest::all_regions.at(Iregion)(
                                                latitude.at(Ilat), longitude.at(Ilon)) ) 
                                    {
                                        int_val +=  postprocess_fields.at(Ifield)->at(index) * dA;
                                    }
                                }
                            }
                        }
                    }
                    field_averages.at(Ifield).at(Iregion).at(int_index) = int_val / reg_area;
                }
            }
        }
    }

    // Now that we have region averages, get region standard deviations
    if (wRank == 0) { fprintf(stdout, "  .. computing region standard deviations\n"); }
    fflush(stdout);

    // Now that we have region averages, get region standard deviations
    for (Ifield = 0; Ifield < num_fields; ++Ifield) {
        for (Iregion = 0; Iregion < num_regions; ++Iregion) {
            for (Itime = 0; Itime < Ntime; ++Itime) {
                for (Idepth = 0; Idepth < Ndepth; ++Idepth) {
                    int_index = Index(Itime, Idepth, 0, 0,
                                      Ntime, Ndepth, 1, 1);

                    reg_index = Index(0, Itime, Idepth, Iregion,
                                      1, Ntime, Ndepth, num_regions);
                    reg_area = region_areas.at(reg_index);

                    int_val = 0.;

                    #pragma omp parallel default(none)\
                    private(Ilat, Ilon, index, dA, area_index )\
                    shared(Ifield, Iregion, Itime, Idepth, int_index, mask, areas,\
                            latitude, longitude, postprocess_fields, field_averages) \
                    reduction(+ : int_val)
                    { 
                        #pragma omp for collapse(2) schedule(guided, chunk_size)
                        for (Ilat = 0; Ilat < Nlat; ++Ilat) {
                            for (Ilon = 0; Ilon < Nlon; ++Ilon) {

                                index = Index(Itime, Idepth, Ilat, Ilon,
                                              Ntime, Ndepth, Nlat, Nlon);

                                if (mask.at(index) == 1) { // Skip land areas

                                    area_index = Index(0, 0, Ilat, Ilon,
                                                       1, 1, Nlat, Nlon);

                                    dA = areas.at(area_index);

                                    if ( RegionTest::all_regions.at(Iregion)(
                                                latitude.at(Ilat), longitude.at(Ilon)) ) 
                                    {
                                        int_val +=
                                            pow(   field_averages.at(Ifield).at(Iregion).at(int_index)
                                                 - postprocess_fields.at(Ifield)->at(index),
                                                2.) * dA;
                                    }
                                }
                            }
                        }
                    }
                    field_std_devs.at(Ifield).at(Iregion).at(int_index) = sqrt( int_val / reg_area );
                }
            }
        }
    }

    if (wRank == 0) { fprintf(stdout, "  .. writing region averages and deviations\n"); }
    fflush(stdout);

    // Dimension order: time - depth - region
    size_t start[3], count[3];
    start[0] = Stime;
    count[0] = Ntime;

    start[1] = Sdepth;
    count[1] = Ndepth;

    // start[2] will be set in the write_integral_to_post function
    count[2] = 1;

    for (Ifield = 0; Ifield < num_fields; ++Ifield) {
        write_integral_to_post(field_averages.at(Ifield), vars_to_process.at(Ifield), 
                "_avg", start, count, filename);
        write_integral_to_post(field_std_devs.at(Ifield), vars_to_process.at(Ifield), 
                "_std", start, count, filename);

        if (wRank == 0) { fprintf(stdout, "  .. .. wrote field index %d\n", Ifield); }
        fflush(stdout);
    }
    if (wRank == 0) { fprintf(stdout, "  .. done writing fields %d\n", Ifield); }
    fflush(stdout);

    // Write the region areas (needed for reference, etc)
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
    //// Time averages
    //

    if (wRank == 0) { fprintf(stdout, "  .. computing time-averages of fields\n"); }
    fflush(stdout);

    // Time-average the fields
    std::vector<std::vector<double>> 
        time_average_loc(num_fields), time_average(num_fields), 
        time_std_dev_loc(num_fields), time_std_dev(num_fields);
    for (Ifield = 0; Ifield < num_fields; ++Ifield) {
        time_average_loc.at(Ifield).resize( Ndepth * Nlat * Nlon, 0. );
        time_average.at(    Ifield).resize( Ndepth * Nlat * Nlon, 0. );

        time_std_dev_loc.at(Ifield).resize( Ndepth * Nlat * Nlon, 0. );
        time_std_dev.at(    Ifield).resize( Ndepth * Nlat * Nlon, 0. );
    }

    #pragma omp parallel default(none)\
    private(Ifield, Ilat, Ilon, Itime, Idepth, \
            index, space_index )\
    shared(latitude, longitude, postprocess_fields, \
            areas, mask, always_masked, mask_count, \
            time_average_loc, full_Ntime)
    { 
        #pragma omp for collapse(3) schedule(guided, chunk_size)
        for (Ilat = 0; Ilat < Nlat; ++Ilat){
            for (Ilon = 0; Ilon < Nlon; ++Ilon){
                for (Idepth = 0; Idepth < Ndepth; ++Idepth){
                    space_index = Index(0, 0, Ilat, Ilon,
                                        1, 1, Nlat, Nlon);
                    if (not(always_masked.at(space_index))) { // Skip land areas
                        for (Itime = 0; Itime < Ntime; ++Itime){

                            // get some indices
                            index = Index(Itime, Idepth, Ilat, Ilon,
                                          Ntime, Ndepth, Nlat, Nlon);

                            if (mask.at(index) == 1) {
                                for (Ifield = 0; Ifield < num_fields; ++Ifield) {

                                    // compute the time average for 
                                    // the part on this processor
                                    time_average_loc.at(Ifield).at(space_index) 
                                        += postprocess_fields.at(Ifield)->at(index) 
                                                / mask_count.at(space_index);
                                }
                            }
                        }
                    } else {
                        for (Ifield = 0; Ifield < num_fields; ++Ifield) {
                            time_average_loc.at(Ifield).at(space_index) 
                                = constants::fill_value;
                        }
                    }
                }
            }
        }
    }

    // Now communicate with other processors to get the full time average
    if (wRank == 0) { fprintf(stdout, "  .. .. reducing across processors\n"); }
    fflush(stdout);

    for (Ifield = 0; Ifield < num_fields; ++Ifield) {
        MPI_Allreduce(&(time_average_loc.at(Ifield)[0]),
                      &(time_average.at(    Ifield)[0]),
                      Ndepth * Nlat * Nlon, MPI_DOUBLE, MPI_SUM, comm);
    }

    // Compute the standard deviation
    if (wRank == 0) { fprintf(stdout, "  .. computing time standard deviations\n"); }
    fflush(stdout);

    #pragma omp parallel default(none)\
    private(Ifield, Ilat, Ilon, Itime, Idepth, \
            index, space_index )\
    shared(latitude, longitude, full_Ntime, wSize, postprocess_fields, \
            areas, mask, always_masked, mask_count, \
            time_average_loc, time_std_dev_loc, time_average)
    { 
        #pragma omp for collapse(3) schedule(guided, chunk_size)
        for (Ilat = 0; Ilat < Nlat; ++Ilat){
            for (Ilon = 0; Ilon < Nlon; ++Ilon){
                for (Idepth = 0; Idepth < Ndepth; ++Idepth){
                    space_index = Index(0, Idepth, Ilat, Ilon,
                                        1, Ndepth, Nlat, Nlon);
                    if (not(always_masked.at(space_index))) { // Skip land areas
                        for (Itime = 0; Itime < Ntime; ++Itime){

                            // get some indices
                            index = Index(Itime, Idepth, Ilat, Ilon,
                                          Ntime, Ndepth, Nlat, Nlon);

                            if (mask.at(index) == 1) { // Skip land areas
                                for (Ifield = 0; Ifield < num_fields; ++Ifield) {

                                    // compute the time std. dev. for 
                                    // the part on this processor
                                    time_std_dev_loc.at(Ifield).at(space_index) += 
                                        pow(   postprocess_fields.at(Ifield)->at(index) 
                                             - time_average.at(Ifield).at(space_index), 
                                            2.) / mask_count.at(space_index);
                                }
                            }
                        }
                    } else {
                        for (Ifield = 0; Ifield < num_fields; ++Ifield) {
                            // We'll handle these later to avoid over-flow issues
                            time_std_dev_loc.at(Ifield).at(space_index) = 0.; 
                        }
                    }
                }
            }
        }
    }

    if (wRank == 0) { fprintf(stdout, "  .. writing time averages and deviations\n"); }
    fflush(stdout);

    // Now communicate with other processors to get the full time average
    for (Ifield = 0; Ifield < num_fields; ++Ifield) {
        MPI_Allreduce(&(time_std_dev_loc.at(Ifield)[0]),
                      &(time_std_dev.at(    Ifield)[0]),
                      Ndepth * Nlat * Nlon, MPI_DOUBLE, MPI_SUM, comm);
    }

    // Re-mask to fill in land areas / apply sqrt to water
    for (int Ifield = 0; Ifield < num_fields; ++Ifield) {
        #pragma omp parallel default(none)\
        private(index, space_index)\
        shared(Ifield, time_average, time_std_dev, full_Ntime, always_masked)
        { 
            #pragma omp for collapse(3) schedule(guided, chunk_size)
            for (Idepth = 0; Idepth < Ndepth; ++Idepth){
                for (Ilat = 0; Ilat < Nlat; ++Ilat){
                    for (Ilon = 0; Ilon < Nlon; ++Ilon){
                        space_index = Index(0, Idepth, Ilat, Ilon,
                                            1, Ndepth, Nlat, Nlon);
                        if (not(always_masked.at(space_index))) { // Skip land areas
                            time_std_dev.at(Ifield).at(index) = 
                                sqrt(time_std_dev.at(Ifield).at(index));
                        } else {
                            time_std_dev.at(Ifield).at(index) = constants::fill_value;
                            time_average.at(Ifield).at(index) = constants::fill_value;
                        }
                    }
                }
            }
        }
    }


    // Write the time averages
    //   dimension order: depth - lat - lon
    start[0] = Sdepth;
    count[0] = Ndepth;

    start[1] = Slat;
    count[1] = Nlat;

    start[2] = Slon;
    count[2] = Nlon;

    for (int Ifield = 0; Ifield < num_fields; ++Ifield) {
        write_time_average_to_post(time_average.at(Ifield), vars_to_process.at(Ifield), 
                "_time_average", start, count, filename);
        write_time_average_to_post(time_std_dev.at(Ifield), vars_to_process.at(Ifield), 
                "_time_std_dev", start, count, filename);
    }
}
