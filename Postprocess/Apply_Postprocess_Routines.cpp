#include <math.h>
#include <mpi.h>
#include <omp.h>
#include <vector>

#include "../constants.hpp"
#include "../functions.hpp"
#include "../postprocess.hpp"
#include "../netcdf_io.hpp"

void Apply_Postprocess_Routines(
        const std::vector<double> & coarse_u_r, 
        const std::vector<double> & coarse_u_lon, 
        const std::vector<double> & coarse_u_lat, 
        const std::vector<double> & coarse_vort_r, 
        const std::vector<double> & coarse_vort_lon, 
        const std::vector<double> & coarse_vort_lat,
        const std::vector<double> & div_J, 
        const std::vector<double> & energy_transfer, 
        const std::vector<double> & lambda_m, 
        const std::vector<double> & PEtoKE,
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

    const int Stime  = myStarts.at(0);
    const int Sdepth = myStarts.at(1);
    const int Slat   = myStarts.at(2);
    const int Slon   = myStarts.at(3);

    int wRank=-1, wSize=-1;
    MPI_Comm_rank( MPI_COMM_WORLD, &wRank );
    MPI_Comm_size( MPI_COMM_WORLD, &wSize );

    //
    //// Start by initializing the postprocess file
    //
    
    // Filename
    char filename[50];
    snprintf(filename, 50, "postprocess_%.6gkm.nc", filter_scale/1e3);

    // Variables to integrate
    std::vector<std::string> vars_to_process;
    vars_to_process.push_back("KE");
    vars_to_process.push_back("Pi");
    vars_to_process.push_back("Transport");
    if (constants::COMP_BC_TRANSFERS) {
        vars_to_process.push_back("lambda_m");
        vars_to_process.push_back("PEtoKE");
    }

    initialize_postprocess_file(
            time, depth, latitude, longitude, 
            RegionTest::region_names,
            vars_to_process,
            filename,
            filter_scale
            );


    //
    //// Domain integrals
    //
    double initial_value = 0;
    const int num_fields  = vars_to_process.size();
    const int num_regions = RegionTest::all_regions.size();

    std::vector<std::vector<std::vector<double>>> field_averages(num_fields),
        field_std_devs(num_fields);

    for (int Ifield = 0; Ifield < num_fields; ++Ifield) {
        field_averages.at(Ifield).resize(num_regions,
            std::vector<double>(Ntime * Ndepth, initial_value));
        field_std_devs.at(Ifield).resize(num_regions,
            std::vector<double>(Ntime * Ndepth, initial_value));
    }

    std::vector<double> field_values(num_fields, 0.);

    // Compute the area of each region
    double dA, tmp, local_area;
    int Ifield, Ilat, Ilon, Itime, Idepth, Iregion, 
        area_index, int_index, ave_index;
    size_t index;
    std::vector<double> region_areas(num_regions);
    for (Iregion = 0; Iregion < num_regions; ++Iregion) {

        local_area = 0.;

        #pragma omp parallel default(none)\
        private(Ilat, Ilon, Itime, Idepth, index, dA, area_index )\
        shared(latitude, longitude, areas, mask, Iregion) \
        reduction(+ : local_area)
        { 

            Itime  = 0;
            Idepth = 0;

            #pragma omp for collapse(2) schedule(dynamic)
            for (Ilat = 0; Ilat < Nlat; ++Ilat) {
                for (Ilon = 0; Ilon < Nlon; ++Ilon) {

                    index = Index(Itime, Idepth, Ilat, Ilon,
                                  Ntime, Ndepth, Nlat, Nlon);

                    if (mask.at(index) == 1) { // Skip land areas

                        area_index = Index(0,     0,      Ilat, Ilon,
                                           Ntime, Ndepth, Nlat, Nlon);

                        dA = areas.at(area_index);

                        if ( RegionTest::all_regions.at(Iregion)(
                                    latitude.at(Ilat), longitude.at(Ilon)) ) 
                        { local_area += dA; }
                    }
                }
            }
        }

        region_areas.at(Iregion) = local_area;
    }


    // Compute region averages
    #pragma omp parallel default(none)\
    private(Ifield, Ilat, Ilon, Itime, Idepth, Iregion, \
            index, tmp, dA, area_index, int_index )\
    shared(latitude, longitude, energy_transfer, \
            coarse_u_r, coarse_u_lon, coarse_u_lat, \
            areas, mask, field_averages, field_values, region_areas, \
            lambda_m, PEtoKE, div_J, vars_to_process)
    { 

        #pragma omp for collapse(1) schedule(dynamic)
        for (index = 0; index < energy_transfer.size(); index++) {
            if (mask.at(index) == 1) { // Skip land areas

                Index1to4(index, Itime, Idepth, Ilat, Ilon,
                                 Ntime, Ndepth, Nlat, Nlon);

                int_index = Index(Itime, Idepth, 0, 0,
                                  Ntime, Ndepth, 1, 1);

                area_index = Index(0,     0,      Ilat, Ilon,
                                   Ntime, Ndepth, Nlat, Nlon);

                for (Ifield = 0; Ifield < num_fields; ++Ifield) {

                    tmp = 0.;

                    if (vars_to_process.at(Ifield) == "KE") {
                      tmp = 0.5 * constants::rho0 * (
                             pow(coarse_u_r.at(  index), 2) 
                           + pow(coarse_u_lon.at(index), 2)
                           + pow(coarse_u_lat.at(index), 2) 
                           );
                    } else if (vars_to_process.at(Ifield) == "Pi") {
                        tmp = energy_transfer.at(index);
                    } else if (vars_to_process.at(Ifield) == "Transport") {
                        tmp = div_J.at(index);
                    } else if (vars_to_process.at(Ifield) == "lambda_m") {
                        tmp = lambda_m.at(index);
                    } else if (vars_to_process.at(Ifield) == "PEtoKE") {
                        tmp = PEtoKE.at(index);
                    }

                    field_values.at(Ifield) = tmp;

                }

                dA = areas.at(area_index);

                for (Iregion = 0; Iregion < num_regions; ++Iregion) {
                    if ( RegionTest::all_regions.at(Iregion)(
                                latitude.at(Ilat), longitude.at(Ilon)) ) 
                    {
                        for (Ifield = 0; Ifield < num_fields; ++Ifield) {
                            field_averages.at(Ifield).at(Iregion).at(int_index) +=
                                field_values.at(Ifield) * dA / region_areas.at(Iregion);
                        }
                    }
                }
            }
        }
    }

    // Now that we have region averages, get region standard deviations
    #pragma omp parallel default(none)\
    private(Ifield, Ilat, Ilon, Itime, Idepth, Iregion, \
            index, tmp, dA, area_index, int_index )\
    shared(latitude, longitude, energy_transfer, \
            coarse_u_r, coarse_u_lon, coarse_u_lat, \
            areas, mask, field_averages, field_std_devs, field_values, \
            lambda_m, PEtoKE, div_J, region_areas, vars_to_process) 
    { 
        #pragma omp for collapse(1) schedule(dynamic)
        for (index = 0; index < energy_transfer.size(); index++) {
            if (mask.at(index) == 1) { // Skip land areas

                Index1to4(index, Itime, Idepth, Ilat, Ilon,
                                 Ntime, Ndepth, Nlat, Nlon);

                int_index = Index(Itime, Idepth, 0, 0,
                                  Ntime, Ndepth, 1, 1);

                area_index = Index(0,     0,      Ilat, Ilon,
                                   Ntime, Ndepth, Nlat, Nlon);

                for (Ifield = 0; Ifield < num_fields; ++Ifield) {

                    tmp = 0.;

                    if (vars_to_process.at(Ifield) == "KE") {
                      tmp = 0.5 * constants::rho0 * (
                             pow(coarse_u_r.at(  index), 2) 
                           + pow(coarse_u_lon.at(index), 2)
                           + pow(coarse_u_lat.at(index), 2) 
                           );
                    } else if (vars_to_process.at(Ifield) == "Pi") {
                        tmp = energy_transfer.at(index);
                    } else if (vars_to_process.at(Ifield) == "Transport") {
                        tmp = div_J.at(index);
                    } else if (vars_to_process.at(Ifield) == "lambda_m") {
                        tmp = lambda_m.at(index);
                    } else if (vars_to_process.at(Ifield) == "PEtoKE") {
                        tmp = PEtoKE.at(index);
                    }

                    field_values.at(Ifield) = tmp;

                }

                dA = areas.at(area_index);

                for (Iregion = 0; Iregion < num_regions; ++Iregion) {
                    if ( RegionTest::all_regions.at(Iregion)(
                                latitude.at(Ilat), longitude.at(Ilon)) ) 
                    {
                        for (Ifield = 0; Ifield < num_fields; ++Ifield) {
                            field_std_devs.at(Ifield).at(Iregion).at(int_index) +=
                                pow(   field_averages.at(Ifield).at(Iregion).at(int_index)
                                     - field_values.at(Ifield),
                                    2.) * dA / region_areas.at(Iregion);
                        }
                    }
                }
            }
        }
        #pragma omp for collapse(4) schedule(static)
        for (Ifield = 0; Ifield < num_fields; ++Ifield) {
            for (Iregion = 0; Iregion < num_regions; ++Iregion) {
                for (Itime = 0; Itime < Ntime; ++Itime){
                    for (Idepth = 0; Idepth < Ndepth; ++Idepth){
                        int_index = Index(Itime, Idepth, 0, 0,
                                          Ntime, Ndepth, 1, 1);
                        field_std_devs.at(Ifield).at(Iregion).at(int_index) = 
                            sqrt( field_std_devs.at(Ifield).at(Iregion).at(int_index) );
                    }
                }
            }
        }
    }


    // Dimension order: time - depth - region
    size_t start[3], count[3];
    start[0] = Stime;
    count[0] = Ntime;

    start[1] = Sdepth;
    count[1] = Ndepth;

    count[2] = 1;

    for (int Ifield = 0; Ifield < num_fields; ++Ifield) {
        write_integral_to_post(field_averages.at(Ifield), vars_to_process.at(Ifield), 
                "_avg", start, count, filename);
        write_integral_to_post(field_std_devs.at(Ifield), vars_to_process.at(Ifield), 
                "_std", start, count, filename);
    }

    // Write the region areas (needed for reference, etc)
    size_t start_r[1], count_r[1];
    start_r[0] = 0;
    count_r[0] = (size_t) num_regions;
    write_field_to_output(region_areas, "region_areas", start_r, count_r, filename, NULL);

    // Write region names
    //   this has to be done separately for reasons
    write_regions_to_post(filename);


    //
    //// Time averages
    //

    // Time-average the domain integrals
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
            index, tmp, ave_index )\
    shared(latitude, longitude, energy_transfer, \
            coarse_u_r, coarse_u_lon, coarse_u_lat, \
            areas, mask, time_average_loc,\
            lambda_m, PEtoKE, div_J, vars_to_process)
    { 
        #pragma omp for collapse(4) schedule(dynamic)
        for (Itime = 0; Itime < Ntime; ++Itime){
            for (Idepth = 0; Idepth < Ndepth; ++Idepth){
                for (Ilat = 0; Ilat < Nlat; ++Ilat){
                    for (Ilon = 0; Ilon < Nlon; ++Ilon){

                        // get some indices
                        index = Index(Itime, Idepth, Ilat, Ilon,
                                      Ntime, Ndepth, Nlat, Nlon);
                        ave_index = Index(0, Idepth, Ilat, Ilon,
                                          1, Ndepth, Nlat, Nlon);

                        if (mask.at(index) == 1) { // Skip land areas
                            for (Ifield = 0; Ifield < num_fields; ++Ifield) {

                                // get local field value
                                tmp = 0.;

                                if (vars_to_process.at(Ifield) == "KE") {
                                    tmp = 0.5 * constants::rho0 * (
                                              pow(coarse_u_r.at(  index), 2) 
                                            + pow(coarse_u_lon.at(index), 2)
                                            + pow(coarse_u_lat.at(index), 2) 
                                            );
                                } else if (vars_to_process.at(Ifield) == "Pi") {
                                    tmp = energy_transfer.at(index);
                                } else if (vars_to_process.at(Ifield) == "Transport") {
                                    tmp = div_J.at(index);
                                } else if (vars_to_process.at(Ifield) == "lambda_m") {
                                    tmp = lambda_m.at(index);
                                } else if (vars_to_process.at(Ifield) == "PEtoKE") {
                                    tmp = PEtoKE.at(index);
                                }

                                // compute the time average for the part on this processor
                                time_average_loc.at(Ifield).at(ave_index) += tmp;
                            }
                        } else {
                            for (Ifield = 0; Ifield < num_fields; ++Ifield) {
                                time_average_loc.at(Ifield).at(ave_index) = 
                                    Ntime * constants::fill_value;
                            }
                        }
                    }
                }
            }
        }
    }

    // Now communicate with other processors to get the full time average
    for (Ifield = 0; Ifield < num_fields; ++Ifield) {
        MPI_Allreduce(&(time_average_loc.at(Ifield)[0]),
                      &(time_average.at(    Ifield)[0]),
                      Ndepth * Nlat * Nlon, MPI_DOUBLE, MPI_SUM, comm);
    }
    int full_Ntime;
    MPI_Allreduce(&Ntime, &full_Ntime, 1, MPI_INT, MPI_SUM, comm);


    // Compute the standard deviation
    #pragma omp parallel default(none)\
    private(Ifield, Ilat, Ilon, Itime, Idepth, \
            index, tmp, ave_index )\
    shared(latitude, longitude, full_Ntime, wSize, energy_transfer, \
            coarse_u_r, coarse_u_lon, coarse_u_lat, \
            areas, mask, time_average_loc, time_std_dev_loc, time_average, \
            lambda_m, PEtoKE, div_J, vars_to_process)
    { 
        #pragma omp for collapse(4) schedule(dynamic)
        for (Itime = 0; Itime < Ntime; ++Itime){
            for (Idepth = 0; Idepth < Ndepth; ++Idepth){
                for (Ilat = 0; Ilat < Nlat; ++Ilat){
                    for (Ilon = 0; Ilon < Nlon; ++Ilon){

                        // get some indices
                        index = Index(Itime, Idepth, Ilat, Ilon,
                                      Ntime, Ndepth, Nlat, Nlon);
                        ave_index = Index(0, Idepth, Ilat, Ilon,
                                          1, Ndepth, Nlat, Nlon);

                        if (mask.at(index) == 1) { // Skip land areas
                            for (Ifield = 0; Ifield < num_fields; ++Ifield) {

                                // get local field value
                                tmp = 0.;

                                if (vars_to_process.at(Ifield) == "KE") {
                                    tmp = 0.5 * constants::rho0 * (
                                              pow(coarse_u_r.at(  index), 2) 
                                            + pow(coarse_u_lon.at(index), 2)
                                            + pow(coarse_u_lat.at(index), 2) 
                                            );
                                } else if (vars_to_process.at(Ifield) == "Pi") {
                                    tmp = energy_transfer.at(index);
                                } else if (vars_to_process.at(Ifield) == "Transport") {
                                    tmp = div_J.at(index);
                                } else if (vars_to_process.at(Ifield) == "lambda_m") {
                                    tmp = lambda_m.at(index);
                                } else if (vars_to_process.at(Ifield) == "PEtoKE") {
                                    tmp = PEtoKE.at(index);
                                }

                                // compute the time average for the part on this processor
                                time_std_dev_loc.at(Ifield).at(ave_index) += 
                                    pow(tmp - time_average.at(Ifield).at(ave_index), 
                                        2.) / full_Ntime;
                            }
                        } else {
                            for (Ifield = 0; Ifield < num_fields; ++Ifield) {
                                time_std_dev_loc.at(Ifield).at(ave_index) = 
                                    pow(constants::fill_value, 2.) / wSize;
                            }
                        }
                    }
                }
            }
        }
    }

    // Now communicate with other processors to get the full time average
    for (Ifield = 0; Ifield < num_fields; ++Ifield) {
        MPI_Allreduce(&(time_std_dev_loc.at(Ifield)[0]),
                      &(time_std_dev.at(    Ifield)[0]),
                      Ndepth * Nlat * Nlon, MPI_DOUBLE, MPI_SUM, comm);
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
        // Scale by full_Ntime
        #pragma omp parallel default(none)\
        private(index)\
        shared(Ifield, time_average, time_std_dev, full_Ntime, mask)
        { 
            #pragma omp for collapse(1) schedule(dynamic)
            for (index = 0; index < time_average.at(Ifield).size(); ++index) {
                time_average.at(Ifield).at(index) *= 1./full_Ntime;
                time_std_dev.at(Ifield).at(index) *= 
                    sqrt(time_std_dev.at(Ifield).at(index));
            }
        }
        write_time_average_to_post(time_average.at(Ifield), vars_to_process.at(Ifield), 
                "_time_average", start, count, filename);
        write_time_average_to_post(time_std_dev.at(Ifield), vars_to_process.at(Ifield), 
                "_time_std_dev", start, count, filename);
    }
}
