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

    const int Ntime   = myCounts.at(0);
    const int Ndepth  = myCounts.at(1);
    const int Nlat    = myCounts.at(2);
    const int Nlon    = myCounts.at(3);

    const int Stime   = myStarts.at(0);
    const int Sdepth  = myStarts.at(1);

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
    std::vector<std::string> integrated_vars;
    integrated_vars.push_back("Pi");
    integrated_vars.push_back("KE");
    if (constants::COMP_BC_TRANSFERS) {
        integrated_vars.push_back("lambda_m");
        integrated_vars.push_back("PEtoKE");
    }


    initialize_postprocess_file(
            time, depth, latitude, longitude, 
            RegionTest::region_names,
            integrated_vars,
            filename,
            filter_scale
            );


    //
    //// Domain integrals
    //
    
    double initial_value = 0;
    const int num_regions = RegionTest::all_regions.size();

    std::vector<std::vector<double>> Pi_integrals, KE_integrals, 
        Lambda_integrals, PEtoKE_integrals;
    Pi_integrals.resize(num_regions, std::vector<double>(Ntime * Ndepth, initial_value));
    KE_integrals.resize(num_regions, std::vector<double>(Ntime * Ndepth, initial_value));
    if (constants::COMP_BC_TRANSFERS) {
        Lambda_integrals.resize(num_regions, 
                std::vector<double>(Ntime * Ndepth, initial_value));
        PEtoKE_integrals.resize(num_regions, 
                std::vector<double>(Ntime * Ndepth, initial_value));
    }

    std::vector<double> region_areas(num_regions);

    double curr_Pi, curr_KE, curr_Lambda, curr_PEtoKE, dA;
    int Ilat, Ilon, Itime, Idepth, Iregion, area_index, int_index;
    size_t index;
    #pragma omp parallel default(none)\
    private(Ilat, Ilon, Itime, Idepth, Iregion, \
            index, curr_Pi, curr_KE, curr_Lambda, curr_PEtoKE, \
            dA, area_index, int_index )\
    shared(latitude, longitude, energy_transfer, \
            coarse_u_r, coarse_u_lon, coarse_u_lat, \
            Pi_integrals, KE_integrals, areas, mask,\
            Lambda_integrals, PEtoKE_integrals,\
            region_areas, lambda_m, PEtoKE) \
    firstprivate(wRank)
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

                curr_Pi = energy_transfer.at(index);
                curr_KE = (  pow(coarse_u_r.at(  index), 2) 
                           + pow(coarse_u_lon.at(index), 2)
                           + pow(coarse_u_lat.at(index), 2) );
                if (constants::COMP_BC_TRANSFERS) {
                    curr_Lambda = lambda_m.at(index);
                    curr_PEtoKE = PEtoKE.at(index);
                }
                dA = areas.at(area_index);

                for (Iregion = 0; Iregion < num_regions; ++Iregion) {
                    if ( RegionTest::all_regions.at(Iregion)(
                                latitude.at(Ilat), longitude.at(Ilon)) ) 
                    {
                        Pi_integrals.at(Iregion).at(int_index) += curr_Pi * dA;
                        KE_integrals.at(Iregion).at(int_index) += curr_KE * dA;

                        if (constants::COMP_BC_TRANSFERS) {
                            Lambda_integrals.at(Iregion).at(int_index) += curr_Lambda * dA;
                            PEtoKE_integrals.at(Iregion).at(int_index) += curr_PEtoKE * dA;
                        }

                        if ( (Itime == 0) and (Idepth == 0) and (wRank == 0)) {
                            region_areas.at(Iregion) += dA;
                        }
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

    write_integral_to_post(Pi_integrals, "Pi", start, count, filename);
    write_integral_to_post(KE_integrals, "KE", start, count, filename);
    if (constants::COMP_BC_TRANSFERS) {
        write_integral_to_post(Lambda_integrals, "lambda_m", start, count, filename);
        write_integral_to_post(PEtoKE_integrals, "PEtoKE",   start, count, filename);
    }

    // Write the region areas (needed for normalization)
    size_t start_r[1], count_r[1];
    start_r[0] = 0;
    count_r[0] = (size_t) num_regions;
    write_field_to_output(region_areas, "region_areas", start_r, count_r, filename, NULL);

    // Write region names
    //   this has to be done separately for reasons
    write_regions_to_post(filename);
}
