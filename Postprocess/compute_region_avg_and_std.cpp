#include <math.h>
#include <mpi.h>
#include <omp.h>
#include <vector>

#include "../constants.hpp"
#include "../functions.hpp"
#include "../postprocess.hpp"

void compute_region_avg_and_std(
        std::vector< std::vector< double > > & field_averages,
        std::vector< std::vector< double > > & field_std_devs,
        const std::vector<const std::vector<double>*> & postprocess_fields,
        const std::vector<double> & region_areas,
        const std::vector<double> & areas,
        const std::vector<bool> & mask,
        const std::vector<double> & latitude,
        const std::vector<double> & longitude,
        const int num_fields,
        const int num_regions,
        const int Ntime,
        const int Ndepth,
        const int Nlat,
        const int Nlon
        ) {

    const int chunk_size = get_omp_chunksize(Nlat, Nlon);

    double int_val, reg_area, dA;

    int Ifield, Iregion, Itime, Idepth, Ilat, Ilon;
    size_t int_index, area_index, index;

    for (Ifield = 0; Ifield < num_fields; ++Ifield) {
        for (Iregion = 0; Iregion < num_regions; ++Iregion) {
            for (Itime = 0; Itime < Ntime; ++Itime) {
                for (Idepth = 0; Idepth < Ndepth; ++Idepth) {

                    int_index = Index(0, Itime, Idepth, Iregion,
                                      1, Ntime, Ndepth, num_regions);
                    reg_area = region_areas.at(int_index);

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

                                if ( mask.at(index) ) { // Skip land areas

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
                    field_averages.at(Ifield).at(int_index) = (reg_area == 0) ? 0. : int_val / reg_area;
                }
            }
        }
    }

    // Now that we have region averages, get region standard deviations
    for (Ifield = 0; Ifield < num_fields; ++Ifield) {
        for (Iregion = 0; Iregion < num_regions; ++Iregion) {
            for (Itime = 0; Itime < Ntime; ++Itime) {
                for (Idepth = 0; Idepth < Ndepth; ++Idepth) {

                    int_index = Index(0, Itime, Idepth, Iregion,
                                      1, Ntime, Ndepth, num_regions);
                    reg_area = region_areas.at(int_index);

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

                                if ( mask.at(index) ) { // Skip land areas

                                    area_index = Index(0, 0, Ilat, Ilon,
                                                       1, 1, Nlat, Nlon);

                                    dA = areas.at(area_index);

                                    if ( RegionTest::all_regions.at(Iregion)(
                                                latitude.at(Ilat), longitude.at(Ilon)) ) 
                                    {
                                        int_val +=
                                            pow(   field_averages.at(Ifield).at(int_index)
                                                 - postprocess_fields.at(Ifield)->at(index),
                                                2.) * dA;
                                    }
                                }
                            }
                        }
                    }
                    field_std_devs.at(Ifield).at(int_index) = (reg_area == 0) ? 0. : sqrt( int_val / reg_area );
                }
            }
        }
    }

}
