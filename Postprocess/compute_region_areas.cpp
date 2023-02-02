#include <math.h>
#include <mpi.h>
#include <omp.h>
#include <vector>

#include "../constants.hpp"
#include "../functions.hpp"
#include "../postprocess.hpp"

void compute_region_areas(
        std::vector<double> & region_areas,
        const std::vector<double> & areas,
        const std::vector<bool> & mask,
        const std::vector<double> & latitude,
        const std::vector<double> & longitude,
        const int num_regions,
        const int Ntime,
        const int Ndepth,
        const int Nlat,
        const int Nlon
        ) {

    double dA, local_area;
    int Iregion, Itime, Idepth, Ilat, Ilon;
    size_t reg_index, index, area_index;

    for (Iregion = 0; Iregion < num_regions; ++Iregion) {
        for (Itime = 0; Itime < Ntime; ++Itime) {
            for (Idepth = 0; Idepth < Ndepth; ++Idepth) {

                reg_index = Index(0, Itime, Idepth, Iregion,
                                  1, Ntime, Ndepth, num_regions);

                local_area = 0.;

                #pragma omp parallel default(none)\
                private(Ilat, Ilon, index, dA, area_index )\
                shared(latitude, longitude, areas, mask, Iregion, Itime, Idepth) \
                firstprivate( Nlon, Nlat, Ndepth, Ntime, RegionTest::all_regions ) \
                reduction(+ : local_area)
                { 
                    #pragma omp for collapse(2) schedule(guided)
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
                                { local_area += dA; }
                            }
                        }
                    }
                }

                region_areas.at(reg_index) = local_area;
            }
        }
    }
}
