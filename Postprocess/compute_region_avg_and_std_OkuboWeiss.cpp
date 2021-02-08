#include <math.h>
#include <mpi.h>
#include <omp.h>
#include <vector>

#include "../constants.hpp"
#include "../functions.hpp"
#include "../postprocess.hpp"

void compute_region_avg_and_std_OkuboWeiss(
        std::vector< std::vector< double > > & field_averages,
        std::vector< std::vector< double > > & field_std_devs,
        std::vector< double > & OkuboWeiss_areas,
        const std::vector<const std::vector<double>*> & postprocess_fields,
        const std::vector<double> & OkuboWeiss,
        const std::vector<double> & region_areas,
        const std::vector<double> & areas,
        const std::vector<bool> & mask,
        const std::vector<double> OkuboWeiss_bounds,
        const std::vector<double> & latitude,
        const std::vector<double> & longitude,
        const int num_fields,
        const int num_regions,
        const int Ntime,
        const int Ndepth,
        const int Nlat,
        const int Nlon,
        const int NOkubo
        ) {

    double reg_area, dA;

    int Ifield, Iregion, Itime, Idepth, Ilat, Ilon, IOkubo;
    size_t int_index, reg_index, area_index, index;

    std::vector<double> int_vals, OW_areas;

    #pragma omp parallel default(none)\
    private( Ifield, Iregion, Itime, Idepth, Ilat, Ilon, IOkubo,\
            index, reg_index, area_index, int_index, \
            int_vals, OW_areas, dA, reg_area )\
    shared( mask, areas, latitude, longitude, postprocess_fields, OkuboWeiss,\
            field_averages, OkuboWeiss_areas, region_areas )
    { 
    #pragma omp for collapse(4) schedule(static)
        for (Ifield = 0; Ifield < num_fields; ++Ifield) {
            for (Iregion = 0; Iregion < num_regions; ++Iregion) {
                for (Itime = 0; Itime < Ntime; ++Itime) {
                    for (Idepth = 0; Idepth < Ndepth; ++Idepth) {

                        reg_index = Index(0, Itime, Idepth, Iregion,
                                          1, Ntime, Ndepth, num_regions);
                        reg_area = region_areas.at(reg_index);

                        int_vals.clear(); int_vals.resize(NOkubo, 0.);
                        OW_areas.clear(); OW_areas.resize(NOkubo, 0.);

                        for (Ilat = 0; Ilat < Nlat; ++Ilat) {
                            for (Ilon = 0; Ilon < Nlon; ++Ilon) {

                                index = Index(Itime, Idepth, Ilat, Ilon,
                                              Ntime, Ndepth, Nlat, Nlon);

                                if ( mask.at(index) ) { // Skip land areas

                                    IOkubo = std::lower_bound( OkuboWeiss_bounds.begin(), OkuboWeiss_bounds.end(), OkuboWeiss.at(index) ) 
                                             - OkuboWeiss_bounds.begin();
                                    if (IOkubo >= NOkubo) { IOkubo = NOkubo - 1; }

                                    if ( RegionTest::all_regions.at(Iregion)( latitude.at(Ilat), longitude.at(Ilon)) ) {
                                        area_index = Index(0, 0, Ilat, Ilon,
                                                           1, 1, Nlat, Nlon);
                                        dA = areas.at(area_index);

                                        int_vals.at(IOkubo) += postprocess_fields.at(Ifield)->at(index) * dA;
                                        OW_areas.at(IOkubo) += dA;
                                    }
                                }
                            }
                        }

                        for (IOkubo = 0; IOkubo < NOkubo; ++IOkubo) {
                            int_index = Index(Itime, Idepth, IOkubo, Iregion,
                                              Ntime, Ndepth, NOkubo, num_regions);

                            field_averages.at(Ifield).at(int_index) = (reg_area == 0) ? 0. : int_vals.at(IOkubo) / reg_area;
                            OkuboWeiss_areas.at(int_index) = OW_areas.at(IOkubo);
                        }
                    }
                }
            }
        }
    }
}
