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
        const dataset & source_data,
        const std::vector<const std::vector<double>*> & postprocess_fields,
        const std::vector<double> & OkuboWeiss,
        const std::vector<double> OkuboWeiss_bounds,
        const int NOkubo
        ) {

    const int   Ntime   = source_data.Ntime,
                Ndepth  = source_data.Ndepth,
                Nlat    = source_data.Nlat,
                Nlon    = source_data.Nlon;

    const int   num_regions   = source_data.region_names.size(),
                num_fields    = postprocess_fields.size();

    double reg_area, dA;

    int Ifield, Iregion, Itime, Idepth, Ilat, Ilon, IOkubo;
    size_t int_index, reg_index, area_index, index;

    std::vector<double> int_vals, OW_areas;

    #pragma omp parallel default(none)\
    private( Ifield, Iregion, Itime, Idepth, Ilat, Ilon, IOkubo,\
            index, reg_index, area_index, int_index, \
            int_vals, OW_areas, dA, reg_area )\
    shared( source_data, postprocess_fields, OkuboWeiss, field_averages, OkuboWeiss_areas ) \
    firstprivate( Nlon, Nlat, Ndepth, Ntime, NOkubo, num_regions, num_fields, OkuboWeiss_bounds )
    { 
    #pragma omp for collapse(4) schedule(static)
        for (Ifield = 0; Ifield < num_fields; ++Ifield) {
            for (Iregion = 0; Iregion < num_regions; ++Iregion) {
                for (Itime = 0; Itime < Ntime; ++Itime) {
                    for (Idepth = 0; Idepth < Ndepth; ++Idepth) {

                        reg_index = Index(0, Itime, Idepth, Iregion,
                                          1, Ntime, Ndepth, num_regions);
                        reg_area = source_data.region_areas.at(reg_index);

                        int_vals.clear(); int_vals.resize(NOkubo, 0.);
                        OW_areas.clear(); OW_areas.resize(NOkubo, 0.);

                        for (Ilat = 0; Ilat < Nlat; ++Ilat) {
                            for (Ilon = 0; Ilon < Nlon; ++Ilon) {

                                index = Index(Itime, Idepth, Ilat, Ilon,
                                              Ntime, Ndepth, Nlat, Nlon);

                                if ( source_data.mask.at(index) ) { // Skip land areas

                                    area_index = Index(0, 0, Ilat, Ilon,
                                                       1, 1, Nlat, Nlon);

                                    IOkubo = std::lower_bound( OkuboWeiss_bounds.begin(), OkuboWeiss_bounds.end(), OkuboWeiss.at(index) ) 
                                             - OkuboWeiss_bounds.begin();
                                    if (IOkubo >= NOkubo) { IOkubo = NOkubo - 1; }

                                    //if ( RegionTest::all_regions.at(Iregion)( latitude.at(Ilat), longitude.at(Ilon)) ) {
                                    if ( source_data.regions.at( source_data.region_names.at( Iregion ) ).at( area_index ) ) {
                                        dA = source_data.areas.at(area_index);

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
