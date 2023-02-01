#include <math.h>
#include <mpi.h>
#include <omp.h>
#include <vector>
#include <algorithm>

#include "../constants.hpp"
#include "../functions.hpp"
#include "../postprocess.hpp"


// Pragma-magic to allow reduction over vector operator
//      thanks to: https://stackoverflow.com/questions/43168661/openmp-and-reduction-on-stdvector
#pragma omp declare reduction(vec_double_plus : std::vector<double> : \
                              std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<double>())) \
                    initializer(omp_priv = decltype(omp_orig)(omp_orig.size()))


void compute_region_avg_and_std(
        std::vector< std::vector< double > > & field_averages,
        std::vector< std::vector< double > > & field_std_devs,
        const dataset & source_data,
        const std::vector<const std::vector<double>*> & postprocess_fields,
        const MPI_Comm comm
        ) {

    const int   Ntime   = source_data.Ntime,
                Ndepth  = source_data.Ndepth,
                Nlat    = source_data.Nlat,
                Nlon    = source_data.Nlon;

    const int   num_regions   = source_data.region_names.size(),
                num_fields    = postprocess_fields.size();

    double reg_area, dA, increment;

    int Ifield, Iregion, Itime, Idepth, Ilat, Ilon;
    size_t int_index, area_index, index;

    #if DEBUG >= 1
    int wRank;
    MPI_Comm_rank( comm, &wRank );
    #endif

    #if DEBUG >= 1
    if (wRank == 0) { fprintf(stdout, "  Computing region means\n"); }
    #endif
    for (Ifield = 0; Ifield < num_fields; ++Ifield) {
        #if DEBUG >= 2
        if (wRank == 0) { fprintf(stdout, "    processing field %d of %d means\n", Ifield + 1, num_fields); }
        #endif

        std::vector<double> field_integrals(Ntime*Ndepth*num_regions,0.);

        #pragma omp parallel default(none)\
        private(Ilat, Ilon, index, dA, area_index, increment, int_index, \
                Idepth, Itime, Iregion )\
        shared( source_data, Ifield, postprocess_fields) \
        reduction(vec_double_plus : field_integrals)
        { 
            #pragma omp for collapse(5) schedule(static)
            for (Iregion = 0; Iregion < num_regions; ++Iregion) {
                for (Itime = 0; Itime < Ntime; ++Itime) {
                    for (Idepth = 0; Idepth < Ndepth; ++Idepth) {
                        for (Ilat = 0; Ilat < Nlat; ++Ilat) {
                            for (Ilon = 0; Ilon < Nlon; ++Ilon) {
                                int_index = Index(0, Itime, Idepth, Iregion, 1, Ntime, Ndepth, num_regions);
                                index = Index(Itime, Idepth, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon);

                                if ( constants::FILTER_OVER_LAND or source_data.mask.at(index) ) { // Skip land areas

                                    area_index = Index(0, 0, Ilat, Ilon, 1, 1, Nlat, Nlon);
                                    dA = source_data.areas.at(area_index);
                                    increment = postprocess_fields.at(Ifield)->at(index) * dA;

                                    if ( source_data.regions.at( source_data.region_names.at(Iregion) ).at(area_index) )
                                    {
                                        field_integrals.at(int_index) += increment;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        for (int_index = 0; int_index < field_integrals.size(); int_index++) {
            reg_area = source_data.region_areas.at(int_index);
            field_averages.at(Ifield).at(int_index) = (reg_area == 0) ? 0. : field_integrals.at(int_index) / reg_area;
        }
    }

    /*
     *  We've silenced standard deviation outputs for now, so don't waste time computing them
     *
     
    // Now that we have region averages, get region standard deviations
    #if DEBUG >= 1
    if (wRank == 0) { fprintf(stdout, "  Computing region standard deviations\n"); }
    #endif
    for (Ifield = 0; Ifield < num_fields; ++Ifield) {
        #if DEBUG >= 2
        if (wRank == 0) { fprintf(stdout, "    processing field %d of %d standard deviations\n", Ifield + 1, num_fields); }
        #endif
        for (Iregion = 0; Iregion < num_regions; ++Iregion) {
            for (Itime = 0; Itime < Ntime; ++Itime) {
                for (Idepth = 0; Idepth < Ndepth; ++Idepth) {

                    int_index = Index(0, Itime, Idepth, Iregion, 1, Ntime, Ndepth, num_regions);
                    reg_area = source_data.region_areas.at(int_index);

                    int_val = 0.;

                    #pragma omp parallel default(none)\
                    private(Ilat, Ilon, index, dA, area_index )\
                    shared( source_data, Ifield, Iregion, Itime, Idepth, int_index,\
                            postprocess_fields, field_averages) \
                    reduction(+ : int_val)
                    { 
                        #pragma omp for collapse(2) schedule(dynamic)
                        for (Ilat = 0; Ilat < Nlat; ++Ilat) {
                            for (Ilon = 0; Ilon < Nlon; ++Ilon) {

                                index = Index(Itime, Idepth, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon);

                                if ( source_data.mask.at(index) ) { // Skip land areas

                                    area_index = Index(0, 0, Ilat, Ilon, 1, 1, Nlat, Nlon);

                                    dA = source_data.areas.at(area_index);

                                    if ( source_data.regions.at( source_data.region_names.at(Iregion) ).at(area_index) )
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

    *
    * 
    */

}
