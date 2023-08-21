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

void compute_region_avg_and_std_OkuboWeiss(
        std::vector< std::vector< double > > & field_averages,
        std::vector< std::vector< double > > & field_std_devs,
        std::vector< double > & OkuboWeiss_areas,
        const dataset & source_data,
        const std::vector<const std::vector<double>*> & postprocess_fields,
        const std::vector<double> & OkuboWeiss,
        const std::vector<double> & OkuboWeiss_bounds,
        const int NOkubo,
        const MPI_Comm comm
        ) {

    //fprintf( stdout, " %zu, %zu, %g, %g, %zu\n ", OkuboWeiss.size(), 
    //        OkuboWeiss_bounds.size(), OkuboWeiss_bounds.front(), OkuboWeiss_bounds.back(), source_data.areas.size() );
    //fflush( stdout );

    const int   Ntime   = source_data.Ntime,
                Ndepth  = source_data.Ndepth,
                Nlat    = source_data.Nlat,
                Nlon    = source_data.Nlon;

    const int   num_regions   = source_data.region_names.size(),
                num_fields    = postprocess_fields.size();

    double dA;

    int Ifield, Iregion, Itime, Idepth, Ilat, Ilon, IOkubo;
    size_t int_index, area_index, index;

    #if DEBUG >= 1
    int wRank;
    MPI_Comm_rank( comm, &wRank );
    #endif

    #if DEBUG >= 1
    if (wRank == 0) { fprintf(stdout, "  Computing area OkuboWeiss histograms\n"); }
    #endif

    std::vector<double> area_sums( Ntime * Ndepth * NOkubo * num_regions, 0. );

    //fprintf( stdout, "%d, %d, %d, %d, %d, %d - %'zu\n", 
    //         Ntime, Ndepth, Nlat, Nlon, NOkubo, num_regions, area_sums.size());

    #pragma omp parallel default(none)\
    private( Iregion, Itime, Idepth, Ilat, Ilon, IOkubo, index, area_index, int_index, dA )\
    shared( source_data, OkuboWeiss, OkuboWeiss_bounds ) \
    firstprivate( Nlon, Nlat, Ndepth, Ntime, NOkubo, num_regions ) \
    reduction(vec_double_plus : area_sums)
    { 
        #pragma omp for collapse(5) schedule(static)
        for (Iregion = 0; Iregion < num_regions; ++Iregion) {
            for (Itime = 0; Itime < Ntime; ++Itime) {
                for (Idepth = 0; Idepth < Ndepth; ++Idepth) {
                    for (Ilat = 0; Ilat < Nlat; ++Ilat) {
                        for (Ilon = 0; Ilon < Nlon; ++Ilon) {

                            index = Index(Itime, Idepth, Ilat, Ilon,
                                          Ntime, Ndepth, Nlat, Nlon);

                            if ( constants::FILTER_OVER_LAND or source_data.mask.at(index) ) { // Skip land areas
                                IOkubo = std::lower_bound(  OkuboWeiss_bounds.begin(), 
                                                            OkuboWeiss_bounds.end(), 
                                                            OkuboWeiss.at(index) ) 
                                         - OkuboWeiss_bounds.begin();
                                if (IOkubo >= NOkubo) { IOkubo = NOkubo - 1; }
                                int_index = Index(Itime, Idepth, IOkubo, Iregion,
                                                  Ntime, Ndepth, NOkubo, num_regions);
                                area_index = Index(0, 0, Ilat, Ilon, 1, 1, Nlat, Nlon);

                                if ( source_data.regions.at( source_data.region_names.at(Iregion) ).at(area_index) ) {
                                    dA = source_data.areas.at(area_index);
                                    area_sums.at(int_index) += dA;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    #pragma omp parallel default(none) \
    private( int_index ) \
    shared( OkuboWeiss_areas, area_sums ) \
    firstprivate( Ifield )
    {
        #pragma omp for collapse(1) schedule(static)
        for (int_index = 0; int_index < area_sums.size(); int_index++) {
            OkuboWeiss_areas.at(int_index) = area_sums.at(int_index);
        }
    }

    //
    //// Now apply to the actual data
    //

    #if DEBUG >= 1
    if (wRank == 0) { fprintf(stdout, "  Computing region OkuboWeiss histograms\n"); }
    #endif

    for (Ifield = 0; Ifield < num_fields; ++Ifield) {
        #if DEBUG >= 2
        if (wRank == 0) { fprintf(stdout, "    processing field %d of %d\n", Ifield + 1, num_fields); }
        #endif

        std::vector<double> field_integrals( Ntime * Ndepth * NOkubo * num_regions, 0. );

        #pragma omp parallel default(none)\
        private( Iregion, Itime, Idepth, Ilat, Ilon, IOkubo,\
                index, area_index, int_index, dA )\
        shared( source_data, postprocess_fields, OkuboWeiss, OkuboWeiss_bounds ) \
        firstprivate( Ifield, Nlon, Nlat, Ndepth, Ntime, NOkubo, num_regions ) \
        reduction(vec_double_plus : field_integrals)
        { 
            #pragma omp for collapse(5) schedule(static)
            for (Iregion = 0; Iregion < num_regions; ++Iregion) {
                for (Itime = 0; Itime < Ntime; ++Itime) {
                    for (Idepth = 0; Idepth < Ndepth; ++Idepth) {
                        for (Ilat = 0; Ilat < Nlat; ++Ilat) {
                            for (Ilon = 0; Ilon < Nlon; ++Ilon) {

                                index = Index(Itime, Idepth, Ilat, Ilon,
                                              Ntime, Ndepth, Nlat, Nlon);

                                if ( constants::FILTER_OVER_LAND or source_data.mask.at(index) ) { // Skip land areas
                                    IOkubo = std::lower_bound(  OkuboWeiss_bounds.begin(), 
                                                                OkuboWeiss_bounds.end(), 
                                                                OkuboWeiss.at(index) ) 
                                             - OkuboWeiss_bounds.begin();
                                    if (IOkubo >= NOkubo) { IOkubo = NOkubo - 1; }
                                    int_index = Index(Itime, Idepth, IOkubo, Iregion,
                                                      Ntime, Ndepth, NOkubo, num_regions);
                                    area_index = Index(0, 0, Ilat, Ilon, 1, 1, Nlat, Nlon);

                                    if ( source_data.regions.at( source_data.region_names.at(Iregion) ).at(area_index) ) {
                                        dA = source_data.areas.at(area_index);
                                        field_integrals.at(int_index) += postprocess_fields.at(Ifield)->at(index) * dA;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        #if DEBUG >= 2
        if (wRank == 0) { fprintf(stdout, "    copying results for field %d of %d to outputs\n", Ifield + 1, num_fields); }
        #endif

        #pragma omp parallel default(none) \
        private( int_index ) \
        shared( field_integrals, field_averages ) \
        firstprivate( Ifield )
        {
            #pragma omp for collapse(1) schedule(static)
            for (int_index = 0; int_index < field_integrals.size(); int_index++) {
                field_averages.at(Ifield).at(int_index) = field_integrals.at(int_index);
            }
        }
    }
}
