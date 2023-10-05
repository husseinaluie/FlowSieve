#include <math.h>
#include <mpi.h>
#include <omp.h>
#include <vector>
#include <algorithm>

#include "../constants.hpp"
#include "../functions.hpp"
#include "../postprocess.hpp"


void compute_coarsened_map(
        std::vector< std::vector< double > > & coarsened_maps,
        const dataset & source_data,
        const std::vector<const std::vector<double>*> & postprocess_fields,
        const MPI_Comm comm
        ) {

    const int   Ntime   = source_data.Ntime,
                Ndepth  = source_data.Ndepth,
                Nlat    = source_data.Nlat,
                Nlon    = source_data.Nlon;

    const std::vector<double>   &latitude   = source_data.latitude,
                                &longitude  = source_data.longitude,
                                &coarse_latitude = source_data.coarse_map_lat,
                                &coarse_longitude = source_data.coarse_map_lon;

    const int   Nlat_coarse = source_data.coarse_map_lat.size(),
                Nlon_coarse = source_data.coarse_map_lon.size();

    const int   num_fields = postprocess_fields.size();

    double dA, dA_coarse, increment;

    int Ifield, Itime, Idepth, Ilat, Ilon, Ilat_coarse, Ilon_coarse,
        lat_LB, lat_UB, lon_LB, lon_UB;
    size_t coarse_index, index;
    bool is_water;

    std::vector<double> coarsened_areas( postprocess_fields[0]->size(), 0. );

    #pragma omp parallel default(none)\
    private(Ilat, Ilon, index, coarse_index, dA, is_water, \
            Idepth, Itime, lat_LB, lat_UB, lon_LB, lon_UB )\
    shared( source_data, latitude, longitude, coarse_latitude, coarse_longitude, coarsened_areas ) \
    firstprivate( Nlon, Nlat, Ndepth, Ntime, Nlat_coarse, Nlon_coarse )
    { 
        #pragma omp for collapse(4) schedule(static)
        for (Itime = 0; Itime < Ntime; ++Itime) {
            for (Idepth = 0; Idepth < Ndepth; ++Idepth) {
                for (Ilat_coarse = 0; Ilat_coarse < Nlat_coarse; ++Ilat_coarse) {
                    for (Ilon_coarse = 0; Ilon_coarse < Nlon_coarse; ++Ilon_coarse) {

                        lat_LB = std::lower_bound( latitude.begin(), latitude.end(), coarse_latitude.at(Ilat_coarse) ) - latitude.begin();
                        lat_UB = ( Ilat_coarse + 1 < Nlat_coarse )
                                 ? ( std::lower_bound( latitude.begin(), latitude.end(), coarse_latitude.at(Ilat_coarse+1) ) ) - latitude.begin()
                                 : Nlat - 1;

                        lat_LB = ( lat_LB < 0 ) ? 0 : ( lat_LB >= Nlat ) ? Nlat - 1 : lat_LB;
                        lat_UB = ( lat_UB < 0 ) ? 0 : ( lat_UB >= Nlat ) ? Nlat - 1 : lat_UB;

                        lon_LB = std::lower_bound( longitude.begin(), longitude.end(), coarse_longitude.at(Ilon_coarse) ) - longitude.begin();
                        lon_UB = ( Ilon_coarse + 1 < Nlon_coarse )
                                 ? ( std::lower_bound( longitude.begin(), longitude.end(), coarse_longitude.at(Ilon_coarse+1) ) ) - longitude.begin()
                                 : Nlon - 1;

                        lon_LB = ( lon_LB < 0 ) ? 0 : ( lon_LB >= Nlon ) ? Nlon - 1 : lon_LB;
                        lon_UB = ( lon_UB < 0 ) ? 0 : ( lon_UB >= Nlon ) ? Nlon - 1 : lon_UB;

                        coarse_index = Index(Itime, Idepth, Ilat_coarse, Ilon_coarse, Ntime, Ndepth, Nlat_coarse, Nlon_coarse);

                        for ( Ilat = lat_LB; Ilat < lat_UB; Ilat++ ) {
                            for ( Ilon = lon_LB; Ilon < lon_UB; Ilon++ ) {
                                index = Index(Itime, Idepth, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon);
                                is_water = source_data.mask.at(index);

                                if (is_water) {
                                    index = Index(0, 0, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon);
                                    dA = source_data.areas.at(index);
                                    coarsened_areas.at(coarse_index) += dA;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    #if DEBUG >= 1
    int wRank;
    MPI_Comm_rank( comm, &wRank );
    #endif

    #if DEBUG >= 1
    if (wRank == 0) { fprintf(stdout, "  Computing coarsened maps\n"); }
    fflush(stdout);
    #endif
    for (Ifield = 0; Ifield < num_fields; ++Ifield) {
        #if DEBUG >= 2
        if (wRank == 0) { fprintf(stdout, "    processing field %d of %d means\n", Ifield + 1, num_fields); }
        fflush(stdout);
        #endif

        #pragma omp parallel default(none)\
        private(Ilat, Ilon, index, coarse_index, dA, dA_coarse, increment, is_water, \
                Idepth, Itime, lat_LB, lat_UB, lon_LB, lon_UB )\
        shared( source_data, Ifield, coarsened_maps, postprocess_fields, latitude, longitude, \
                coarse_latitude, coarse_longitude, coarsened_areas ) \
        firstprivate( Nlon, Nlat, Ndepth, Ntime, Nlat_coarse, Nlon_coarse )
        { 
            #pragma omp for collapse(4) schedule(static)
            for (Itime = 0; Itime < Ntime; ++Itime) {
                for (Idepth = 0; Idepth < Ndepth; ++Idepth) {
                    for (Ilat_coarse = 0; Ilat_coarse < Nlat_coarse; ++Ilat_coarse) {
                        for (Ilon_coarse = 0; Ilon_coarse < Nlon_coarse; ++Ilon_coarse) {

                            lat_LB = std::lower_bound( latitude.begin(), latitude.end(), coarse_latitude.at(Ilat_coarse) ) - latitude.begin();
                            lat_UB = ( Ilat_coarse + 1 < Nlat_coarse )
                                     ? ( std::lower_bound( latitude.begin(), latitude.end(), coarse_latitude.at(Ilat_coarse+1) ) ) - latitude.begin()
                                     : Nlat - 1;

                            lat_LB = ( lat_LB < 0 ) ? 0 : ( lat_LB >= Nlat ) ? Nlat - 1 : lat_LB;
                            lat_UB = ( lat_UB < 0 ) ? 0 : ( lat_UB >= Nlat ) ? Nlat - 1 : lat_UB;

                            lon_LB = std::lower_bound( longitude.begin(), longitude.end(), coarse_longitude.at(Ilon_coarse) ) - longitude.begin();
                            lon_UB = ( Ilon_coarse + 1 < Nlon_coarse )
                                     ? ( std::lower_bound( longitude.begin(), longitude.end(), coarse_longitude.at(Ilon_coarse+1) ) ) - longitude.begin()
                                     : Nlon - 1;

                            lon_LB = ( lon_LB < 0 ) ? 0 : ( lon_LB >= Nlon ) ? Nlon - 1 : lon_LB;
                            lon_UB = ( lon_UB < 0 ) ? 0 : ( lon_UB >= Nlon ) ? Nlon - 1 : lon_UB;

                            coarse_index = Index(Itime, Idepth, Ilat_coarse, Ilon_coarse, Ntime, Ndepth, Nlat_coarse, Nlon_coarse);
                            dA_coarse = coarsened_areas.at(coarse_index);

                            for ( Ilat = lat_LB; Ilat < lat_UB; Ilat++ ) {
                                for ( Ilon = lon_LB; Ilon < lon_UB; Ilon++ ) {
                                    index = Index(Itime, Idepth, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon);
                                    is_water = source_data.mask.at(index);

                                    if (is_water) {
                                        index = Index(0, 0, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon);
                                        dA = source_data.areas.at(index);

                                        index = Index(Itime, Idepth, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon);

                                        increment = (dA_coarse > 0) ? (postprocess_fields.at(Ifield)->at(index) * dA / dA_coarse) : 0.;

                                        coarsened_maps.at(Ifield).at( coarse_index ) += increment;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

}
