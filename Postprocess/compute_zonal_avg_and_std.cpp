#include <math.h>
#include <mpi.h>
#include <omp.h>
#include <vector>

#include "../constants.hpp"
#include "../functions.hpp"
#include "../postprocess.hpp"


void compute_zonal_avg_and_std(
        std::vector<std::vector<double>> & zonal_average,
        std::vector<std::vector<double>> & zonal_std_dev,
        const dataset & source_data,
        const std::vector<const std::vector<double>*> & postprocess_fields,
        const MPI_Comm comm
        ){

    int wRank=-1, wSize=-1;
    MPI_Comm_rank( MPI_COMM_WORLD, &wRank );
    MPI_Comm_size( MPI_COMM_WORLD, &wSize );

    const int   Ntime  = source_data.Ntime,
                Ndepth = source_data.Ndepth,
                Nlat   = source_data.Nlat,
                Nlon   = source_data.Nlon;

    const int num_fields = postprocess_fields.size();

    int Ifield, Itime, Idepth, Ilat, Ilon;
    size_t index, area_index, int_index;
    double dA;

    #if DEBUG >= 2
    if (wRank == 0) { fprintf( stdout, "Computing zonal areas\n" ); }
    #endif

    // First, get the zonal areas
    std::vector<double> zonal_areas( Ntime * Ndepth * Nlat, 0. );
    for (Itime = 0; Itime < Ntime; ++Itime) {
        for (Idepth = 0; Idepth < Ndepth; ++Idepth) {
            for (Ilat = 0; Ilat < Nlat; ++Ilat) {
                int_index = Index( Itime, Idepth, Ilat, 0, Ntime, Ndepth, Nlat, 1 );
                for (Ilon = 0; Ilon < Nlon; ++Ilon) {
                    index = Index(Itime, Idepth, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon);

                    // Sum up the area of water cells
                    if ( ( (constants::FILTER_OVER_LAND) and ( source_data.reference_mask.at(index) ) )
                         or
                         ( not(constants::FILTER_OVER_LAND) and ( source_data.mask.at(index) ) )
                       ) {
                        area_index = Index(0, 0, Ilat, Ilon, 1, 1, Nlat, Nlon);
                        dA = source_data.areas.at(area_index);

                        zonal_areas.at( int_index ) += dA;
                    }
                }
            }
        }
    }


    #if DEBUG >= 2
    if (wRank == 0) { fprintf( stdout, "Computing zonal sums\n" ); }
    #endif

    #pragma omp parallel default(none)\
    private( Ifield, Ilat, Ilon, Itime, Idepth, index, int_index, area_index, dA )\
    shared( postprocess_fields, source_data, zonal_average, zonal_areas ) \
    firstprivate( Nlon, Nlat, Ndepth, Ntime, num_fields )
    { 
        #pragma omp for collapse(3) schedule(dynamic)
        for (Itime = 0; Itime < Ntime; ++Itime){
            for (Ilat = 0; Ilat < Nlat; ++Ilat){
                for (Idepth = 0; Idepth < Ndepth; ++Idepth){

                    int_index = Index( 0, Itime, Idepth, Ilat, 1, Ntime, Ndepth, Nlat );

                    for (Ilon = 0; Ilon < Nlon; ++Ilon){

                        // get some indices
                        index = Index(Itime, Idepth, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon);

                        if ( source_data.mask.at(index) ) {

                            area_index = Index(0, 0, Ilat, Ilon, 1, 1, Nlat, Nlon);
                            dA = source_data.areas.at(area_index);

                            for (Ifield = 0; Ifield < num_fields; ++Ifield) {

                                // compute the time average for the part on this processor
                                zonal_average.at(Ifield).at(int_index) 
                                        += dA * postprocess_fields.at(Ifield)->at(index);
                            }
                        }
                    }
                }
            }
        }
    }

    #if DEBUG >= 2
    if (wRank == 0) { fprintf( stdout, "Normalizing to zonal means\n" ); }
    #endif

    // Finally, normalize by zonal area
    for (Itime = 0; Itime < Ntime; ++Itime) {
        for (Idepth = 0; Idepth < Ndepth; ++Idepth) {
            for (Ilat = 0; Ilat < Nlat; ++Ilat) {
                int_index = Index( Itime, Idepth, Ilat, 0, Ntime, Ndepth, Nlat, 1 );

                for (Ifield = 0; Ifield < num_fields; ++Ifield) {
                    if ( zonal_areas.at( int_index ) > 0 ) {
                        zonal_average.at(Ifield).at(int_index) = zonal_average.at(Ifield).at(int_index) / zonal_areas.at( int_index );
                    } else {
                        zonal_average.at(Ifield).at(int_index) = 0.;
                    }
                }
            }
        }
    }
}
