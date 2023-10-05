#include <algorithm>
#include <math.h>
#include <mpi.h>
#include <omp.h>
#include <vector>

#include "../constants.hpp"
#include "../functions.hpp"
#include "../postprocess.hpp"


void compute_zonal_median(
        std::vector<std::vector<double>> & zonal_median,
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
    size_t index, int_index, Ipt;

    std::vector<double> lon_slice;

    #if DEBUG >= 2
    if (wRank == 0) { fprintf( stdout, "Computing zonal medians\n" ); }
    #endif


    #pragma omp parallel default(none)\
    private( Ifield, Ilat, Ilon, Itime, Idepth, index, int_index, Ipt, lon_slice )\
    shared( postprocess_fields, source_data, zonal_median ) \
    firstprivate( Nlon, Nlat, Ndepth, Ntime, num_fields )
    { 
        lon_slice.resize(Nlon);

        #pragma omp for collapse(4) schedule(dynamic)
        for (Ifield = 0; Ifield < num_fields; ++Ifield) {
            for (Itime = 0; Itime < Ntime; ++Itime){
                for (Ilat = 0; Ilat < Nlat; ++Ilat){
                    for (Idepth = 0; Idepth < Ndepth; ++Idepth){

                        std::fill( lon_slice.begin(), lon_slice.end(), 0.);
                        Ipt = 0;

                        for (Ilon = 0; Ilon < Nlon; ++Ilon){
                            index = Index(Itime, Idepth, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon);

                            if ( source_data.mask.at(index) ) {

                                lon_slice[Ipt] = postprocess_fields.at(Ifield)->at(index);
                                Ipt++;

                            }
                        }
                        std::nth_element( lon_slice.begin(), lon_slice.begin() + Ipt/2, lon_slice.begin() + Ipt );
                        int_index = Index( 0, Itime, Idepth, Ilat, 1, Ntime, Ndepth, Nlat );
                        zonal_median[Ifield][int_index] = lon_slice[ Ipt/2 ];
                    }
                }
            }
        }
    }
}
