#include <math.h>
#include <mpi.h>
#include <omp.h>
#include <vector>

#include "../constants.hpp"
#include "../functions.hpp"
#include "../postprocess.hpp"


void compute_time_avg_std(
        std::vector<std::vector<double>> & time_average,
        std::vector<std::vector<double>> & time_std_dev,
        const dataset & source_data,
        const std::vector<const std::vector<double>*> & postprocess_fields,
        const std::vector<int> & mask_count,
        const std::vector<bool> & always_masked,
        const int full_Ntime,
        ){

    //MPI_Comm comm = source_data.MPI_Comm_Global;
    MPI_Comm comm = source_data.MPI_subcomm_samedepths;

    int wRank=-1, wSize=-1;
    MPI_Comm_rank( comm, &wRank );
    MPI_Comm_size( comm, &wSize );

    const int   Ntime  = source_data.Ntime,
                Ndepth = source_data.Ndepth,
                Nlat   = source_data.Nlat,
                Nlon   = source_data.Nlon;

    const int num_fields = postprocess_fields.size();

    const int chunk_size = get_omp_chunksize(Nlat, Nlon);

    int Ifield, Itime, Idepth, Ilat, Ilon;
    size_t index, space_index;

    // storage arrays for local values (before MPI reducing)
    std::vector<std::vector<double>> time_average_loc(num_fields), time_std_dev_loc(num_fields);
    for (Ifield = 0; Ifield < num_fields; ++Ifield) {
        time_average_loc.at(Ifield).resize( Ndepth * Nlat * Nlon, 0. );
        time_std_dev_loc.at(Ifield).resize( Ndepth * Nlat * Nlon, 0. );
    }

    #pragma omp parallel default(none)\
    private(Ifield, Ilat, Ilon, Itime, Idepth, index, space_index )\
    shared(postprocess_fields, source_data, always_masked, mask_count, time_average_loc)
    { 
        #pragma omp for collapse(3) schedule(dynamic, chunk_size)
        for (Ilat = 0; Ilat < Nlat; ++Ilat){
            for (Ilon = 0; Ilon < Nlon; ++Ilon){
                for (Idepth = 0; Idepth < Ndepth; ++Idepth){

                    space_index = Index(0, Idepth, Ilat, Ilon, 1, Ndepth, Nlat, Nlon);

                    if (not(always_masked.at(space_index))) { // Skip land areas
                        for (Itime = 0; Itime < Ntime; ++Itime){

                            // get some indices
                            index = Index(Itime, Idepth, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon);

                            if ( source_data.mask.at(index) ) {
                                for (Ifield = 0; Ifield < num_fields; ++Ifield) {

                                    // compute the time average for the part on this processor
                                    time_average_loc.at(Ifield).at(space_index) 
                                        += postprocess_fields.at(Ifield)->at(index) / mask_count.at(space_index);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // Now communicate with other processors to get the full time average
    #if DEBUG >= 2
    if (wRank == 0) { fprintf(stdout, "  .. .. reducing across processors\n"); }
    fflush(stdout);
    #endif

    for (Ifield = 0; Ifield < num_fields; ++Ifield) {
        MPI_Allreduce(&(time_average_loc.at(Ifield)[0]),
                      &(time_average.at(    Ifield)[0]),
                      Ndepth * Nlat * Nlon, MPI_DOUBLE, MPI_SUM, comm);
    }

    /*
     *  We've silenced time standard deviation outputs for now, so don't waste time computing them
     *

    // Compute the standard deviation
    #if DEBUG >= 2
    if (wRank == 0) { fprintf(stdout, "  .. computing time standard deviations\n"); }
    fflush(stdout);
    #endif

    #pragma omp parallel default(none)\
    private(Ifield, Ilat, Ilon, Itime, Idepth, index, space_index )\
    shared(wSize, postprocess_fields, source_data, always_masked, mask_count, \
            time_average_loc, time_std_dev_loc, time_average)
    { 
        #pragma omp for collapse(3) schedule(guided, chunk_size)
        for (Ilat = 0; Ilat < Nlat; ++Ilat){
            for (Ilon = 0; Ilon < Nlon; ++Ilon){
                for (Idepth = 0; Idepth < Ndepth; ++Idepth){
                    space_index = Index(0, Idepth, Ilat, Ilon,
                                        1, Ndepth, Nlat, Nlon);
                    if (not(always_masked.at(space_index))) { // Skip land areas
                        for (Itime = 0; Itime < Ntime; ++Itime){

                            // get some indices
                            index = Index(Itime, Idepth, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon);

                            if ( source_data.mask.at(index) ) { // Skip land areas
                                for (Ifield = 0; Ifield < num_fields; ++Ifield) {

                                    // compute the time std. dev. for 
                                    // the part on this processor
                                    time_std_dev_loc.at(Ifield).at(space_index) += 
                                        pow(   postprocess_fields.at(Ifield)->at(index) 
                                             - time_average.at(Ifield).at(space_index), 
                                            2.) / mask_count.at(space_index);
                                }
                            }
                        }
                    } else {
                        for (Ifield = 0; Ifield < num_fields; ++Ifield) {
                            // We'll handle these later to avoid over-flow issues
                            time_std_dev_loc.at(Ifield).at(space_index) = 0.; 
                        }
                    }
                }
            }
        }
    }

    #if DEBUG >= 2
    if (wRank == 0) { fprintf(stdout, "  .. .. reducing across processors and taking square root\n"); }
    fflush(stdout);
    #endif

    // Now communicate with other processors to get the full time average
    for (Ifield = 0; Ifield < num_fields; ++Ifield) {
        MPI_Allreduce(&(time_std_dev_loc.at(Ifield)[0]),
                      &(time_std_dev.at(    Ifield)[0]),
                      Ndepth * Nlat * Nlon, MPI_DOUBLE, MPI_SUM, comm);
    }

    // Re-mask to fill in land areas / apply sqrt to water
    for (int Ifield = 0; Ifield < num_fields; ++Ifield) {
        #pragma omp parallel default(none)\
        private(index, space_index)\
        shared(Ifield, time_average, time_std_dev, always_masked)
        { 
            #pragma omp for collapse(3) schedule(dynamic, chunk_size)
            for (Idepth = 0; Idepth < Ndepth; ++Idepth){
                for (Ilat = 0; Ilat < Nlat; ++Ilat){
                    for (Ilon = 0; Ilon < Nlon; ++Ilon){
                        space_index = Index(0, Idepth, Ilat, Ilon, 1, Ndepth, Nlat, Nlon);
                        if ( not( always_masked.at(space_index) ) ) {
                            time_std_dev.at(Ifield).at(index) = sqrt(time_std_dev.at(Ifield).at(index));
                        }
                    }
                }
            }
        }
    }

    *
    * 
    */

}
