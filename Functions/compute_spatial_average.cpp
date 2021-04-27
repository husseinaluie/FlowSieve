#include "../constants.hpp"
#include "../functions.hpp"
#include "../differentiation_tools.hpp"
#include <algorithm>
#include <vector>
#include <omp.h>
#include <math.h>

/*!
 * \brief Given a field, compute the horizontal average.
 *
 * The result, means, is a function of time and depth.
 *
 *  @param[in,out]  means                       where to store horizontal spatial averages
 *  @param[in]      field                       array of vector to be averaged
 *  @param[in]      areas                       2D array of cell areas
 *  @param[in]      Ntime,Ndepth,Nlat,Nlon      (MPI-local) dimension sizes
 *  @param[in]      mask                        2D array to distinguish land from water
 *
 */
void compute_spatial_average(
        std::vector<double> & means,
        const std::vector<double> & field,
        const std::vector<double> & areas,
        const int Ntime,
        const int Ndepth,
        const int Nlat,
        const int Nlon,
        const std::vector<bool> & mask
        ) {

    const int OMP_chunksize = get_omp_chunksize(Nlat,Nlon);

    double integrated_area=0., integrated_sum=0.;
    int index, sub_index, Itime, Idepth, Ilat, Ilon;

    means.resize(Ntime*Ndepth);

    for (Itime = 0; Itime < Ntime; ++Itime) {
        for (Idepth = 0; Idepth < Ndepth; ++Idepth) {

            integrated_area = 0.;
            integrated_sum  = 0.;

            #pragma omp parallel \
            default(none) \
            shared( field, mask, areas, Itime, Idepth ) \
            private( Ilat, Ilon, index, sub_index )\
            reduction(+ : integrated_area,integrated_sum)
            {
                #pragma omp for collapse(2) schedule(dynamic, OMP_chunksize)
                for (Ilat = 0; Ilat < Nlat; ++Ilat) {
                    for (Ilon = 0; Ilon < Nlon; ++Ilon) {

                        index = Index( Itime, Idepth, Ilat, Ilon,
                                       Ntime, Ndepth, Nlat, Nlon );
                        sub_index = Index( 0, 0, Ilat, Ilon,
                                           1, 1, Nlat, Nlon );

                        if ( mask.at(index) ) { // Skip land areas
                            integrated_area += areas.at(sub_index);
                            integrated_sum  += areas.at(sub_index) * field.at(index);
                        }
                    }
                }
            }
            sub_index = Index( 0, 0, Itime, Idepth,
                               1, 1, Ntime, Ndepth );

            means.at(sub_index) = integrated_sum / integrated_area;
        }
    }
}
