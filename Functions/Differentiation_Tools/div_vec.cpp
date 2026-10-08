#include "../../constants.hpp"
#include "../../functions.hpp"
#include "../../differentiation_tools.hpp"
#include <algorithm>
#include <vector>
#include <omp.h>
#include <math.h>

void div_vec(
        std::vector<double> & out_arr,
        const std::vector<double> & u_lon,
        const std::vector<double> & u_lat,
        const dataset & source_data,
        const std::vector<short int> & mask
        ) {

    // compute grad_T( cos(lat) * u, cos(lat) * v )
    grad_scalar_T( out_arr, u_lon, u_lat, source_data, mask, true );

    // div = -(1/cos(lat)) * grad_T( cos(lat) * u, cos(lat) * v )

    const size_t Npts  = u_lon.size();

    size_t index;
    const std::vector<double> &latitude = source_data.latitude;

    #pragma omp parallel \
    default(none) \
    shared( latitude, source_data, out_arr )\
    private( index ) firstprivate( Npts )
    {
        #pragma omp for collapse(1) schedule(guided)
        for (index = 0; index < Npts; ++index) {

            int Itime, Idepth, Ilat, Ilon;
            if (constants::GRID_TYPE == constants::GridType::MeshGrid) {
                source_data.index1to4_local(index, Itime, Idepth, Ilat, Ilon);
            }

            const bool pre_comp = 
                constants::PRECOMPUTE_SIN_COS_LATITUDE 
                and source_data.sin_cos_lat_pre_comp;

            double cos_lat;
            if ( pre_comp ) {
                #ifdef NDEBUG
                cos_lat = source_data.cos_lat[Ilat];
                #else
                cos_lat = source_data.cos_lat.at(Ilat);
                #endif
            } else {
                #ifdef NDEBUG
                double local_lat = (constants::GRID_TYPE == constants::GridType::MeshGrid) 
                    ? latitude[Ilat]
                    : latitude[index];
                #else
                double local_lat = (constants::GRID_TYPE == constants::GridType::MeshGrid) 
                    ? latitude.at(Ilat)
                    : latitude.at(index);
                #endif

                cos_lat = cos(local_lat);
            }

            out_arr.at(index) = -1. * out_arr.at(index) / cos_lat;
        }
    } // end pragma
}
