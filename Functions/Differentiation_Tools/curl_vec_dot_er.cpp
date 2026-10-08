#include "../../constants.hpp"
#include "../../functions.hpp"
#include "../../differentiation_tools.hpp"
#include <algorithm>
#include <vector>
#include <omp.h>
#include <math.h>

void curl_vec_dot_er(
        std::vector<double> & out_arr,
        const std::vector<double> & u_lon,
        const std::vector<double> & u_lat,
        const dataset & source_data,
        const std::vector<short int> & mask
        ) {

    // compute perp_grad_T( cos(lat) * u, cos(lat) * v )
    perp_grad_scalar_T( out_arr, u_lon, u_lat, source_data, mask, true );

    // vort = - (1/cos(lat)) * perp_grad_T( cos(lat) * u, cos(lat) * v )

    const size_t Npts  = u_lon.size();
    size_t index;
    const std::vector<double>   &latitude   = source_data.latitude;

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

            double local_lat = (constants::GRID_TYPE == constants::GridType::MeshGrid) 
                            ? latitude.at(Ilat)
                            : latitude.at(index);

            out_arr[index] = - out_arr[index] / cos(local_lat);
        }
    } // end pragma
}
