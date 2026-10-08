#include "../../constants.hpp"
#include "../../functions.hpp"
#include "../../differentiation_tools.hpp"
#include <algorithm>
#include <vector>
#include <omp.h>
#include <math.h>

void grad_scalar_T(
        std::vector<double> & out_arr,
        const std::vector<double> & u_lon,
        const std::vector<double> & u_lat,
        const dataset & source_data,
        const std::vector<short int> & mask,
        const bool scale_vels_by_cos
        ) {

    #if DEBUG>= 2
    fprintf( stdout, "entering grad_scalar_T\n" );
    #endif

    // ret = ddlon(vel_lon) / cos_lat + ddlat( u_lat * cos_lat ) / cos_lat 

    const size_t Npts  = u_lon.size();

    size_t index;
    const std::vector<double>   &latitude   = source_data.latitude,
                                &longitude  = source_data.longitude;

    // Store u_lon / cos_lat in an array
    std::vector<double> scaled_vel( Npts, 0. );

    #if DEBUG>= 2
    fprintf( stdout, "   scaling the velocities\n" );
    #endif

    #pragma omp parallel \
    default(none) \
    shared( latitude, source_data, u_lon, u_lat, scaled_vel )\
    private( index ) firstprivate( Npts, scale_vels_by_cos )
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

            #ifdef NDEBUG
            if ( scale_vels_by_cos ) {
                scaled_vel[index] = u_lat[index] * cos_lat;
            } else {
                scaled_vel[index] = u_lon[index] / cos_lat;
            }
            #else
            if ( scale_vels_by_cos ) {
                scaled_vel.at(index) = u_lat.at(index) * cos_lat;
            } else {
                scaled_vel.at(index) = u_lon.at(index) / cos_lat;
            }
            #endif
        }
    } // end pragma

    std::vector<const std::vector<double>*> deriv_fields = { 
        scale_vels_by_cos ? &u_lon : &scaled_vel, 
        scale_vels_by_cos ? &scaled_vel : &u_lat 
    };

    #if DEBUG>= 2
    fprintf( stdout, "   computing derivatives\n" );
    #endif
    
    #pragma omp parallel \
    default(none) \
    shared( out_arr, latitude, longitude, mask, source_data, deriv_fields, stdout )\
    private( index ) \
    firstprivate( Npts )
    {

        double dulon_dlon, dulat_dlat;

        std::vector<double*> lon_deriv_vals = { &dulon_dlon, NULL }, 
                             lat_deriv_vals = { NULL, &dulat_dlat };

        #pragma omp for collapse(1) schedule(guided)
        for (index = 0; index < Npts; ++index) {

            int Itime, Idepth, Ilat, Ilon;

            if (constants::GRID_TYPE == constants::GridType::MeshGrid) {
                source_data.index1to4_local(index, Itime, Idepth, Ilat, Ilon);

                longitude_derivative_transpose_at_point(
                        lon_deriv_vals, deriv_fields,
                        source_data, Itime, Idepth, Ilat, Ilon, mask);

                latitude_derivative_transpose_at_point(
                        lat_deriv_vals, deriv_fields,
                        source_data, Itime, Idepth, Ilat, Ilon, mask, true);
            } else {
                longitude_derivative_transpose_at_point(
                        lon_deriv_vals, deriv_fields, 
                        source_data, 0, 0, index, index, mask);

                latitude_derivative_transpose_at_point(
                        lat_deriv_vals, deriv_fields,
                        source_data, 0, 0, index, index, mask);
            }

            double tmp = (dulon_dlon + dulat_dlat) / constants::R_earth;
            //fprintf( stdout, "    calculated val = %f\n", tmp );

            #ifndef NDEBUG
            if ( isnan(tmp) ) { throw std::runtime_error("NaN value detected.\n"); }
            #endif

            #ifdef NDEBUG
            out_arr[index] = tmp;
            #else
            out_arr.at(index) = tmp;
            #endif

            //fprintf( stdout, "    written to index %'zu\n", index );

        }
    } // end pragma

    #if DEBUG>= 2
    fprintf( stdout, "exiting grad_scalar_T\n" );
    #endif
}
