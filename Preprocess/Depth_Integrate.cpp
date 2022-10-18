#include "../constants.hpp"
#include "../functions.hpp"
#include "../netcdf_io.hpp"
#include "../preprocess.hpp"
#include "../differentiation_tools.hpp"
#include <algorithm>
#include <vector>
#include <omp.h>
#include <math.h>
#include "../ALGLIB/stdafx.h"
#include "../ALGLIB/linalg.h"


void depth_integrate(
        std::vector<double> & depth_integral,
        const std::vector<double> & field_to_integrate,
        const dataset & source_data,
        const MPI_Comm comm
        ) {

    int wRank, wSize;
    MPI_Comm_rank( comm, &wRank );
    MPI_Comm_size( comm, &wSize );

    // Create some tidy names for variables
    const std::vector<double>   &time       = source_data.time,
                                &depth      = source_data.depth,
                                &latitude   = source_data.latitude,
                                &longitude  = source_data.longitude;

    const std::vector<int>  &myCounts = source_data.myCounts,
                            &myStarts = source_data.myStarts;

    const int   Ntime   = source_data.Ntime,    // this is the MPI-local Ntime, not the full Ntime
                Ndepth  = source_data.Ndepth,   // this is the MPI-local Ndepth, not the full Ndepth
                Nlat    = source_data.Nlat,
                Nlon    = source_data.Nlon;

    const size_t Npts = Ntime * Ndepth * Nlat * Nlon;
    size_t index, int_index, index_below;
    int Itime, Idepth, Ilat, Ilon, Iz;


    #if DEBUG >= 0
    fprintf( stdout, "Looping over space to integrate\n");
    #endif
    #pragma omp parallel \
    default(none) \
    shared( field_to_integrate, depth_integral, depth )\
    private( Itime, Idepth, Ilat, Ilon, index, index_below, Iz )
    {
        #pragma omp for collapse(3) schedule(static)
        for ( Itime = 0; Itime < Ntime; ++Itime ) {
            for ( Ilat  = 0; Ilat < Nlat;  ++Ilat  ) {
                for ( Ilon  = 0; Ilon < Nlon;  ++Ilon ) {

                    // For all points but the bottom, compute the integral.
                    //      the bottom point is assumed to have zero value
                    // Integrate from the bottom up using trapezoid rule
                    index = Index( Itime, Ndepth-1, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon );
                    depth_integral.at(index) = 0.;
                    for ( Iz = Ndepth - 2; Iz > -1; --Iz ) {
                        index       = Index( Itime, Iz,   Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon );
                        index_below = Index( Itime, Iz+1, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon );
                        depth_integral.at(index) = depth_integral.at( index_below ) + 
                              ( depth.at(Iz) - depth.at(Iz+1) ) 
                            * ( 0.5 * ( field_to_integrate.at( index ) + field_to_integrate.at( index_below ) ) );
                    }
                }
            }
        }
    }

}
