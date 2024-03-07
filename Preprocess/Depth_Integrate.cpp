#include "../constants.hpp"
#include "../functions.hpp"
#include "../netcdf_io.hpp"
#include "../preprocess.hpp"
#include "../differentiation_tools.hpp"
#include <algorithm>
#include <vector>
#include <omp.h>
#include <math.h>
#include <assert.h>
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
    const std::vector<double>   &depth      = source_data.depth;

    const int   Ntime   = source_data.Ntime,    // this is the MPI-local Ntime, not the full Ntime
                Ndepth  = source_data.Ndepth,   // this is the MPI-local Ndepth, not the full Ndepth
                Nlat    = source_data.Nlat,
                Nlon    = source_data.Nlon;

    size_t index, index_below;
    int Itime, Ilat, Ilon, Iz;

    const bool  IS_ELEV = source_data.depth_is_elevation,
                IS_INCR = source_data.depth_is_increasing;


    #if DEBUG >= 0
    fprintf( stdout, "Looping over space to integrate\n");
    #endif
    #pragma omp parallel \
    default(none) \
    shared( field_to_integrate, depth_integral, depth, source_data )\
    private( Itime, Ilat, Ilon, index, index_below, Iz ) \
    firstprivate( Nlon, Nlat, Ndepth, Ntime, IS_INCR, IS_ELEV )
    {
        #pragma omp for collapse(3) schedule(static)
        for ( Itime = 0; Itime < Ntime; ++Itime ) {
            for ( Ilat  = 0; Ilat < Nlat;  ++Ilat  ) {
                for ( Ilon  = 0; Ilon < Nlon;  ++Ilon ) {
                    if ( (IS_ELEV and IS_INCR) or (not(IS_ELEV) and not(IS_INCR)) ) {  
                        // For all points but the bottom, compute the integral.
                        //      the bottom point is assumed to have zero value
                        // Integrate from the bottom up using trapezoid rule
                        index = Index( Itime, 0, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon );
                        depth_integral.at(index) = 0.;
                        for ( Iz = 1; Iz < Ndepth; ++Iz ) {
                            index       = Index( Itime, Iz,   Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon );
                            index_below = Index( Itime, Iz-1, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon );
                            depth_integral.at(index) = depth_integral.at( index_below ) + 
                                ( depth.at(Iz-1) - depth.at(Iz) ) 
                                * ( 0.5 * ( field_to_integrate.at( index ) + field_to_integrate.at( index_below ) ) );
                        }
                    } else {
                        // For all points but the bottom, compute the integral.
                        //      the bottom point is assumed to have zero value
                        // Integrate from the bottom up using trapezoid rule
                        index = Index( Itime, Ndepth-1, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon );
                        depth_integral.at(index) = 0.;
                        for ( Iz = Ndepth - 2; Iz > -1; --Iz ) {
                            index       = Index( Itime, Iz,   Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon );
                            index_below = Index( Itime, Iz+1, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon );
                            depth_integral.at(index) = depth_integral.at( index_below ) + 
                                ( depth.at(Iz+1) - depth.at(Iz) ) 
                                * ( 0.5 * ( field_to_integrate.at( index ) + field_to_integrate.at( index_below ) ) );
                        }
                    }
                }
            }
        }
    }
}
