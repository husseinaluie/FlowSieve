#include "../constants.hpp"
#include "../functions.hpp"
#include "../preprocess.hpp"
#include <algorithm>
#include <vector>
#include <omp.h>
#include <math.h>


void map_grid_to_grid_MeshGrid(
        const dataset & source_data,
        dataset & target_data,
        std::vector<std::string> vars_to_map,
        const MPI_Comm comm
        ) {

    if (not(constants::GRID_TYPE == constants::GridType::MeshGrid)) {
        throw std::runtime_error("map_grid_to_grid_MeshGrid called with a non-MeshGrid grid.");
    }

    const size_t Nvars = vars_to_map.size();

    // Confirm that the vars are in both source and target
    for ( size_t Ivar = 0; Ivar < Nvars; Ivar++ ) {
        if ( not( source_data.variables.count( vars_to_map.at(Ivar) ) ) ) {
            fprintf( stderr, "Source does not contain variable %s.\n", vars_to_map.at(Ivar).c_str() );
            throw std::invalid_argument( "Source does not contain requested variable.\n" );
        }
        if ( not( target_data.variables.count( vars_to_map.at(Ivar) ) ) ) {
            fprintf( stderr, "Target does not contain variable %s.\n", vars_to_map.at(Ivar).c_str() );
            throw std::invalid_argument( "Target does not contain requested variable.\n" );
        }
        #if DEBUG >= 1
        fprintf( stdout, "Variable %'zu has source size %'zu and target size %'zu.\n", 
                Ivar,
                source_data.variables.at( vars_to_map.at(Ivar) ).size(),
                target_data.variables.at( vars_to_map.at(Ivar) ).size()
               );
        #endif
    }

    const size_t Npts_target = target_data.mask.size();
    const size_t Nlat_target = target_data.Nlat;
    const size_t Nlon_target = target_data.Nlon;

    const size_t Npts_source = source_data.mask.size();
    const size_t Nlat_source = source_data.Nlat;
    const size_t Nlon_source = source_data.Nlon;

    const bool mapping_from_fine_to_coarse = ( Npts_source > Npts_target );

    #if DEBUG >= 1
    fprintf( stdout, "Mapping from %'zu points to %'zu.\n", Npts_source, Npts_target );
    #endif

    // we're mapping onto a coarse grid
    //  for this, we'll filter using 2*delta as the
    //  target scale
    const double filter_ell = 2 * ( M_PI * constants::R_earth / std::fmin(Nlat_target, Nlat_source) );

    fprintf( stdout, "Mapping between grids will coarsen at ell = %.0f km\n", filter_ell/1e3 );

    int Ilat_target, Ilon_target;
    #pragma omp parallel \
    default(none) \
    shared( source_data, target_data, vars_to_map ) \
    private( Ilat_target, Ilon_target ) \
    firstprivate( Npts_target,Nlat_target, Nlon_target, Npts_source, Nlat_source, Nlon_source, \
            Nvars, filter_ell, mapping_from_fine_to_coarse )
    { 
        #pragma omp for collapse(2) schedule(static)
        for ( Ilat_target = 0; Ilat_target < Nlat_target; Ilat_target++ ) {
            for ( Ilon_target = 0; Ilon_target < Nlon_target; Ilon_target++ ) {

                double target_lat = target_data.latitude[ Ilat_target];
                double target_lon = target_data.longitude[Ilon_target];

                //
                //// Get nearest lat index in the source
                //
                int centre_Ilat = std::upper_bound( source_data.latitude.begin(), source_data.latitude.end(), target_lat )
                    - source_data.latitude.begin();
                if ( centre_Ilat >= Nlat_source ) { centre_Ilat = Nlat_source - 1; }

                //
                //// Get nearest lon index
                //
                //int centre_Ilon = std::upper_bound( source_data.longitude.begin(), source_data.longitude.end(), target_lon )
                //    - source_data.longitude.begin();
                //if ( centre_Ilon >= Nlon_source ) { centre_Ilon = Nlon_source - 1; }


                // We're fixing the scale at 2 delta, so just go +/- two latitude indices
                std::vector<double> numer_accum( Nvars );
                const int del_lat_logical = 
                    mapping_from_fine_to_coarse ?
                        (int) ceil( (double)Nlat_source / Nlat_target )
                        : 2; // if source is the coarse grid, just do a couple lats
                double denom_accum = 0., water_area_accum = 0.;
                for ( int source_Ilat = std::fmax( centre_Ilat-2*del_lat_logical, 0 );
                        source_Ilat < std::fmin( centre_Ilat+2*del_lat_logical+1, Nlat_source-1 );
                        source_Ilat++
                    ) {
                    for ( int source_Ilon = 0; source_Ilon < Nlon_source; source_Ilon++ ) {
                        size_t source_index = source_data.local_index( 0, 0, source_Ilat, source_Ilon );

                        double dist = source_data.distance_between_datasets_from_IlatIlon( 
                                source_Ilat, source_Ilon, 
                                Ilat_target, Ilon_target,
                                target_data );
                        double kern = kernel(dist, filter_ell);

                        double dA = source_data.areas[source_index];
                        double coeff = dA * kern;

                        for ( int Ivar = 0; Ivar < Nvars; Ivar++ ) {
                            numer_accum[Ivar]+= coeff * source_data.variables.at(vars_to_map[Ivar])[source_index];
                        }

                        denom_accum += coeff;

                        if ( source_data.mask[source_index] ) {
                            water_area_accum += coeff;
                        }

                    }
                }

                size_t target_index = target_data.local_index( 0, 0, Ilat_target, Ilon_target );
                if ( mapping_from_fine_to_coarse ) {
                    //target_data.mask[target_index] = ( water_area_accum >= 0.5 * denom_accum );
                    //target_data.mask[target_index] = ( water_area_accum >= 0.25 * denom_accum );
                    target_data.mask[target_index] = ( water_area_accum >= 0.10 * denom_accum );
                    //target_data.mask[target_index] = ( water_area_accum >= 0.05 * denom_accum );
                    //target_data.mask[target_index] = ( water_area_accum > 0. );
                } else {
                    // if we're refining, don't worry about land
                    //target_data.mask[target_index] = true;
                }
                for ( int Ivar = 0; Ivar < Nvars; Ivar++ ) {
                    target_data.variables.at(vars_to_map[Ivar]).at(target_index) = 
                        ( denom_accum > 0 ) ?
                        numer_accum[Ivar] / denom_accum :
                        0.;
                }

            }
        }
    }

}
