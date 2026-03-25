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

    int Ilat_target, Ilon_target;
    #pragma omp parallel \
    default(none) \
    shared( source_data, target_data, vars_to_map ) \
    private( Ilat_target, Ilon_target ) \
    firstprivate( Npts_target,Nlat_target, Nlon_target, Npts_source, Nlat_source, Nlon_source, \
           Nvars, mapping_from_fine_to_coarse )
    { 
        #pragma omp for collapse(2) schedule(static)
        for ( Ilat_target = 0; Ilat_target < Nlat_target; Ilat_target++ ) {
            for ( Ilon_target = 0; Ilon_target < Nlon_target; Ilon_target++ ) {

                double target_lat = target_data.latitude[ Ilat_target];
                double target_lon = target_data.longitude[Ilon_target];

                //
                //// Get nearest lat index in the source
                //
                int Ilat_source_LB, Ilat_source_UB;
                double perc_lat;
                int test_Ilat = std::upper_bound( source_data.latitude.begin(), source_data.latitude.end(), target_lat )
                    - source_data.latitude.begin();
                if ( test_Ilat == 0 ) {
                    // if target point is out of bounds, just use nearest
                    Ilat_source_LB = 0;
                    Ilat_source_UB = 0;
                    perc_lat = 0.;
                } else if ( test_Ilat >= Nlat_source ) {
                    // if target point is out of bounds, just use nearest
                    Ilat_source_LB = Nlat_source - 1;
                    Ilat_source_UB = Nlat_source - 1;
                    perc_lat = 0.;
                } else {
                    // Otherwise, we'll linearly interpolate between lats
                    double Ldist = std::fabs( source_data.latitude[test_Ilat-1] - target_lat );
                    double Rdist = std::fabs( source_data.latitude[test_Ilat]   - target_lat );
                    perc_lat = Ldist / (Ldist+Rdist);
                    Ilat_source_LB = test_Ilat-1;
                    Ilat_source_UB = test_Ilat;
                }

                //
                //// Get nearest lon index
                //
                int Ilon_source_LB, Ilon_source_UB;
                double perc_lon;
                int test_Ilon = std::upper_bound( source_data.longitude.begin(), source_data.longitude.end(), target_lon )
                    - source_data.longitude.begin();
                int Lind, Rind;
                
                // account for periodicity
                if ( test_Ilon == 0 ) {
                    Lind = Nlon_source - 1;
                    Rind = 0;
                } else if ( test_Ilon >= Nlon_source ) {
                    Lind = Nlon_source - 1;
                    Rind = 0;
                } else {
                    Lind = test_Ilon - 1;
                    Rind = test_Ilon;
                }

                // and set up to linearly interpolate
                double Ldist = std::fabs( source_data.longitude[Lind] - target_lon );
                if ( Ldist > M_PI ) { Ldist = 2*M_PI - Ldist; }
                double Rdist = std::fabs( source_data.longitude[Rind] - target_lon );
                if ( Rdist > M_PI ) { Rdist = 2*M_PI - Rdist; }
                perc_lon = Ldist / (Ldist + Rdist);
                Ilon_source_LB = Lind;
                Ilon_source_UB = Rind;


                // Get (source) indices
                size_t source_LL = source_data.local_index( 0, 0, Ilat_source_LB, Ilon_source_LB ),
                       source_UL = source_data.local_index( 0, 0, Ilat_source_UB, Ilon_source_LB ),
                       source_LU = source_data.local_index( 0, 0, Ilat_source_LB, Ilon_source_UB ),
                       source_UU = source_data.local_index( 0, 0, Ilat_source_UB, Ilon_source_UB );

                size_t target_index = target_data.local_index( 0, 0, Ilat_target, Ilon_target );


                // If we're coarsening, determine the new land mask
                if ( mapping_from_fine_to_coarse ) {

                    unsigned int num_mask = 0;

                    num_mask += source_data.mask[source_LL];
                    num_mask += source_data.mask[source_UL];
                    num_mask += source_data.mask[source_LU];
                    num_mask += source_data.mask[source_UU];

                    // if at least two neighbours land, then we're land
                    target_data.mask[target_index] = num_mask >= 2;
                }

                // If we're refining (no land, since working with Psi/Phi)
                // or if the target point is water, compute the interpolation.
                if ( (not(mapping_from_fine_to_coarse)) or (target_data.mask[target_index]) ) {

                    // loop over all vars to process
                    for ( size_t Ivar = 0; Ivar < Nvars; Ivar++ ) {

                        // Linear interpolation between bounding box, in lat-lon coords
                        double val = 0;
                        double val_LL = source_data.variables.at(vars_to_map[Ivar]).at(source_LL),
                               val_UL = source_data.variables.at(vars_to_map[Ivar]).at(source_UL),
                               val_LU = source_data.variables.at(vars_to_map[Ivar]).at(source_LU),
                               val_UU = source_data.variables.at(vars_to_map[Ivar]).at(source_UU);

                        // First, perform latitude interpolation on "west" side
                        double val_L;
                        bool use_L = true;
                        if ( source_data.mask[source_LL] and source_data.mask[source_UL] ) {
                            val_L = val_LL * (1-perc_lat)  +  val_UL * perc_lat;
                        } else if ( source_data.mask[source_LL] ) {
                            val_L = val_LL;
                        } else if ( source_data.mask[source_UL] ) {
                            val_L = val_UL;
                        } else {
                            val_L = 0;
                            use_L = false;
                        }

                        // Then perform latitude interpolation on "east" side
                        double val_U;
                        bool use_U = true;
                        if ( source_data.mask[source_LU] and source_data.mask[source_UU] ) {
                            val_U = val_LU * (1-perc_lat)  +  val_UU * perc_lat;
                        } else if ( source_data.mask[source_LU] ) {
                            val_U = val_LU;
                        } else if ( source_data.mask[source_UU] ) {
                            val_U = val_UU;
                        } else {
                            val_U = 0;
                            use_U = false;
                        }

                        // Then perform the longitude interpolation
                        if ( use_L and use_U ) {
                            val = val_L * (1-perc_lon)  +  val_U * perc_lon;
                        } else if ( use_L ) {
                            val = val_L;
                        } else if ( use_U ) {
                            val = val_U;
                        } else {
                            throw std::runtime_error("Unable to perform interpolation.");
                        }

                        target_data.variables.at(vars_to_map[Ivar]).at(target_index) = val;
                    }
                }

            }
        }

    }

}
