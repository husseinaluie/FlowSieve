#include "../constants.hpp"
#include "../functions.hpp"
#include "../preprocess.hpp"
#include "../differentiation_tools.hpp"
#include <algorithm>
#include <vector>
#include <omp.h>
#include <math.h>
#include <random>
#include <deque>
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>

// Pragma-magic to allow reduction over vector operator
//      thanks to: https://stackoverflow.com/questions/43168661/openmp-and-reduction-on-stdvector
#pragma omp declare reduction(vec_double_plus : std::vector<double> : \
                              std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<double>())) \
                    initializer(omp_priv = decltype(omp_orig)(omp_orig.size()))
#pragma omp declare reduction(vec_int_plus : std::vector<int> : \
                              std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<int>())) \
                    initializer(omp_priv = decltype(omp_orig)(omp_orig.size()))

void map_grid_to_grid(
        const dataset & source_data,
        dataset & target_data,
        std::vector<std::string> vars_to_map,
        const MPI_Comm comm
        ) {


    const size_t Nvars = vars_to_map.size();
    const bool mapping_from_fine_to_coarse = 
        ( source_data.latitude.size() > target_data.latitude.size() );

    // Confirm that vars are in both source and target
    for ( size_t Ivar = 0; Ivar < Nvars; Ivar++ ) {
        if ( not( source_data.variables.count( vars_to_map[Ivar] ) ) ) {
            fprintf( stderr, "Source does not contain variable %s.\n", vars_to_map[Ivar].c_str() );
            throw std::invalid_argument( "Source does not contain requested variable.\n" );
        }
        if ( not( target_data.variables.count( vars_to_map[Ivar] ) ) ) {
            fprintf( stderr, "Target does not contain variable %s.\n", vars_to_map[Ivar].c_str() );
            throw std::invalid_argument( "Source does not contain requested variable.\n" );
        }
    }

    const size_t Nlatlon_target = target_data.mask.size();
    const size_t Nlatlon_source = source_data.mask.size();

    #if DEBUG >= 1
    fprintf( stdout, "Mapping from %'zu points to %'zu.\n", Nlatlon_source, Nlatlon_target );
    #endif

    size_t source_index, target_index;

    if (constants::GRID_TYPE == constants::GridType::MeshGrid) {
        // The Mesh and Unstructured grid routines are different enough that we'll just
        // wrap the MeshGrid version in a different function
        map_grid_to_grid_MeshGrid( source_data, target_data, vars_to_map, comm );
        return;
    }

    const double typical_spacing = ( mapping_from_fine_to_coarse ) 
        ? sqrt( 4 * M_PI * pow(constants::R_earth,2.) / Nlatlon_target )
        : sqrt( 4 * M_PI * pow(constants::R_earth,2.) / Nlatlon_source );

    // First, we need to find the nearest point in source for each target
    std::vector<size_t> nearest_source(Nlatlon_target, Nlatlon_source);
    double max_dist = 0;
    #pragma omp parallel default(none) private( target_index ) \
    shared( source_data, target_data, nearest_source ) \
    firstprivate( Nlatlon_source, Nlatlon_target, typical_spacing, Nvars ) \
    reduction( max:max_dist )
    {
        max_dist = 0.;

        #pragma omp for collapse(1) schedule(guided)
        for ( target_index = 0; target_index < Nlatlon_target; target_index++ ) {

            double target_lat = target_data.latitude.at(target_index),
                   target_lon = target_data.longitude.at(target_index);

            double best_dist = 2 * M_PI * constants::R_earth;

            size_t best_index = 0;

            // start with 1000 random guesses to get a less terrible starting guess
            // Make a random number generator
            std::random_device random_device; // obtain a random number from hardware
            std::mt19937_64 random_generator(random_device()); // seed the generator
            std::uniform_int_distribution<size_t> random_index(0, Nlatlon_source-1); // define the range
            for ( int II = 0; II < 1000; II++ ) {
                size_t source_index = random_index(random_generator);

                double source_lat = source_data.latitude.at(source_index),
                       source_lon = source_data.longitude.at(source_index);

                double test_dist = distance( source_lon, source_lat, target_lon, target_lat );

                if ( test_dist < best_dist ) {
                    best_dist = test_dist;
                    best_index = source_index;
                }
            }

            // also check our neighbours, see if they already have their nearest source set
            for ( size_t II = 0; II < target_data.num_neighbours; II++ ) {
                size_t neighbour_index = target_data.adjacency_indices.at(target_index).at(II);

                if ( nearest_source[neighbour_index] < Nlatlon_source ) {
                    size_t source_index = nearest_source[neighbour_index];

                    double source_lat = source_data.latitude.at(source_index),
                           source_lon = source_data.longitude.at(source_index);

                    double test_dist = distance( source_lat, source_lat, target_lon, target_lat );

                    if ( test_dist < best_dist ) {
                        best_dist = test_dist;
                        best_index = source_index;
                    }
                }
            }

            // and finally loop over space
            for ( size_t source_index = 0; source_index < Nlatlon_source; source_index++ ) {
                double source_lat = source_data.latitude.at(source_index),
                       source_lon = source_data.longitude.at(source_index);

                // short-circuit quickly if obviously too far
                double lat_dist = constants::R_earth * std::fabs( source_lat - target_lat );
                if ( lat_dist > best_dist ) { continue; }

                double del_lon = std::fabs( source_lon - target_lon );
                if ( del_lon > M_PI ) { del_lon -= 2 * M_PI; }
                double approx_cos = cos( std::fmax( std::fabs(target_lat), std::fabs(source_lat) ) );
                double lon_dist = constants::R_earth * del_lon * approx_cos;
                if ( lon_dist > best_dist ) { continue; }

                // otherwise compute the actual distance
                double test_dist = distance( source_lon, source_lat, target_lon, target_lat );

                if ( test_dist < best_dist ) {
                    best_dist = test_dist;
                    best_index = source_index;
                }
            }

            nearest_source[target_index] = best_index;
            max_dist = std::fmax( max_dist, best_dist );

        }
    }

    // To ensure a smooth mapping, we're just going to coarse grain across grids
    // at 10 times the 'typical' spacing of the coarse grid.
    // This will help to ensure a well-behaved mapping, which is more
    // important for seeding than a rough-but-accurate mapping
    const double filter_scale = 2 * 5 * std::fmin(typical_spacing, max_dist);


    #if DEBUG >= 0
    fprintf( stdout, "Greatest distance between source and target is %.2ekm. Typical spacing set as %.2ekm and filter scale at %.2ekm.\n",
          max_dist/1e3, typical_spacing/1e3, filter_scale/1e3 );
    #endif


    // Now that we have the mapping from target to the nearest source, we can proceed with filtering
    #pragma omp parallel \
    default(none) \
    shared( source_data, target_data, vars_to_map, nearest_source ) \
    private( target_index, source_index ) \
    firstprivate( Nlatlon_target, Nlatlon_source, Nvars, filter_scale, mapping_from_fine_to_coarse, stdout )
    {

        #pragma omp for collapse(1) schedule(dynamic)
        for ( target_index = 0; target_index < Nlatlon_target; target_index++) {
            std::vector<bool>   was_rejected(Nlatlon_source, false), 
                                was_accepted(Nlatlon_source, false), 
                                planned_for_testing(Nlatlon_source, false);
            // intentionally using the secretly-a-bitset vector<bool>. bitset itself
            // doesn't allow dynamic sizing

            std::vector<double> value_sums( Nvars, 0. );

            double target_lat = target_data.latitude.at(  target_index );
            double target_lon = target_data.longitude.at( target_index );

            // Seed using info from the nearest source
            source_index = nearest_source[target_index];
            double source_lat = source_data.latitude.at(  source_index ),
                   source_lon = source_data.longitude.at( source_index ),
                   local_dist = distance( target_lon, target_lat, source_lon, source_lat ),
                   kern_val = kernel( local_dist, filter_scale );

            double dArea = source_data.areas.at( source_index ),
                   kernel_normalization = kern_val * dArea,
                   land_area = 0, 
                   water_area = 0;

            #if DEBUG >= 2
            if ( kern_val < 0.99 ) {
                fprintf( stdout, "Nearest to (%.1f, %.1f) is (%.1f, %.1f)), with kernel value of %.2e.\n",
                        target_lat, target_lon, source_lat, source_lon, kern_val  );
            }
            #endif

            if (source_data.mask[source_index]) {
                for (size_t Ivar = 0; Ivar < Nvars; Ivar++) {
                    double local_val = source_data.variables.at(vars_to_map[Ivar]).at(source_index);
                    value_sums.at(Ivar) = kern_val * dArea * local_val;
                    // not +=, here we re-set the values to the local value
                    // before looping to accumulate over the kernel
                }
                water_area += kern_val * dArea;
            } else {
                land_area += kern_val * dArea;
            }

            // Next, seed the 'points to test' with the adjacent points
            std::deque<size_t> points_to_test;
            for (size_t II = 0; II < source_data.num_neighbours; II++ ) {
                size_t neighbour_index = source_data.adjacency_indices.at(source_index)[II];
                points_to_test.push_back( neighbour_index );
                planned_for_testing[ neighbour_index ] = true;
            }

            // So long as we still have points to test, keep testing!
            // This implements a breadth-first search through the adjacency matrix to build
            //      filtering kernel. If a point is within distance, add it to the kernel,
            //      and then test it's neighbours. If those are in, test their neighbours,
            //      and so on. If a point is too far away, we do not test it's neighbours.
            // This then assumes that for any point X within distance L of Y, that there
            //      is an adjacency path from Y to X strictly using points within distance
            //      L of Y.
            // Along the way, test for double inclusing to make sure that we don't test 
            //      points repeatedly etc.
            while ( points_to_test.size() > 0 ) {

                // Pull out the most-recently-added point, and remove it from the 'to test' list,
                // since we're testing it now.
                // Since we're pulling out the most-recently-added, that effectively makes this a
                // depth-first search to build the kernel.
                source_index = points_to_test.front();
                points_to_test.pop_front();
                planned_for_testing[ source_index ] = false;

                source_lat = source_data.latitude.at(  source_index );
                source_lon = source_data.longitude.at( source_index );
                local_dist = distance( target_lon, target_lat, source_lon, source_lat );

                kern_val = kernel( local_dist, filter_scale );

                if ( ( kern_val > 1e-10 ) or ( local_dist <= 0.5 * filter_scale )  ) {
                    was_accepted[source_index] = true;

                    // Accumulate the coarse values
                    dArea = source_data.areas.at(source_index);
                    kernel_normalization += kern_val * dArea;
                    if (source_data.mask[source_index]) {
                        for ( size_t Ivar = 0; Ivar < Nvars; Ivar++) {
                            double local_val = source_data.variables.at(vars_to_map[Ivar]).at(source_index);
                            value_sums.at(Ivar) += kern_val * dArea * local_val;
                        }
                        water_area += kern_val * dArea;
                    } else {
                        land_area += kern_val * dArea;
                    }

                    // Since this point is in the kernel
                    // add its neighbours to the 'to search' list
                    for ( size_t II = 0; II < source_data.num_neighbours; II++ ) {
                        size_t neighbour_index = source_data.adjacency_indices.at(source_index)[II];

                        // but first check if that neighbour has been rejected already
                        if ( was_rejected[neighbour_index] ) { continue; }

                        // then check if that neighbour is already accepted
                        if ( was_accepted[neighbour_index] ) { continue; }

                        // then check if that neighbour is already on the search list
                        if ( planned_for_testing[neighbour_index] ) { continue; }

                        // if it's not on any of those list already
                        // then add it to the 'points to test'
                        points_to_test.push_back( neighbour_index );
                        planned_for_testing[neighbour_index] = true;

                    }
                } else {
                    // Otherwise, record this point as rejected, and move on
                    was_rejected[source_index] = true;
                }

            } // loop through to make the kernel

            // Store the filtered values in the appropriate arrays
            if ( mapping_from_fine_to_coarse ) {
                // Only set the mask if we're coarsening
                //target_data.mask[target_index] = source_data.mask[ nearest_source[target_index] ];
                target_data.mask[target_index] = ( water_area > 0 );
                //target_data.mask[target_index] = ( water_area > land_area );
            }
            for ( size_t Ivar = 0; Ivar < Nvars; Ivar++ ) {
                target_data.variables.at(vars_to_map[Ivar]).at(target_index) = 
                    ( kernel_normalization == 0 )
                    ?
                    0.
                    :
                    value_sums[Ivar] / kernel_normalization;
            }
        }
    }
}
