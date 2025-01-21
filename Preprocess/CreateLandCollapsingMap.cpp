#include "../constants.hpp"
#include "../functions.hpp"
#include "../preprocess.hpp"
#include "../differentiation_tools.hpp"
#include <algorithm>
#include <deque>
#include <vector>
#include <omp.h>
#include <math.h>

// Create land-collapsing map [i.e. map continguous land to single point]
void HelmholtzDataClass::CreateLandCollapsingMap(
        const dataset & data
        ) {

    
    const int Ntime  = data.Ntime,
              Ndepth = data.Ndepth,
              Nlat   = data.Nlat,
              Nlon   = data.Nlon;
    const size_t Npts = data.mask.size();

    size_t num_mapped_onto = 0, num_mapped_points = 0, Ipt, Ineighbour, Itest;
    std::deque<size_t> points_to_test;

    const bool collapse_coast = false;

    pt_maps_to.resize(Npts, Npts);
    num_mapped_before_col.resize( Npts, 0 ); 
    num_mapped_before_row.resize( Npts, 0 );
    num_land_before.resize( Npts, 0 );
    num_mapped_points = 0;
    num_mapped_onto = 0;

    #if DEBUG >= 2
    fprintf( stdout, "CreateLandMap: Npts = %'zu\n", Npts );
    #endif
    size_t num_land = 0;

    if (constants::GRID_TYPE == constants::GridType::MeshGrid) {

        int Itime, Idepth, Ilat, Ilon, Ineigh;
        std::vector<double> diff_vec;

        for (Ipt = 0; Ipt < Npts; Ipt++) {

            if (data.mask[Ipt]) { 
                // Water points map to themselves
                pt_maps_to[Ipt] = Ipt; 
            } else {
                num_land++;
                // Check if we've already mapped this point
                if (pt_maps_to[Ipt] < Npts) { continue; }

                // Begin a depth-first search through the adjacency
                // matrix, mapping all land-only neighbours to this one
                // This point will become the representative for this 'island'
                points_to_test.clear();
                pt_maps_to[Ipt] = Ipt; 
                num_mapped_onto++;


                // Get coordinates
                Index1to4( Ipt, Itime, Idepth, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon );

                // Add neighbours to to-test list
                Ineigh = ( (Ilon+1) % Nlon + Nlon ) % Nlon;
                points_to_test.push_back( Index(0, 0, Ilat, Ineigh, 1, 1, Nlat, Nlon) );

                Ineigh = ( (Ilon-1) % Nlon + Nlon ) % Nlon;
                points_to_test.push_back( Index(0, 0, Ilat, Ineigh, 1, 1, Nlat, Nlon) );

                if (Ilat > 0) {
                    points_to_test.push_back( Index(0, 0, Ilat-1, Ilon, 1, 1, Nlat, Nlon) );
                }
                if (Ilat < Nlat-1) {
                    points_to_test.push_back( Index(0, 0, Ilat+1, Ilon, 1, 1, Nlat, Nlon) );
                }

                // So long as we still have points to test, keep testing!
                while ( points_to_test.size() > 0 ) {
                    Itest = points_to_test.front();
                    points_to_test.pop_front();
                    if ( pt_maps_to[Itest] < Npts ) { continue; } // already mapped, skip
                    else if (data.mask[Itest]) { continue; } // water, so skip
                    else {
                        // Otherwise, map it to Ipt, and add its neighbours to the test list
                        pt_maps_to[Itest] = Ipt;
                        num_mapped_points++;

                        // Get coordinates
                        Index1to4( Itest, Itime, Idepth, Ilat, Ilon, Ntime, Ndepth, Nlat, Nlon );

                        Ineigh = ( (Ilon+1) % Nlon + Nlon ) % Nlon;
                        Ineighbour = Index(0, 0, Ilat, Ineigh, 1, 1, Nlat, Nlon);
                        if ( pt_maps_to[Ineighbour] == Npts ) {
                            // Only add the neighbour if we haven't already looked at it
                            points_to_test.push_back( Ineighbour );
                        }

                        Ineigh = ( (Ilon-1) % Nlon + Nlon ) % Nlon;
                        Ineighbour = Index(0, 0, Ilat, Ineigh, 1, 1, Nlat, Nlon);
                        if ( pt_maps_to[Ineighbour] == Npts ) {
                            // Only add the neighbour if we haven't already looked at it
                            points_to_test.push_back( Ineighbour );
                        }

                        if (Ilat > 0) {
                            Ineighbour = Index(0, 0, Ilat-1, Ilon, 1, 1, Nlat, Nlon);
                            if ( pt_maps_to[Ineighbour] == Npts ) {
                                // Only add the neighbour if we haven't already looked at it
                                points_to_test.push_back( Ineighbour );
                            }
                        }
                        if (Ilat < Nlat-1) {
                            Ineighbour = Index(0, 0, Ilat+1, Ilon, 1, 1, Nlat, Nlon);
                            if ( pt_maps_to[Ineighbour] == Npts ) {
                                // Only add the neighbour if we haven't already looked at it
                                points_to_test.push_back( Ineighbour );
                            }
                        }
                    }
                }
            }
        }
        #if DEBUG >= 0
        fprintf( stdout, "Mapping %'zu land points onto %'zu 'islands' (contiguous land masses). %'zu points are coastal.\n",
                num_mapped_points+num_mapped_onto, num_mapped_onto, num_coastal );
        #endif
        Npts_mapped = Npts - num_mapped_points;
        Ncol = Npts - num_mapped_points - 1;
        Nrow = Npts - (num_mapped_points+num_mapped_onto) + num_coastal;
        #if DEBUG >= 2
        fprintf(stdout, "Total Number of Land Points: %'zu\n", num_land);
        fprintf(stdout, "Resulting Nrow and Ncol are %'zu and %'zu\n", Nrow, Ncol);
        #endif

        // Get the counts for how many rows and columns were removed
        size_t col_counter = 1, 
               row_counter = 0,
               land_counter = 0;
        for (Ipt = 1; Ipt < Npts; Ipt++) {
            if ( pt_maps_to[Ipt-1] != (Ipt-1) ) {
                col_counter++; 
            }

            if ( (pt_maps_to[Ipt-1] != (Ipt-1) ) or (Ipt-1 == 0) ) { 
                if ( all_land_neighbours[Ipt-1] == 1) { 
                    row_counter++; 
                }
            }

            if ( not(data.mask[Ipt-1]) ) { 
                land_counter++;
            }

            num_mapped_before_col[Ipt] = col_counter;
            num_mapped_before_row[Ipt] = row_counter;
            num_land_before[Ipt] = land_counter;
        }


    } else if (constants::GRID_TYPE == constants::GridType::LLC) {

        const size_t num_neighbours = data.num_neighbours;
        size_t neighbour_ind;

        for (Ipt = 0; Ipt < Npts; Ipt++) {
            bool is_land = not( data.mask[Ipt] );
            bool is_coastal = is_land and (all_land_neighbours[Ipt] == 0);

            // Water+coast points map to themselves
            if ( not(is_land) ) { 
                pt_maps_to[Ipt] = Ipt; 
            } else if ( not(collapse_coast) and is_coastal ) {
                pt_maps_to[Ipt] = Ipt; 
            } else {
                // Check if we've already mapped this point
                if (pt_maps_to[Ipt] == Npts) {
                    // Begin a depth-first search through the adjacency
                    // matrix, mapping all land-only neighbours to this one
                    // This point will become the representative for this 'island'
                    points_to_test.clear();
                    pt_maps_to[Ipt] = Ipt; 
                    num_mapped_onto++;
                    for (neighbour_ind = 0; neighbour_ind < num_neighbours; neighbour_ind++ ) {
                        points_to_test.push_back( data.adjacency_indices.at(Ipt)[neighbour_ind] );
                    }
                    // So long as we still have points to test, keep testing!
                    while ( points_to_test.size() > 0 ) {
                        Itest = points_to_test.front();
                        points_to_test.pop_front();

                        bool test_is_land = not( data.mask[Itest] );
                        bool test_is_coastal = test_is_land and (all_land_neighbours[Itest] == 0);
                            
                        if ( pt_maps_to[Itest] < Npts ) { continue; } // already mapped, skip
                        else if ( not(collapse_coast) and test_is_coastal ) { 
                            pt_maps_to[Itest] = Itest; // flag it as mapping to itself
                            continue; 
                        } 
                        else if ( not(test_is_land) ) {
                            pt_maps_to[Itest] = Itest; // flag it as mapping to itself
                            continue; 
                        } 
                        else {
                            // Otherwise, map it to Ipt, and add its neighbours to the test list
                            pt_maps_to[Itest] = Ipt;
                            num_mapped_points++;
                            if ( not(test_is_land) and ( collapse_coast or not(test_is_coastal) ) ) {
                                throw std::runtime_error( "Invalid point getting mapped" );
                            }
                            for (neighbour_ind = 0; neighbour_ind < num_neighbours; neighbour_ind++ ) {
                                Ineighbour = data.adjacency_indices.at(Itest)[neighbour_ind];
                                if ( pt_maps_to[Ineighbour] == Npts ) {
                                    // Only add the neighbour if we haven't already looked at it
                                    points_to_test.push_back( Ineighbour );
                                }
                            }
                        }
                    }
                }
            }
        }

        Npts_mapped = Npts - num_mapped_points;
        Ncol = Npts - (num_mapped_points + 1); // +1 accounts for the first point [south pole], removes all land except for a 
                                               // 'representative' of each continguous mass
        Nrow = Npts - (num_mapped_points+num_mapped_onto) + ( collapse_coast ? num_coastal : 0); // removes all land, but adds coastal back in
        #if DEBUG >= 0
        fprintf( stdout, "Mapping %'zu land points onto %'zu 'islands' (contiguous land masses). %'zu points are coastal.\n",
                num_mapped_points+num_mapped_onto, num_mapped_onto, num_coastal );
        fprintf(stdout, "Resulting Nrow and Ncol are %'zu and %'zu\n", Nrow, Ncol);
        #endif

        size_t col_counter = 1, 
               row_counter = 0,
               land_counter = 0,
               mapped_counter = 0,
               mapped_onto_counter = 0,
               coastal_counter = 0,
               mapped_coastal_counter = 0;
        for (Ipt = 1; Ipt < Npts; Ipt++) {

            bool prev_is_land = not(data.mask[Ipt-1]);
            bool prev_is_mapped = pt_maps_to[Ipt-1] != (Ipt-1);
            bool prev_is_coastal = prev_is_land and (all_land_neighbours[Ipt-1] == 0);
            bool prev_is_mapped_onto = prev_is_land and not(prev_is_mapped) and (not(prev_is_coastal) or collapse_coast);

            // Counting of land, coast, and islands before this one
            if ( prev_is_land ) { 
                land_counter++;
                if ( prev_is_coastal ) { coastal_counter++; }
                if ( prev_is_mapped ) { mapped_counter++; }
                if ( prev_is_mapped_onto ) { mapped_onto_counter++; }
                if ( prev_is_coastal and prev_is_mapped) { 
                    mapped_coastal_counter++; 
                    if ( not(collapse_coast) ) {
                        fprintf( stdout, " current number of mapped coastal points is %'zu\n", mapped_coastal_counter );
                    }
                }
            }

            //row_counter = (mapped_counter + island_counter) - coastal_counter;
            //row_counter = land_counter - coastal_counter;
            row_counter = mapped_counter + mapped_onto_counter;
            col_counter = mapped_counter + 1;

            // Some sanity checks / bounds checking
            bool curr_is_mapped = pt_maps_to[Ipt] != Ipt;
            if ( not(curr_is_mapped) ) {
                if ( row_counter > Ipt ) { throw std::range_error( "Row counter has exceeded point index." ); }
                if ( Ipt - row_counter >= Nrow ) { 
                    fprintf(stdout, "Ipt = %'zu, num_mapped = %'zu, mapped_counter = %'zu, mapped_onto_counter = %'zu\n", Ipt, num_mapped_points, mapped_counter, mapped_onto_counter );
                    throw std::range_error( "(Ipt - Row counter) exceeds Nrow." ); 
                }
                if ( col_counter > Ipt ) { throw std::range_error( "Column counter has exceeded point index." ); }
                if ( Ipt - col_counter >= Ncol ) { 
                    fprintf(stdout, "Ipt = %'zu, num_mapped = %'zu, mapped_counter = %'zu\n", Ipt, num_mapped_points, mapped_counter);
                    throw std::range_error( "(Ipt - Column counter) exceeds Ncol." ); 
                }
            }

            num_mapped_before_col[Ipt] = col_counter;
            num_mapped_before_row[Ipt] = row_counter;
            num_land_before[Ipt] = land_counter;

        }

        // If the last point is mapped, it doesn't impact anything, but still account for it
        // for the sake of enabling sanity checks
        bool last_is_land = not(data.mask[Npts-1]);
        bool last_is_mapped = pt_maps_to[Npts-1] != (Npts-1);
        bool last_is_coastal = last_is_land and (all_land_neighbours[Npts-1] == 0);
        bool last_is_mapped_onto = last_is_land and not(last_is_mapped) and (not(last_is_coastal) or collapse_coast);
        if ( last_is_mapped or last_is_mapped_onto ) { row_counter++; }
        if ( last_is_mapped ) { col_counter++; }
        

        // Sanity checks: make sure we accounted for all of the points
        if ( row_counter + Nrow != Npts ) {
            fprintf(stdout, "row_counter, Nrow, Npt = %'zu, %'zu, %'zu\n", row_counter, Nrow, Npts );
            throw std::runtime_error( "row_counter not accounting for all skipped row." );
        }
        if ( col_counter + Ncol != Npts ) {
            fprintf(stdout, "col_counter, Ncol, Npt = %'zu, %'zu, %'zu\n", col_counter, Ncol, Npts );
            throw std::runtime_error( "col_counter not accounting for all skipped columns." );
        }
    } else {
        // raise an error
        throw std::invalid_argument("Mesh grid option not recognized.");
    }

}
