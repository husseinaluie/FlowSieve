
#include <algorithm>
#include <math.h>
#include "../functions.hpp"
#include <vector>
#include "../preprocess.hpp"


void get_coast(
        dataset & source_data
        ) {

    std::vector<short int> &mask = source_data.mask;
    std::vector<short int> &coast = source_data.coast;

    const size_t Npts = mask.size();

    coast.assign( mask.size(), -1 );
    source_data.coastal_indices.clear();
    source_data.inner_coastal_indices.clear();

    size_t Ipt, Jpt;

    for (Ipt = 0; Ipt < Npts; Ipt++) {
        if ( mask[Ipt] ) { coast[Ipt] = 0; }
    }


    const int Nlat = source_data.Nlat;
    const int Nlon = source_data.Nlon;

    int Itime, Idepth, Ilat, Ilon, Ineigh;


    for ( int coast_level = 1; coast_level <= constants::DiffOrd; coast_level++ ) {
        int num_coast = 0;
        for (Ipt = 0; Ipt < Npts; Ipt++) {

            // if water or other coast, skip it
            if ( coast[Ipt] >= 0 ) { continue; }

            bool is_coast = false;

            source_data.index1to4_local( Ipt, Itime, Idepth, Ilat, Ilon );

            // check each neighbour
            Ineigh = apply_periodicity( Ilon+1, Nlon );
            Jpt = source_data.local_index(0, 0, Ilat, Ineigh);
            if ( coast[Jpt] == coast_level - 1 ) { is_coast = true; }

            Ineigh = apply_periodicity( Ilon-1, Nlon );
            Jpt = source_data.local_index(0, 0, Ilat, Ineigh);
            if ( coast[Jpt] == coast_level - 1 ) { is_coast = true; }

            if (Ilat > 0) {
                Jpt = source_data.local_index(0, 0, Ilat-1, Ilon);
                if ( coast[Jpt] == coast_level - 1 ) { is_coast = true; }
            }

            if (Ilat < Nlat-1) {
                Jpt = source_data.local_index(0, 0, Ilat+1, Ilon);
                if ( coast[Jpt] == coast_level - 1 ) { is_coast = true; }
            }

            if (is_coast) { 
                num_coast++;
                coast[Ipt] = coast_level;
            }
        }
        fprintf( stdout, "Identified %d points at coastal depth %d.\n", num_coast, coast_level );
    }

    /*
    // First, we do a large loop to simply
    //   count how many coast points we have.
    int num_coast = 0;
    for (Ipt = 0; Ipt < Npts; Ipt++) {
        
        if ( mask[Ipt] ) { continue; }

        bool is_coast = false;

        source_data.index1to4_local( Ipt, Itime, Idepth, Ilat, Ilon );

        // check each neighbour
        Ineigh = apply_periodicity( Ilon+1, Nlon );
        Jpt = source_data.local_index(0, 0, Ilat, Ineigh);
        if ( mask[Jpt] ) { is_coast = true; }

        Ineigh = apply_periodicity( Ilon-1, Nlon );
        Jpt = source_data.local_index(0, 0, Ilat, Ineigh);
        if ( mask[Jpt] ) { is_coast = true; }

        if (Ilat > 0) {
            Jpt = source_data.local_index(0, 0, Ilat-1, Ilon);
            if ( mask[Jpt] ) { is_coast = true; }
        }

        if (Ilat < Nlat-1) {
            Jpt = source_data.local_index(0, 0, Ilat+1, Ilon);
            if ( mask[Jpt] ) { is_coast = true; }
        }

        if (is_coast) { 
            num_coast++;
            source_data.coastal_indices.push_back(Ipt);
            coast[Ipt] = true;
        }
    }

    #if DEBUG>0
    fprintf(stdout, "      ... found %d coast points (%.4g percent)\n", num_coast, (100.*num_coast) / Npts);
    #endif

    if ( num_coast != source_data.coastal_indices.size() ) {
        throw std::runtime_error("Mismatch in size of coast");
    }

    //
    //// And now we do it all again to get the next ring of coastal points
    //
    num_coast = 0;
    for (Ipt = 0; Ipt < Npts; Ipt++) {
        
        if ( mask[Ipt] or coast[Ipt] ) { continue; }

        bool is_coast = false;

        source_data.index1to4_local( Ipt, Itime, Idepth, Ilat, Ilon );

        // check each neighbour
        Ineigh = apply_periodicity( Ilon+1, Nlon );
        Jpt = source_data.local_index(0, 0, Ilat, Ineigh);
        if ( coast[Jpt] ) { is_coast = true; }

        Ineigh = apply_periodicity( Ilon-1, Nlon );
        Jpt = source_data.local_index(0, 0, Ilat, Ineigh);
        if ( coast[Jpt] ) { is_coast = true; }

        if (Ilat > 0) {
            Jpt = source_data.local_index(0, 0, Ilat-1, Ilon);
            if ( coast[Jpt] ) { is_coast = true; }
        }

        if (Ilat < Nlat-1) {
            Jpt = source_data.local_index(0, 0, Ilat+1, Ilon);
            if ( coast[Jpt] ) { is_coast = true; }
        }

        if (is_coast) { 
            num_coast++;
            source_data.inner_coastal_indices.push_back(Ipt);
        }
    }

    #if DEBUG>0
    fprintf(stdout, "      ... found %d inner coast points (%.4g percent)\n", num_coast, (100.*num_coast) / Npts);
    #endif

    if ( num_coast != source_data.inner_coastal_indices.size() ) {
        throw std::runtime_error("Mismatch in size of inner coast");
    }
    */

}
