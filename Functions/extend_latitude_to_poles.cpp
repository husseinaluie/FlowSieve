#include <math.h>
#include <algorithm>
#include <vector>
#include <cassert>
#include "../functions.hpp"
#include "../constants.hpp"
#include <mpi.h>


void extend_latitude_to_poles(
        const std::vector<double> & original_latitude,
        std::vector<double> & extended_latitude,
        int & orig_Ilat_start,
        const bool IS_DEGREES,
        const MPI_Comm comm
        ) {

    const int Nlat_original = original_latitude.size();

    // Assuming a constant spacing for the extension
    const double dlat = original_latitude.at(1) - original_latitude.at(0);

    const double pole_lat = IS_DEGREES ? 90. : (M_PI/2);

    // Determine how many points needed to pad to southern pole
    int pad_south = floor( ( original_latitude.at(0) - (-pole_lat) ) / dlat );
    // If our padding would put a point on the actual pole, pull back by one.
    if ( (pad_south > 0) and ( original_latitude.at(0) - pad_south * dlat <= -0.999*pole_lat ) ) { pad_south--; }
    orig_Ilat_start = pad_south;

    // Determine how many points needed to pad to northern pole
    int pad_north = floor( ( pole_lat - original_latitude.back() ) / dlat );
    // If our padding would put a point on the actual pole, pull back by one.
    if ( (pad_north > 0) and ( original_latitude.back() + pad_north * dlat >= 0.999*pole_lat ) ) { pad_north--; }

    #if DEBUG>=0
    int wRank=-1;
    MPI_Comm_rank( comm, &wRank );
    if (wRank == 0) {
        if ( (pad_south > 0) or (pad_north > 0) ) {
            fprintf( stdout, "NOTE :: Extending latitude grid to the poles. "); 
            if (pad_south > 0) { fprintf( stdout, "%'d points added to south ", pad_south ); }
            if (pad_north > 0) { fprintf( stdout, "%'d points added to north ", pad_north ); }
            fprintf( stdout, "\n");
        }
    }
    #endif

    // Create new latitude grid
    const int Nlat_extended = Nlat_original + pad_north + pad_south;
    extended_latitude.resize( Nlat_extended );

    // Copy in the original values
    //fprintf( stdout, " Copying in the original values\n " );
    for ( int Ilat = 0; Ilat < Nlat_original; Ilat++ ) {
        extended_latitude.at( Ilat + pad_south ) = original_latitude.at( Ilat );
    }

    // Fill in the southern extension
    //fprintf( stdout, " Filling in the southern extension\n " );
    for ( int Ipad = 0; Ipad < pad_south; Ipad++ ) {
        extended_latitude.at( pad_south - (Ipad+1) ) = original_latitude.at( 0 ) - dlat * (Ipad+1);
    }

    // Fill in the northern extension
    //fprintf( stdout, " Filling in the northern extension\n " );
    for ( int Ipad = 0; Ipad < pad_north; Ipad++ ) {
        extended_latitude.at( pad_south + Nlat_original + Ipad ) = original_latitude.back() + dlat * (Ipad+1);
    }

}
