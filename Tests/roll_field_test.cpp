#include <fenv.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <algorithm>
#include <math.h>
#include <vector>
#include <cassert>

#include "../functions.hpp"

void print_field( const std::vector<double> & field, const int & Nlat, const int & Nlon ) {
    for ( int Ilat = 0; Ilat < Nlat; ++Ilat ) {
        for ( int Ilon = 0; Ilon < Nlon; ++Ilon ) {
            fprintf(stdout, " %g ", field.at( Ilat * Nlon + Ilon ));
        }
        fprintf(stdout, "\n");
    }
}

void get_true_roll( std::vector<double> & arr, const int Iroll, const int & Nlat, const int & Nlon ) {
    double cntr = 0;
    int Ilon_adj;
    for ( int Ilat = 0; Ilat < Nlat; ++Ilat ) {
        for ( int Ilon = 0; Ilon < Nlon; ++Ilon ) {
            Ilon_adj = Ilon + Iroll >= Nlon ? Ilon + Iroll - Nlon : Ilon + Iroll;
            arr.at( Ilat * Nlon + Ilon_adj ) = cntr;
            cntr++;
        }
    }
}

int main(int argc, char *argv[]) {

    const int   Ntime = 1,
                Ndepth = 1,
                Nlat = 4,
                Nlon = 4;

    std::vector<double> field( Ntime*Ndepth*Nlat*Nlon, 0 ),
                        answer(Ntime*Ndepth*Nlat*Nlon, 0 );

    // Build starting array
    get_true_roll( field, 0, Nlat, Nlon );

    fprintf(stdout, "Original Field\n");
    print_field( field, Nlat, Nlon ); 
    fprintf(stdout, "\n\n");

    bool are_same;
    for (int Iroll = 1; Iroll < Nlon; ++Iroll) {
        roll_field( field, "lon", 1, Ntime, Ndepth, Nlat, Nlon );
        get_true_roll( answer, Iroll, Nlat, Nlon );

        are_same = true;
        for ( size_t ii = 0; ii < field.size(); ++ii) {
            are_same = are_same and ( std::fabs( field.at(ii) - answer.at(ii) ) < 1e-15 );
        }
        if (are_same) {
            fprintf( stdout, "PASSED FOR IROLL = %d\n", Iroll );
        } else {
            fprintf( stdout, "FAILED FOR IROLL = %d\n", Iroll );

            fprintf(stdout, " Iroll = %d\n", Iroll);
            fprintf(stdout, "\nComputed:\n");
            print_field( field, Nlat, Nlon ); 
            fprintf(stdout, "\nExpected:\n");
            print_field( answer, Nlat, Nlon ); 
            fprintf(stdout, "\n\n");
        }
    }

    return 0;

}
