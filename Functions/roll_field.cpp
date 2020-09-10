#include <math.h>
#include <algorithm>
#include <vector>
#include "../functions.hpp"

void roll_field(
        std::vector<double> & field_to_roll,
        const std::string dimension,
        const int roll_count,
        const int Ntime,
        const int Ndepth,
        const int Nlat,
        const int Nlon
        ) {

    // If roll_count is zero, just stop now
    if (roll_count == 0) {return;}

    size_t index;
    std::vector<double> work_arr( Nlon );

    for (int Itime = 0; Itime < Ntime; Itime++) {
        for (int Idepth = 0; Idepth < Ndepth; Idepth++) {
            for (int Ilat = 0; Ilat < Nlat; Ilat++) {

                // First, pull out the Ilon slice and store it in the work array
                for (int Ilon = 0; Ilon < Nlon; Ilon++) {
                    index = Index(Itime, Idepth, Ilat, Ilon,
                                  Ntime, Ndepth, Nlat, Nlon);
                    work_arr.at(Ilon) = field_to_roll.at(index);
                }

                // Now rotate the work array
                if (roll_count > 0) {
                    std::rotate( work_arr.rbegin(), work_arr.rbegin() + roll_count, work_arr.rend());
                } else {
                    std::rotate( work_arr.begin(),  work_arr.begin() - roll_count,  work_arr.end());
                }

                // Now copy the rotate vector back into the main array
                for (int Ilon = 0; Ilon < Nlon; Ilon++) {
                    index = Index(Itime, Idepth, Ilat, Ilon,
                                  Ntime, Ndepth, Nlat, Nlon);
                    field_to_roll.at(index) = work_arr.at(Ilon);
                }

            }
        }
    }

}
