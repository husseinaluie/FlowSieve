#include <math.h>
#include <algorithm>
#include <vector>
#include <cassert>
#include "../functions.hpp"

/*!
 *
 * \brief Roll field along dimension
 *
 * Currently hard-coded to roll along lon dimension only.
 *
 * @param[in,out]   field_to_roll               field to be rolled
 * @param[in]       dimension                   dimension to roll (currently must be 0)
 * @param[in]       roll_count                  ammount by which to roll (currently must be 1)
 * @param[in]       Ntime,Ndepth,Nlat,Nlon      (MPI-local) dimension sizes
 *
 */
void roll_field(
        std::vector<double> & field_to_roll,
        const std::string dimension,
        const int roll_count,
        const int Ntime,
        const int Ndepth,
        const int Nlat,
        const int Nlon
        ) {

    assert(roll_count == 1);

    // If roll_count is zero, just stop now
    if (roll_count == 0) {return;}

    size_t index;
    const size_t Npts = ( (size_t) Ntime) * ( (size_t) Ndepth) * ( (size_t) Nlat) * ( (size_t) Nlon);
    std::vector<double> work_arr( Nlon );

    std::vector<double>::reverse_iterator rbegin, rend;
    //std::reverse_iterator<std::string::iterator>

    for (int Itime = 0; Itime < Ntime; Itime++) {
        for (int Idepth = 0; Idepth < Ndepth; Idepth++) {
            for (int Ilat = 0; Ilat < Nlat; Ilat++) {

                // Since we're assuming a roll in lon (the last [i.e. fastest] index, just rotate in-place)
                index = Index(Itime, Idepth, Ilat, Nlon - 1,
                              Ntime, Ndepth, Nlat, Nlon     );

                rbegin  = field_to_roll.rbegin() + (Npts - index) - 1 ;
                rend    = rbegin + Nlon;

                std::rotate( rbegin, rbegin + roll_count, rend );

                /*
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
                */

            }
        }
    }

}
