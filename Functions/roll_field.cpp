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
 * @param[in]       dimension                   dimension to roll (currently must be "lon")
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

    for (int Itime = 0; Itime < Ntime; Itime++) {
        for (int Idepth = 0; Idepth < Ndepth; Idepth++) {
            for (int Ilat = 0; Ilat < Nlat; Ilat++) {

                // Since we're assuming a roll in lon (the last [i.e. fastest] index, just rotate in-place)
                index = Index(Itime, Idepth, Ilat, Nlon - 1,
                              Ntime, Ndepth, Nlat, Nlon     );

                rbegin  = field_to_roll.rbegin() + (Npts - index) - 1 ;
                rend    = rbegin + Nlon;

                std::rotate( rbegin, rbegin + roll_count, rend );

            }
        }
    }

}
