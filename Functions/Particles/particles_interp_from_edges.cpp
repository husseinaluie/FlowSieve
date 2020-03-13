#include <math.h>
#include <algorithm>
#include <vector>
#include <omp.h>
#include "../../constants.hpp"
#include "../../functions.hpp"
#include "../../particles.hpp"

double particles_interp_from_edges(
        double ref_lat,
        double ref_lon,
        const std::vector<double> lat,
        const std::vector<double> lon,
        const std::vector<double> field,
        const std::vector<double> mask,
        const int left,
        const int right,
        const int bottom,
        const int top,
        const double time_p,
        const int Itime,
        const int Ntime
        ){

    const int Nlat = lat.size(),
              Nlon = lon.size();

    const double dlon = lon.at(1) - lon.at(0);
    const double dlat = lat.at(1) - lat.at(0);

    double lon_p, lat_p;

    // default values
    double  top_L_pre_val = 0.,
            top_R_pre_val = 0.,
            bot_L_pre_val = 0.,
            bot_R_pre_val = 0.,
            top_L_fut_val = 0.,
            top_R_fut_val = 0.,
            bot_L_fut_val = 0.,
            bot_R_fut_val = 0.,
            interp_val;

    // Get indices for the four corners at both times
    const int   BL_pre_ind = Index(Itime,   0, bottom, left,
                                   Ntime,   1, Nlat,   Nlon ),
                BR_pre_ind = Index(Itime,   0, bottom, right,
                                   Ntime,   1, Nlat,   Nlon ),
                TL_pre_ind = Index(Itime,   0, top,    left,
                                   Ntime,   1, Nlat,   Nlon ),
                TR_pre_ind = Index(Itime,   0, top,    right,
                                   Ntime,   1, Nlat,   Nlon ),
                BL_fut_ind = Index(Itime+1, 0, bottom, left,
                                   Ntime,   1, Nlat,   Nlon ),
                BR_fut_ind = Index(Itime+1, 0, bottom, right,
                                   Ntime,   1, Nlat,   Nlon ),
                TL_fut_ind = Index(Itime+1, 0, top,    left,
                                   Ntime,   1, Nlat,   Nlon ),
                TR_fut_ind = Index(Itime+1, 0, top,    right,
                                   Ntime,   1, Nlat,   Nlon );

    // Get field values at each corner
    if (top >= 0) {
        if ( mask.at(TL_pre_ind) == 1 ) { top_L_pre_val = field.at(TL_pre_ind); }
        if ( mask.at(TR_pre_ind) == 1 ) { top_R_pre_val = field.at(TR_pre_ind); }

        if ( mask.at(TL_fut_ind) == 1 ) { top_L_fut_val = field.at(TL_fut_ind); }
        if ( mask.at(TR_fut_ind) == 1 ) { top_R_fut_val = field.at(TR_fut_ind); }
    } 

    if (bottom >= 0) {
        if ( mask.at(BL_pre_ind) == 1 ) { bot_L_pre_val = field.at(BL_pre_ind); }
        if ( mask.at(BR_pre_ind) == 1 ) { bot_R_pre_val = field.at(BR_pre_ind); }

        if ( mask.at(BL_fut_ind) == 1 ) { bot_L_fut_val = field.at(BL_fut_ind); }
        if ( mask.at(BR_fut_ind) == 1 ) { bot_R_fut_val = field.at(BR_fut_ind); }
    }

    // Do the interpolation in time
    const double TL_I_val = (1. - time_p) * top_L_pre_val + time_p * top_L_fut_val,
                 TR_I_val = (1. - time_p) * top_R_pre_val + time_p * top_R_fut_val,
                 BL_I_val = (1. - time_p) * bot_L_pre_val + time_p * bot_L_fut_val,
                 BR_I_val = (1. - time_p) * bot_R_pre_val + time_p * bot_R_fut_val;


    // Interpolate in longitude
    if (ref_lon > lon.at(left)) {
        lon_p = ( ref_lon - lon.at(left) ) / dlon;
    } else {
        lon_p = 1 - ( lon.at(right) - ref_lon ) / dlon;
    }

    const double top_I_val = (1. - lon_p) * TL_I_val  +  lon_p * TR_I_val;
    const double bot_I_val = (1. - lon_p) * BL_I_val  +  lon_p * BR_I_val;

    // Interpolate in latitude
    //   For now, just say that things near the poles are 'broken'
    if ( (top < 0) or (bottom < 0) ) {
        interp_val = 0.;
    } else {
        lat_p = ( ref_lat - lat.at(bottom) ) / dlat;
        interp_val = (1. - lat_p) * bot_I_val  +  lat_p * top_I_val;
    }

    // Return the computed value
    return interp_val;

}
