#include <math.h>
#include <algorithm>
#include <vector>
#include <omp.h>
#include <cassert>
#include "../../constants.hpp"
#include "../../functions.hpp"
#include "../../particles.hpp"

double particles_interp_from_edges(
        double ref_lat,
        double ref_lon,
        const std::vector<double> & lat,
        const std::vector<double> & lon,
        const std::vector<double> * field,
        const std::vector<bool> & mask,
        const int left,
        const int right,
        const int bottom,
        const int top,
        const double time_p,
        const int Itime,
        const int Ntime
        ){

    const unsigned int  Nlat = lat.size(),
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
    const size_t    BL_pre_ind = Index(Itime,   0, bottom, left,
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

    const size_t cutoff = 1e19;
    if (BL_pre_ind > cutoff) {
        fprintf(stdout, "BL_pre_ind = %'zu: (Itime,l,r,b,t) = (%'d, %'d, %'d, %'d, %'d)\n",
                BL_pre_ind, Itime, left, right, bottom, top);
        assert(false);
    }
    if (BR_pre_ind > cutoff) {
        fprintf(stdout, "BR_pre_ind = %'zu: (Itime,l,r,b,t) = (%'d, %'d, %'d, %'d, %'d)\n",
                BR_pre_ind, Itime, left, right, bottom, top);
        assert(false);
    }
    if (TL_pre_ind > cutoff) {
        fprintf(stdout, "TL_pre_ind = %'zu: (Itime,l,r,b,t) = (%'d, %'d, %'d, %'d, %'d)\n",
                TL_pre_ind, Itime, left, right, bottom, top);
        assert(false);
    }
    if (BR_pre_ind > cutoff) {
        fprintf(stdout, "TR_pre_ind = %'zu: (Itime,l,r,b,t) = (%'d, %'d, %'d, %'d, %'d)\n",
                TR_pre_ind, Itime, left, right, bottom, top);
        assert(false);
    }

    if (BL_fut_ind > cutoff) {
        fprintf(stdout, "BL_fut_ind = %'zu: (Itime,l,r,b,t) = (%'d, %'d, %'d, %'d, %'d)\n",
                BL_fut_ind, Itime, left, right, bottom, top);
        assert(false);
    }
    if (BR_fut_ind > cutoff) {
        fprintf(stdout, "BR_fut_ind = %'zu: (Itime,l,r,b,t) = (%'d, %'d, %'d, %'d, %'d)\n",
                BR_fut_ind, Itime, left, right, bottom, top);
        assert(false);
    }
    if (TL_fut_ind > cutoff) {
        fprintf(stdout, "TL_fut_ind = %'zu: (Itime,l,r,b,t) = (%'d, %'d, %'d, %'d, %'d)\n",
                TL_fut_ind, Itime, left, right, bottom, top);
        assert(false);
    }
    if (BR_fut_ind > cutoff) {
        fprintf(stdout, "TR_fut_ind = %'zu: (Itime,l,r,b,t) = (%'d, %'d, %'d, %'d, %'d)\n",
                TR_fut_ind, Itime, left, right, bottom, top);
        assert(false);
    }

    // Get field values at each corner
    if (top >= 0) {
        if ( mask.at(TL_pre_ind) ) { top_L_pre_val = field->at(TL_pre_ind); }
        if ( mask.at(TR_pre_ind) ) { top_R_pre_val = field->at(TR_pre_ind); }

        if (Ntime > 1) {
            if ( mask.at(TL_fut_ind) ) { top_L_fut_val = field->at(TL_fut_ind); }
            if ( mask.at(TR_fut_ind) ) { top_R_fut_val = field->at(TR_fut_ind); }
        } else {
            top_L_fut_val = top_L_pre_val;
            top_R_fut_val = top_R_pre_val;
        }
    } 

    if (bottom >= 0) {
        if ( mask.at(BL_pre_ind) ) { bot_L_pre_val = field->at(BL_pre_ind); }
        if ( mask.at(BR_pre_ind) ) { bot_R_pre_val = field->at(BR_pre_ind); }

        if (Ntime > 1) {
            if ( mask.at(BL_fut_ind) ) { bot_L_fut_val = field->at(BL_fut_ind); }
            if ( mask.at(BR_fut_ind) ) { bot_R_fut_val = field->at(BR_fut_ind); }
        } else {
            bot_L_fut_val = bot_L_pre_val;
            bot_R_fut_val = bot_R_pre_val;
        }
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
