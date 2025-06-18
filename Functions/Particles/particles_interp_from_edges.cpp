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
        const int top
        ){

    const unsigned int  Nlat = lat.size(),
                        Nlon = lon.size();

    const double dlon = lon.at(1) - lon.at(0);
    const double dlat = lat.at(1) - lat.at(0);

    double lon_p, lat_p, interp_val;

    // Get indices for the four corners at both times
    const size_t    BL_ind = Index(0,   0, bottom, left,
                                   1,   1, Nlat,   Nlon ),
                    BR_ind = Index(0,   0, bottom, right,
                                   1,   1, Nlat,   Nlon ),
                    TL_ind = Index(0,   0, top,    left,
                                   1,   1, Nlat,   Nlon ),
                    TR_ind = Index(0,   0, top,    right,
                                   1,   1, Nlat,   Nlon );

    // Verify that all of the indices are valid
    const size_t cutoff = field->size();
    if (BL_ind >= cutoff) {
        fprintf(stderr, "BL_ind = %'zu: (l,r,b,t) = (%'d, %'d, %'d, %'d)\n",
                BL_ind, left, right, bottom, top);
        throw std::runtime_error("Unable to find bounding box for particle on grid.");
    }
    if (BR_ind >= cutoff) {
        fprintf(stderr, "BR_ind = %'zu: (l,r,b,t) = (%'d, %'d, %'d, %'d)\n",
                BR_ind, left, right, bottom, top);
        throw std::runtime_error("Unable to find bounding box for particle on grid.");
    }
    if (TL_ind >= cutoff) {
        fprintf(stderr, "TL_ind = %'zu: (l,r,b,t) = (%'d, %'d, %'d, %'d)\n",
                TL_ind, left, right, bottom, top);
        throw std::runtime_error("Unable to find bounding box for particle on grid.");
    }
    if (TR_ind >= cutoff) {
        fprintf(stderr, "TR_ind = %'zu: (l,r,b,t) = (%'d, %'d, %'d, %'d)\n",
                TR_ind, left, right, bottom, top);
        throw std::runtime_error("Unable to find bounding box for particle on grid.");
    }

    // Make sure the bounding box is correct
    double left_lon = lon.at(left), 
           right_lon = lon.at(right),
           bottom_lat = lat.at(bottom),
           top_lat = lat.at(top);
    if (not( (right < left) or ( ( left_lon <= ref_lon ) and ( ref_lon <= right_lon ) ) ) ) {
        throw std::runtime_error("Bounding box for particle on grid invalid.");
    }
    if (not( (top <= bottom) or ( ( bottom_lat <= ref_lat ) and ( ref_lat <= top_lat ) ) ) ) {
        throw std::runtime_error("Bounding box for particle on grid invalid.");
    }

    // Interpolate in longitude
    double lon_diff = std::fabs(ref_lon - left_lon);
    if ( lon_diff > M_PI ) { lon_diff = 2 * M_PI - lon_diff; }
    lon_p = lon_diff / dlon;
    if (not( (-1e-3 <= lon_p) and (lon_p <= 1 + 1e-3) ) ) {
        throw std::runtime_error("Invalid interpolant.");
    }

    // Get field values at each corner
    bool top_mask = true;
    double top_I_val = 0, TL_I_val = 0, TR_I_val = 0;
    if (top >= 0) {
        if ( mask.at(TL_ind) ) { TL_I_val = field->at(TL_ind); }
        if ( mask.at(TR_ind) ) { TR_I_val = field->at(TR_ind); }

        if ( mask.at(TL_ind) and mask.at(TR_ind) ) {
            top_I_val = (1. - lon_p) * TL_I_val  +  lon_p * TR_I_val;
        } else if ( mask.at(TL_ind) ) {
            top_I_val = TL_I_val;
        } else if ( mask.at(TR_ind) ) {
            top_I_val = TR_I_val;
        } else {
            top_mask = false;
            top_I_val = constants::fill_value;
        }
    } 

    bool bot_mask = true;
    double bot_I_val = 0, BL_I_val = 0, BR_I_val = 0;
    if (bottom >= 0) {
        if ( mask.at(BL_ind) ) { BL_I_val = field->at(BL_ind); }
        if ( mask.at(BR_ind) ) { BR_I_val = field->at(BR_ind); }

        if ( mask.at(BL_ind) and mask.at(BR_ind) ) {
            bot_I_val = (1. - lon_p) * BL_I_val  +  lon_p * BR_I_val;
        } else if ( mask.at(BL_ind) ) {
            bot_I_val = BL_I_val;
        } else if ( mask.at(BR_ind) ) {
            bot_I_val = BR_I_val;
        } else {
            bot_mask = false;
            bot_I_val = constants::fill_value;
        }
    }

    // Interpolate in latitude
    //   For now, just say that things near the poles are 'broken'
    if ( top == bottom ) {
        interp_val = top_I_val; // top = bottom, so they're the same
                                // this is the case where we're outside
                                // of the lat bounds
    } else {
        if ( top_mask and bot_mask ) {
            lat_p = ( ref_lat - lat.at(bottom) ) / dlat;
            assert( (-1e-3 <= lat_p) and (lat_p <= 1 + 1e-3) );
            interp_val = (1. - lat_p) * bot_I_val  +  lat_p * top_I_val;
        } else if ( top_mask ) {
            interp_val = top_I_val;
        } else if ( bot_mask ) {
            interp_val = bot_I_val;
        } else {
            interp_val = constants::fill_value;
        }
    }

    // Return the computed value
    return interp_val;

}
