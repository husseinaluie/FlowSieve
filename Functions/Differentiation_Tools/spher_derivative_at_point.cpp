#include <vector>
#include <string>
#include <assert.h>
#include "../../differentiation_tools.hpp"
#include "../../constants.hpp"
#include "../../functions.hpp"

void spher_derivative_at_point(
        const std::vector<double*> & deriv_vals,
        const std::vector<const std::vector<double>*> & fields,
        const std::vector<double> & grid,
        const std::string & dim,
        const dataset & source_data,
        const int Itime,
        const int Idepth,
        const int Ilat,
        const int Ilon,
        const std::vector<bool> & mask,
        const int order_of_deriv,
        const int diff_ord
        ) {

    // Confirm that input sizes match
    assert(deriv_vals.size() == fields.size());
    const int num_deriv = deriv_vals.size();

    // Zero out before computing
    for (int ii = 0; ii < num_deriv; ii++) {
        if (deriv_vals.at(ii) != NULL) {
            *(deriv_vals.at(ii)) = 0.;
        }
    }

    // Check which derivative we're taking
    int index = 0;
    const bool do_dep = (dim == "depth");
    const bool do_lat = (dim == "lat");
    const bool do_lon = (dim == "lon");
    assert( do_dep ^ (do_lat ^ do_lon) ); // ^ = xor

    if ( (not(do_dep)) and ( constants::GRID_TYPE == constants::GridType::LLC ) ) {
        // If we're on an LLC grid, then we've already computed
        // the differentiation weights / stencils when we build
        // the adjacency matrix. So just apply those now and exit.

        // Unless we want a depth derivative, in which case just use the
        // previous method
        // NOTE:: THIS REALLY SHOULD BE TESTED TO MAKE SURE THERE ARE NO GREMLINS HIDING

        assert( Ilat == Ilon );
        index = (size_t) Ilat;

        double weight, val;
        size_t val_ind;

        for (int ii = 0; ii < num_deriv; ii++) {
            if (deriv_vals.at(ii) != NULL) {
                for ( size_t II = 0; II < source_data.num_neighbours+1; II++ ) {
                    if ( II < source_data.num_neighbours ) {
                        val_ind = source_data.adjacency_indices.at(index).at(II);
                    } else {
                        val_ind = index;
                    }
                    #if DEBUG >= 1
                    if (do_lat) {
                        weight = source_data.adjacency_ddlat_weights.at(Ilat).at(II);
                    } else {
                        weight = source_data.adjacency_ddlon_weights.at(Ilat).at(II);
                    }
                    val = fields[ii]->at( val_ind );
                    #else
                    if (do_lat) {
                        weight = source_data.adjacency_ddlat_weights[Ilat][II];
                    } else {
                        weight = source_data.adjacency_ddlon_weights[Ilat][II];
                    }
                    val = (*fields[ii])[ val_ind ];
                    #endif
                    *(deriv_vals.at(ii)) += val * weight;
                }
            }
        }

        return;
    }

    const int   Ntime   = source_data.Ntime,    // this is the MPI-local Ntime, not the full Ntime
                Ndepth  = source_data.Ndepth,   // this is the MPI-local Ndepth, not the full Ndepth
                Nlat    = source_data.Nlat,
                Nlon    = source_data.Nlon;

    int Iref = do_dep ? Idepth :
               do_lat ? Ilat :
               do_lon ? Ilon : -1;
    const int Nref = grid.size();

    // If it's a singleton dimension, just return zeros (we zeroed out earlier)
    //      this if for the case of no actual depth values
    if (Nref == 1) { return; }

    // Determine lowest lower bound (LLB) and upperest upper bound (UUB)
    //   for the integration region. This essentially just depends on periodicity.
    const bool periodic = do_dep ? false : 
                          do_lat ? constants::PERIODIC_Y : 
                          do_lon ? constants::PERIODIC_X : false;
    const int LLB = periodic ? Iref - Nref : 0 ;
    const int UUB = periodic ? Iref + Nref : Nref - 1 ;

    // Differentiation vector
    const int num_deriv_pts = diff_ord + order_of_deriv;
    //std::vector<double> ddl(num_deriv_pts);
    std::vector<double> ddl(0);

    // Build outwards to try and build the stencil, but stop when
    //   we either hit a land cell or have gone far enough.
    // lb (lower case) will be the periodicity-adjusted value of LB 
    int LB = Iref, lb;
    while (LB > LLB) {

        // If we have enough points, stop
        if ( (Iref - LB) >= num_deriv_pts ) { break; }
       
        // Check if the next point would be land
        lb = ( ( LB - 1 ) % Nref + Nref ) % Nref;
        index = Index( Itime, do_dep ? lb : Idepth, do_lat ? lb : Ilat, do_lon ? lb : Ilon,
                       Ntime, Ndepth,               Nlat,               Nlon );
        

        #if DEBUG >= 1
        if ( mask.at(index) )   { LB--;  }  // If next point is still water, extend stencil over it
        #else
        if ( mask[index] )      { LB--;  }  // If next point is still water, extend stencil over it
        #endif
        else                    { break; }  // Otherwise, halt [ without extending stencil ]
    }

    // ub (lower case) will be the periodicity-adjusted value of UB 
    int UB = Iref, ub;
    while (UB < UUB) {

        // If we have enough points, stop
        if ( (UB - Iref) >= num_deriv_pts ) { break; }
       
        // Check if the next point would be land
        ub = ( ( UB + 1 ) % Nref + Nref ) % Nref;
        index = Index( Itime, do_dep ? ub : Idepth, do_lat ? ub : Ilat, do_lon ? ub : Ilon,
                       Ntime, Ndepth,               Nlat,               Nlon );

        #if DEBUG >= 1
        if ( mask.at(index) )   { UB++;  }  // If next point is still water, extend stencil over it
        #else
        if ( mask[index] )      { UB++;  }  // If next point is still water, extend stencil over it
        #endif
        else                    { break; }  // Otherwise, halt [ without extending stencil ]
    }

    // NOTE
    //   In the case of periodicity, LB may be negative and UB may exceed Nref
    //     this means that subtraction still gives the correct number of points
    //     in the stencil, but that a periodicity-adjusted value will be needed
    //     when determining logical indices.

    // We've possibly made too large of a stencil, so now collapse it back down
    while (UB - LB + 1 > num_deriv_pts) {
        if ((UB - Iref > Iref - LB) and (UB >= Iref)) { UB--; }
        else { LB++; }
    }

    // We're including LB and UB in our stencil, so the stencil
    //   has UB - LB + 1 points. The requisit number of points is
    //   num_deriv_pts.
    // Again, lower case (ind) will be the periodicty-adjusted value
    int ind; 
    if (UB - LB + 1 == num_deriv_pts) {
        // If we have enough cells for differentiation, do it
        if ( do_lon or (do_lat and constants::UNIFORM_LAT_GRID)) {
            // Since we're on a uniform grid, we can use pre-computed
            //   differentiation coefficients
            const double dl = grid.at(1) - grid.at(0);
            differentiation_vector(ddl, dl, Iref - LB, order_of_deriv, diff_ord);
        } else {
            // We're on a non-uniform grid, so we can't guarantee the
            //   differentiation coefficients a priori, so we need
            //   to actually compute them now.
            // This will get expensive (or ugly...) for higher orders of accuracy.
            // NOTE: This CANNOT handle periodicity
            non_uniform_diff_vector(ddl, grid, Iref, LB, UB, diff_ord);
        }

        for (int IND = LB; IND <= UB; IND++) {

            // Apply periodicity adjustment 
            //   (has no effect for non-periodic, since then LB >= 0 and UB < Nref)
            ind = ( IND % Nref + Nref ) % Nref;

            index = Index( Itime, do_dep ? ind : Idepth, do_lat ? ind : Ilat, do_lon ? ind : Ilon,
                           Ntime, Ndepth,                Nlat,                Nlon );

            for (int ii = 0; ii < num_deriv; ii++) {
                if (deriv_vals.at(ii) != NULL) {
                    #if DEBUG >= 1
                    *(deriv_vals.at(ii)) += fields[ii]->at(index) * ddl.at(IND - LB);
                    #else
                    *(deriv_vals.at(ii)) += (*fields[ii])[index] * ddl[IND - LB];
                    #endif
                }
            }
        }
    } else if (diff_ord > 2) {
        // If we couldn't build a large enough stencil, then 
        //   try again with a lower order. This will allow us
        //   to fill in smaller areas with at least something,
        //   if not the most accurate something.
        spher_derivative_at_point(
                deriv_vals, fields, grid, dim, source_data,
                Itime, Idepth, Ilat, Ilon,
                mask, order_of_deriv, diff_ord - 2);
    }
}
