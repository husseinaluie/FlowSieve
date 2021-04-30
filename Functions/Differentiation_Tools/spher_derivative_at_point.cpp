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
        const int Itime,
        const int Idepth,
        const int Ilat,
        const int Ilon,
        const int Ntime,
        const int Ndepth,
        const int Nlat,
        const int Nlon,
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
    int index, Iref;
    const int Nref = grid.size();
    const bool do_lat = (dim == "lat");
    const bool do_lon = (dim == "lon");
    assert( do_lat ^ do_lon ); // xor

    // Verify input
    if      (do_lon) { Iref = Ilon; }
    else if (do_lat) { Iref = Ilat; }
    else { 
        fprintf(stderr, "Illegal dimension provided! %s given to %s\n", 
                dim.c_str(), __FILE__);
        assert(false);
    }

    // Determine lowest lower bound (LLB) and upperest upper bound (UUB)
    //   for the integration region. This essentially just depends on periodicity.
    const bool periodic = do_lat ? constants::PERIODIC_Y : constants::PERIODIC_X ;
    const int LLB = periodic ? Iref - Nref : 0 ;
    const int UUB = periodic ? Iref + Nref : Nref - 1 ;

    // Differentiation vector
    const int num_deriv_pts = diff_ord + order_of_deriv;
    std::vector<double> ddl(num_deriv_pts);

    // Assuming uniform grid
    const double dl = grid.at(1) - grid.at(0);

    // Build outwards to try and build the stencil, but stop when
    //   we either hit a land cell or have gone far enough.
    // lb (lower case) will be the periodicity-adjusted value of LB 
    int LB = Iref, lb;
    while (LB > LLB) {

        // If we have enough points, stop
        if ( (Iref - LB) > diff_ord ) { break; }
       
        // Check if the next point would be land
        lb = ( ( LB - 1 ) % Nref + Nref ) % Nref;
        if (do_lon) { index = Index(Itime, Idepth, Ilat, lb,   Ntime, Ndepth, Nlat, Nlon); }
        else        { index = Index(Itime, Idepth, lb,   Ilon, Ntime, Ndepth, Nlat, Nlon); }
        
        if ( mask.at(index) )   { LB--;  }
        else                    { break; }
    }

    // ub (lower case) will be the periodicity-adjusted value of UB 
    int UB = Iref, ub;
    while (UB < UUB) {

        // If we have enough points, stop
        if ( (UB - Iref) > diff_ord ) { break; }
       
        // Check if the next point would be land
        ub = ( ( UB + 1 ) % Nref + Nref ) % Nref;
        if (do_lon) { index = Index(Itime, Idepth, Ilat, ub,   Ntime, Ndepth, Nlat, Nlon); }
        else        { index = Index(Itime, Idepth, ub,   Ilon, Ntime, Ndepth, Nlat, Nlon); }

        if ( mask.at(index) )   { UB++;  }
        else                    { break; }
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
        if ( do_lon or (constants::UNIFORM_LAT_GRID)) {
            // Since we're on a uniform grid, we can use pre-computed
            //   differentiation coefficients
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
            ind = ( IND % Nref + Nref ) % Nref;

            if (do_lon) { index = Index(Itime, Idepth, Ilat, ind,  Ntime, Ndepth, Nlat, Nlon); }
            else        { index = Index(Itime, Idepth, ind,  Ilon, Ntime, Ndepth, Nlat, Nlon); }

            for (int ii = 0; ii < num_deriv; ii++) {
                if (deriv_vals.at(ii) != NULL) {
                    *(deriv_vals.at(ii)) += fields[ii]->at(index) * ddl.at(IND - LB);
                }
            }
        }
    } else if (diff_ord > 2) {
        // If we couldn't build a large enough stencil, then 
        //   try again with a lower order. This will allow us
        //   to fill in smaller areas with at least something,
        //   if not the most accurate something.
        spher_derivative_at_point(
                deriv_vals, fields, grid, dim,
                Itime, Idepth, Ilat, Ilon,
                Ntime, Ndepth, Nlat, Nlon,
                mask, order_of_deriv, diff_ord - 2);
    }
}
