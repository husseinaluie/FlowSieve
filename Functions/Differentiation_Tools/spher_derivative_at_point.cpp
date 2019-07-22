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
        const std::vector<double> & mask,
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

    int index, Iref;
    const int Nref = grid.size();
    const bool do_lat = (dim == "lat");
    const bool do_lon = (dim == "lon");
    assert( do_lat ^ do_lon ); // xor

    if      (do_lon) { Iref = Ilon; }
    else if (do_lat) { Iref = Ilat; }
    else { 
        fprintf(stderr, "Illegal dimension provided! %s given to %s\n", 
                dim.c_str(), __FILE__);
        assert(false);
    }


    // Differentiation vector
    std::vector<double> ddl(diff_ord + 1);

    // Assuming uniform grid
    const double dl = grid.at(1) - grid.at(0);

    // Build outwards to try and build the stencil, but stop when
    //   we either hit a land cell or have gone far enough.
    int LB = Iref;
    while (LB > 0) {
        if ( (Iref - LB) > diff_ord ) { break; }
       
        if (do_lon) { index = Index(0, 0, Ilat, LB,   Ntime, Ndepth, Nlat, Nlon); }
        else        { index = Index(0, 0, LB,   Ilon, Ntime, Ndepth, Nlat, Nlon); }
        
        if (mask.at(index) == 0) { LB++; break; }

        LB--;
    }

    int UB = Iref;
    while (UB < Nref-1) {
        if ( (UB - Iref) > diff_ord ) { break; }
       
        if (do_lon) { index = Index(0, 0, Ilat, UB,   Ntime, Ndepth, Nlat, Nlon); }
        else        { index = Index(0, 0, UB,   Ilon, Ntime, Ndepth, Nlat, Nlon); }

        if (mask.at(index) == 0) { UB--; break; }

        UB++;
    }

    // We've possibly made too large of a stencil, so now collapse it back down
    while (UB - LB > diff_ord) {
        if ((UB - Iref > Iref - LB) and (UB >= Iref)) { UB--; }
        else { LB++; }
    }

    // We're including LB and UB in our stencil, so the stencil
    //   has UB - LB + 1 points. The requisit number of points is
    //   diff_ord + 1.
    if (UB - LB + 1 == diff_ord + 1) {
        // If we have enough cells for differentiation, do it
        differentiation_vector(ddl, dl, Iref - LB, diff_ord);
        for (int IND = LB; IND <= UB; IND++) {

            if (do_lon) { index = Index(Itime, Idepth, Ilat, IND,  
                                        Ntime, Ndepth, Nlat, Nlon); }
            else        { index = Index(Itime, Idepth, IND,  Ilon, 
                                        Ntime, Ndepth, Nlat, Nlon); }

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
                mask, diff_ord - 2);
    }
}
