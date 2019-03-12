#include <vector>
#include <string>
#include "../../differentiation_tools.hpp"
#include "../../constants.hpp"
#include "../../functions.hpp"

double spher_derivative_at_point(
        const std::vector<double> & field,  /**< [in] field to differentiate */
        const std::vector<double> & grid,   /**< [in] (1D) latitude  or longitude array */
        const std::string & dim,            /**< [in] "lat" or "lon": dimensional along which to differentate */
        const int Itime,                    /**< [in] time index at which to differentiate */
        const int Idepth,                   /**< [in] depth index at which to differentiate */
        const int Ilat,                     /**< [in] latitude index at which to differentiate */
        const int Ilon,                     /**< [in] longitude at which to differentiate */
        const int Ntime,                    /**< [in] size of time dimension */
        const int Ndepth,                   /**< [in] size of depth dimension */
        const int Nlat,                     /**< [in] size of latitude dimension */
        const int Nlon,                     /**< [in] size of longitude dimension */
        const std::vector<double> & mask    /**< [in] (2D) array to distinguish land/water cells */
        ) {

    int index, Iref;
    double deriv_val = 0.;
    const int Nref = grid.size();

    if      (dim == "lon") { Iref = Ilon; }
    else if (dim == "lat") { Iref = Ilat; }
    else { 
        fprintf(stderr, "Illegal dimensions provided! %s given to %s\n", dim.c_str(), __FILE__);
        Iref = 0;
        index = 0;
    }
    int LB = Iref;
    int UB = Iref;

    // Differentiation vector
    std::vector<double> ddl(constants::DiffOrd + 1);

    // Assuming uniform grid
    const double dl = grid.at(1) - grid.at(0);

    // Build outwards to try and build the stencil, but stop when
    //   we either hit a land cell or have gone far enough.
    while (LB > 0) {
        if ( (Iref - LB) > constants::DiffOrd ) { break; }
       
        if      (dim == "lon") { index = Index(0, 0, Ilat, LB,   Ntime, Ndepth, Nlat, Nlon); }
        else if (dim == "lat") { index = Index(0, 0, LB,   Ilon, Ntime, Ndepth, Nlat, Nlon); }
        
        if (mask.at(index) == 0) { LB++; break; }

        LB--;
    }

    while (UB < Nref-1) {
        if ( (UB - Iref) > constants::DiffOrd ) { break; }
       
        if      (dim == "lon") { index = Index(0, 0, Ilat, UB,   Ntime, Ndepth, Nlat, Nlon); }
        else if (dim == "lat") { index = Index(0, 0, UB,   Ilon, Ntime, Ndepth, Nlat, Nlon); }

        if (mask.at(index) == 0) { UB--; break; }

        UB++;
    }

    // We've possibly made too large of a stencil, so now collapse it back down
    while (UB - LB > constants::DiffOrd) {
        if (UB - Iref > Iref - LB) { UB--; }
        else { LB++; }
    }

    if (UB - LB == constants::DiffOrd) {
        // If we have enough cells for differentiation, do it
        differentiation_vector(ddl, dl, Iref - LB);
        for (int IND = LB; IND <= UB; IND++) {

            if      (dim == "lon") { index = Index(Itime, Idepth, Ilat, IND,  Ntime, Ndepth, Nlat, Nlon); }
            else if (dim == "lat") { index = Index(Itime, Idepth, IND,  Ilon, Ntime, Ndepth, Nlat, Nlon); }

            deriv_val +=  field.at(index) * ddl.at(IND - LB);
        }
    } else {
        // Otherwise, return 0
        deriv_val = 0.;
    }

    return deriv_val;
}
