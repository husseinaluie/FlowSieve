#include <vector>
#include <string>
#include <assert.h>
#include "../../differentiation_tools.hpp"
#include "../../constants.hpp"
#include "../../functions.hpp"

void spher_derivative_transpose_at_point(
        const std::vector<double*> & deriv_vals,
        const std::vector<const std::vector<double>*> & fields,
        const std::vector<double> & grid,
        const std::string & dim,
        const dataset & source_data,
        const int Itime,
        const int Idepth,
        const int Ilat,
        const int Ilon,
        const std::vector<short int> & mask,
        const int order_of_deriv,
        const int diff_ord
        ) {

    // Confirm that input sizes match
    if (deriv_vals.size() != fields.size()) {
        throw std::runtime_error( "Input and output dimension mismatch in derivatives." );
    }
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

    if ( (not(do_dep)) and ( constants::GRID_TYPE != constants::GridType::MeshGrid ) ) {
        throw std::runtime_error("Transpose of derivatives only implemented for MeshGrids.");
    }

    if ( do_dep ) {
        throw std::runtime_error("Transpose of derivatives not implements for depth derivatives.");
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

    const bool periodic = do_dep ? false : 
                          do_lat ? constants::PERIODIC_Y : 
                          do_lon ? constants::PERIODIC_X : false;

    if ( do_lon and periodic ) {
        // if it's longitude and periodic, then we can just take
        // the negative of the regular derivative
        // since the stencils are anti-symmetric
        spher_derivative_at_point( deriv_vals, fields, grid, dim, source_data,
               Itime, Idepth, Ilat, Ilon, mask, order_of_deriv, diff_ord );
        for (int ii = 0; ii < num_deriv; ii++) {
            if (deriv_vals.at(ii) != NULL) {
                *(deriv_vals.at(ii)) = -(*(deriv_vals.at(ii)));
            }
        }
        return;
    }


    // get the transpose of the differentiation vector
    int index_in_stencil;
    const int Idim = periodic ? -1 : Iref;

    // ASSUMING UNIFORM GRID
    const double delta = grid.at(1) - grid.at(0);

    std::vector<double> ddl; 
    differentiation_vector_T(ddl, index_in_stencil, delta, Idim, Nref, diff_ord);

    for ( int II = 0; II < ddl.size(); II++ ) {
        int Ideriv = Iref - index_in_stencil + II;

        // Apply periodicity adjustment 
        //   (has no effect for non-periodic, since then Iderv >= 0 and Ideriv < Nref)
        Ideriv = apply_periodicity( Ideriv, Nref );

        index = source_data.local_index(    Itime, 
                                            do_dep ? Ideriv : Idepth, 
                                            do_lat ? Ideriv : Ilat, 
                                            do_lon ? Ideriv : Ilon
                                        );

        for (int ii = 0; ii < num_deriv; ii++) {
            if (deriv_vals.at(ii) != NULL) {
                #ifdef NDEBUG
                *(deriv_vals.at(ii)) += (*fields[ii])[index] * ddl[II];
                #else
                *(deriv_vals.at(ii)) += fields.at(ii)->at(index) * ddl.at(II);
                #endif
            }
        }
    }
}
