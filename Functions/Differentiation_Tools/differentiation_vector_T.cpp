#include <vector>
#include <math.h>
#include <cassert>
#include "../../differentiation_tools.hpp"
#include "../../constants.hpp"

// transpose of differentiation vector
void differentiation_vector_T(
        std::vector<double> & diff_array,
        int & index_in_stencil,
        const double delta,
        const int Idim,
        const int Ndim,
        const int diff_ord
        ) {

    diff_array.clear();
    diff_array.reserve(32);

    std::vector<double>::iterator i0;
    i0 = diff_array.begin();

    double scale_factor = 1.;

    index_in_stencil = -1;

    if (diff_ord == 2) {

        scale_factor = 1.;
        // Use second order
        if ( Idim == 0 )            { diff_array.insert(i0, { -1.5, -0.5             }); }
        else if ( Idim == 1 )       { diff_array.insert(i0, {  2,    0,   -0.5       }); }
        else if ( Idim == 2 )       { diff_array.insert(i0, { -0.5,  0.5,  0,   -0.5 }); }
        else if ( Idim == Ndim-3 )  { diff_array.insert(i0, {  0.5,  0,   -0.5,  0.5 }); }
        else if ( Idim == Ndim-2 )  { diff_array.insert(i0, {  0.5,  0,   -2         }); }
        else if ( Idim == Ndim-1 )  { diff_array.insert(i0, {  0.5,  1.5             }); }
        else                        { diff_array.insert(i0, {  0.5,  0,   -0.5       }); }

        // now say where the target point is in the stencil
        if ( Idim == 0 )            { index_in_stencil = 0; }
        else if ( Idim == 1 )       { index_in_stencil = 1; }
        else if ( Idim == 2 )       { index_in_stencil = 2; }
        else if ( Idim == Ndim-3 )  { index_in_stencil = 1; }
        else if ( Idim == Ndim-2 )  { index_in_stencil = 1; }
        else if ( Idim == Ndim-1 )  { index_in_stencil = 1; }
        else                        { index_in_stencil = 1; }
    } else if (diff_ord == 4) {

        scale_factor = 12;
        // Use fourth order
        if ( Idim == 0 )            { diff_array.insert(i0, { -25,  -3,  1                  }); }
        else if ( Idim == 1 )       { diff_array.insert(i0, {  48, -10, -8,   1             }); }
        else if ( Idim == 2 )       { diff_array.insert(i0, { -36,  18,  0,  -8,  1         }); }
        else if ( Idim == 3 )       { diff_array.insert(i0, {  16,  -6,  8,   0, -8,   1    }); }
        else if ( Idim == 4 )       { diff_array.insert(i0, {  -3,   1, -1,   8,  0,  -8, 1 }); }
        else if ( Idim == Ndim-5 )  { diff_array.insert(i0, {  -1,   8,  0,  -8,  1,  -1, 3 }); }
        else if ( Idim == Ndim-4 )  { diff_array.insert(i0, {  -1,   8,  0,  -8,  6, -16    }); }
        else if ( Idim == Ndim-3 )  { diff_array.insert(i0, {  -1,   8,  0, -18, 36         }); }
        else if ( Idim == Ndim-2 )  { diff_array.insert(i0, {  -1,   8, 10, -48             }); }
        else if ( Idim == Ndim-1 )  { diff_array.insert(i0, {  -1,   3, 25                  }); }
        else                        { diff_array.insert(i0, {  -1,   8,  0, -8, 1           }); }

        // now say where the target point is in the stencil
        if ( Idim == 0 )            { index_in_stencil = 0; }
        else if ( Idim == 1 )       { index_in_stencil = 1; }
        else if ( Idim == 2 )       { index_in_stencil = 2; }
        else if ( Idim == 3 )       { index_in_stencil = 3; }
        else if ( Idim == 4 )       { index_in_stencil = 4; }
        else if ( Idim == Ndim-5 )  { index_in_stencil = 2; }
        else if ( Idim == Ndim-4 )  { index_in_stencil = 2; }
        else if ( Idim == Ndim-3 )  { index_in_stencil = 2; }
        else if ( Idim == Ndim-2 )  { index_in_stencil = 2; }
        else if ( Idim == Ndim-1 )  { index_in_stencil = 2; }
        else                        { index_in_stencil = 2; }
    }

    #ifndef NDEBUG
    if ( index_in_stencil == -1 ) {
        throw std::runtime_error("Derivative stencil index not set.");
    }

    if ( index_in_stencil >= diff_array.size() ) {
        throw std::runtime_error("Derivative stencil index out of bounds.");
    }
    #endif

    for (size_t II = 0; II < diff_array.size(); II++) {
        #ifdef NDEBUG
        diff_array[II] = diff_array[II] / (delta * scale_factor);
        #else
        diff_array.at(II) = diff_array.at(II) / (delta * scale_factor);
        #endif
    }

}
