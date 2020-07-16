#include <vector>
#include <math.h>
#include <cassert>
#include "../../differentiation_tools.hpp"
#include "../../constants.hpp"

void differentiation_vector(
        std::vector<double> & diff_array,
        const double delta,
        const int index,
        const int order_of_deriv,
        const int diff_ord
        ) {

    // Check that the diff_ord / order_of_deriv choice is valid
    assert( (order_of_deriv == 1) or (order_of_deriv == 2) );
    if (order_of_deriv == 1) {
        // First order derivatives can have 2nd, 4th, or 6th order convergence
        assert( (diff_ord == 2) or (diff_ord == 4) or (diff_ord == 6) );
    }
    if (order_of_deriv == 2) {
        // Second order derivatives can have 2nd, 4th, or 6th order convergence
        assert( (diff_ord == 2) or (diff_ord == 4) or (diff_ord == 6) );
    }

    std::vector<double>::iterator i0;
    i0 = diff_array.begin();

    double scale_factor = 1.;

    if (diff_ord == 2) {
        // Use second order
        switch (order_of_deriv) {
            case 1 :
                scale_factor = 1.;
                switch (index) {
                    case 0 :                // [-1.5,  2., -0.5]
                        diff_array.insert(i0, { -1.5,  2., -0.5});
                        break;
                    case 1 :                // [-0.5,  0.,  0.5]
                        diff_array.insert(i0, { -0.5,  0.,  0.5});
                        break;
                    case 2 :                // [ 0.5, -2.,  1.5]
                        diff_array.insert(i0, {  0.5, -2.,  1.5});
                        break;
                } break;
            case 2 :
                scale_factor = 1.;
                switch (index) {
                    case 0 :                // [ 2., -5.,  4., -1. ]
                        diff_array.insert(i0, {  2., -5.,  4., -1. });
                        break;
                    case 1 :                // [ 1., -2.,  1.,  0. ]
                        diff_array.insert(i0, {  1., -2.,  1.,  0. });
                        break;
                    case 2 :                // [ 0.,  1., -2.,  1. ]
                        diff_array.insert(i0, {  0.,  1., -2.,  1. });
                        break;
                    case 3 :                // [-1.,  4., -5.,  2. ]
                        diff_array.insert(i0, { -1.,  4., -5.,  2. });
                        break;
                } break;
        }
    } else if (diff_ord == 4) {
        // Use fourth order
        switch (order_of_deriv) {
            case 1 :
                scale_factor = 3.;
                switch (index) {
                    case 0 :              // [ -6.25,  12.,  -9.,    4., -0.75]  /  3
                        diff_array.insert(i0, {-6.25,  12.,  -9,     4., -0.75 });
                        break;
                    case 1 :              // [ -0.75, - 2.5,  4.5, - 1.5, 0.25]  /  3
                        diff_array.insert(i0, {-0.75, - 2.5,  4.5, - 1.5, 0.25 });
                        break;
                    case 2 :              // [  0.25, - 2.,   0.,    2., -0.25]  /  3
                        diff_array.insert(i0, { 0.25, - 2.,   0.,    2., -0.25 });
                        break;
                    case 3 :              // [ -0.25,   1.5, -4.5,   2.5, 0.75]  /  3
                        diff_array.insert(i0, {-0.25,   1.5, -4.5,   2.5, 0.75 });
                        break;
                    case 4 :              // [  0.75, - 4.,   9.,  -12.,  6.25]  /  3
                        diff_array.insert(i0, { 0.75, - 4.,   9.,  -12.,  6.25 });
                        break;
                } break;
            case 2 :
                scale_factor = 12.;
                switch (index) {
                    case 0 :              // [   45., -154.,  214., -156.,  61.,  -10. ]
                        diff_array.insert(i0, {  45., -154.,  214., -156.,  61.,  -10. });
                        break;
                    case 1 :              // [   10.,  -15.,   -4.,   14.,  -6.,    1. ]
                        diff_array.insert(i0, {  10.,  -15.,    4.,   14.,  -6.,    1. });
                        break;
                    case 2 :              // [   -1.,   16.,  -30.,   16.,   -1.,   0. ]
                        diff_array.insert(i0, {  -1.,   16.,  -30.,   16.,   -1.,   0. });
                        break;
                    case 3 :              // [   -0.,   -1.,   16.,  -30.,   16.,  -1. ]
                        diff_array.insert(i0, {  -0.,   -1.,   16.,  -30.,   16.,  -1. });
                        break;
                    case 4 :              // [    1.,   -6.,   14.,   -4.,  -15.,  10. ]
                        diff_array.insert(i0, {   1.,   -6.,   14.,   -4.,  -15.,  10. });
                        break;
                    case 5 :              // [  -10.,   61., -156.,  214., -154.,  45. ]
                        diff_array.insert(i0, { -10.,   61., -156.,  214., -154.,  45. });
                        break;
                } break;
        }
    } else if (diff_ord == 6) {
        // Use sixth order
        switch (order_of_deriv) {
            case 1 :
                scale_factor = 6.;
                switch (index) {
                    case 0 :               // [-14.7,  36. , -45. ,  40. , -22.5,   7.2,  -1. ]  /  6
                        diff_array.insert(i0, {-14.7,  36. , -45. ,  40. , -22.5,   7.2,  -1. });
                        break;
                    case 1 :               // [ -1. ,  -7.7,  15. , -10. ,   5. ,  -1.5,   0.2]  /  6
                        diff_array.insert(i0, { -1. ,  -7.7,  15. , -10. ,   5. ,  -1.5,   0.2});
                        break;
                    case 2 :               // [  0.2,  -2.4,  -3.5,   8. ,  -3. ,   0.8,  -0.1]  /  6
                        diff_array.insert(i0, {  0.2,  -2.4,  -3.5,   8. ,  -3. ,   0.8,  -0.1});
                        break;
                    case 3 :               // [ -0.1,   0.9,  -4.5,   0. ,   4.5,  -0.9,   0.1]  /  6
                        diff_array.insert(i0, { -0.1,   0.9,  -4.5,   0. ,   4.5,  -0.9,   0.1});
                        break;
                    case 4 :               // [  0.1,  -0.8,   3. ,  -8. ,   3.5,   2.4,  -0.2]  /  6
                        diff_array.insert(i0, {  0.1,  -0.8,   3. ,  -8. ,   3.5,   2.4,  -0.2});
                        break;
                    case 5 :               // [ -0.2,   1.5,  -5. ,  10. , -15. ,   7.7,   1. ]  /  6
                        diff_array.insert(i0, { -0.2,   1.5,  -5. ,  10. , -15. ,   7.7,   1. });
                        break;
                    case 6 :               // [  1. ,  -7.2,  22.5, -40. ,  45. , -36. ,  14.7]  /  6
                        diff_array.insert(i0, {  1. ,  -7.2,  22.5, -40. ,  45. , -36. ,  14.7});
                        break;
                } break;
            case 2 :
                scale_factor = 9.;
                switch (index) {
                    case 0 :               // [ 46.9, -200.7,  395.55,  -474.5,  369.,   -180.9,    50.95,  -6.3 ]
                        diff_array.insert(i0, { 46.9, -200.7,  395.55,  -474.5,  369.,   -180.9,    50.95,  -6.3  });
                        break;
                    case 1 :               // [  6.3,   -3.5,  -24.3,    42.75,  -33.5,    16.2,    -4.5,    0.55]
                        diff_array.insert(i0, {  6.3,   -3.5,  -24.3,    42.75,  -33.5,    16.2,    -4.5,    0.55 });
                        break;
                    case 2 :               // [ -0.55,  10.7,  -18.9,     6.5,     4.25,   -2.7,     0.8,   -0.1 ]
                        diff_array.insert(i0, { -0.55,  10.7,  -18.9,     6.5,     4.25,   -2.7,     0.8,   -0.1  });
                        break;
                    case 3 :               // [  0.1,   -1.35,  13.5,   -24.5,    13.5,    -1.35,    0.1,    0.  ]
                        diff_array.insert(i0, {  0.1,   -1.35,  13.5,   -24.5,    13.5,    -1.35,    0.1,    0.   });
                        break;
                    case 4 :               // [  0.,     0.1,   -1.35,   13.5,   -24.5,    13.5,    -1.35,   0.1 ]
                        diff_array.insert(i0, {  0.,     0.1,   -1.35,   13.5,   -24.5,    13.5,    -1.35,   0.1  });
                        break;
                    case 5 :               // [ -0.1,    0.8,   -2.7,     4.25,    6.5,   -18.9,    10.7,   -0.55]
                        diff_array.insert(i0, { -0.1,    0.8,   -2.7,     4.25,    6.5,   -18.9,    10.7,   -0.55 });
                        break;
                    case 6 :               // [  0.55,  -4.5,   16.2,   -33.5,    42.75,  -24.3,    -3.5,    6.3 ]
                        diff_array.insert(i0, {  0.55,  -4.5,   16.2,   -33.5,    42.75,  -24.3,    -3.5,    6.3  });
                        break;
                    case 7 :               // [ -6.3,   50.95, -180.9,  369.,    -474.5,  395.55, -200.7,   46.9 ]
                        diff_array.insert(i0, { -6.3,   50.95, -180.9,  369.,    -474.5,  395.55, -200.7,   46.9  });
                        break;
                } break;
        }
    }

    for (size_t II = 0; II < diff_array.size(); II++) {
        diff_array.at(II) = 
            diff_array.at(II) / ( pow(delta, order_of_deriv) * scale_factor);
    }

}
