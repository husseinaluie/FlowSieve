#include <vector>
#include <math.h>
#include "../../differentiation_tools.hpp"
#include "../../constants.hpp"

// differentaition_vector sets the vector needed to compute 
//   a fourth-order derivative on a five-point stencil.
// index specifies the index at which we want the derivative to be computed

void differentiation_vector(
        std::vector<double> & diff_array,
        const double delta,
        const int index,
        const int order_of_deriv,
        const int diff_ord
        ) {

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
        scale_factor = 6.;
        // Use sixth order
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
        }
    } else if (diff_ord == 8) {
        scale_factor = 42.;
        // Use eighth order
        //   Not well tested. Do not use.
        switch (index) {
            case 0 :               // [-1.1415e+02,  3.3600e+02, -5.8800e+02,  7.8400e+02, -7.3500e+02,  4.7040e+02, -1.9600e+02,  4.8000e+01, -5.2500e+00]
                diff_array.insert(i0, {-1.1415e+02,  3.3600e+02, -5.8800e+02,  7.8400e+02, -7.3500e+02,  4.7040e+02, -1.9600e+02,  4.8000e+01, -5.2500e+00});
                break;
            case 1 :               // [-5.2500e+00, -6.6900e+01,  1.4700e+02, -1.4700e+02,  1.2250e+02, -7.3500e+01,  2.9400e+01, -7.0000e+00,  7.5000e-01]
                diff_array.insert(i0, {-5.2500e+00, -6.6900e+01,  1.4700e+02, -1.4700e+02,  1.2250e+02, -7.3500e+01,  2.9400e+01, -7.0000e+00,  7.5000e-01});
                break;
            case 2 :               // [ 7.5000e-01, -1.2000e+01, -3.9900e+01,  8.4000e+01, -5.2500e+01,  2.8000e+01, -1.0500e+01,  2.4000e+00, -2.5000e-01]
                diff_array.insert(i0, { 7.5000e-01, -1.2000e+01, -3.9900e+01,  8.4000e+01, -5.2500e+01,  2.8000e+01, -1.0500e+01,  2.4000e+00, -2.5000e-01});
                break;
            case 3 :               // [-2.5000e-01,  3.0000e+00, -2.1000e+01, -1.8900e+01,  5.2500e+01, -2.1000e+01,  7.0000e+00, -1.5000e+00,  1.5000e-01]
                diff_array.insert(i0, {-2.5000e-01,  3.0000e+00, -2.1000e+01, -1.8900e+01,  5.2500e+01, -2.1000e+01,  7.0000e+00, -1.5000e+00,  1.5000e-01});
                break;
            case 4 :               // [ 1.5000e-01, -1.6000e+00,  8.4000e+00, -3.3600e+01, -2.1000e-13,  3.3600e+01, -8.4000e+00,  1.6000e+00, -1.5000e-01]
                diff_array.insert(i0, { 1.5000e-01, -1.6000e+00,  8.4000e+00, -3.3600e+01, -2.1000e-13,  3.3600e+01, -8.4000e+00,  1.6000e+00, -1.5000e-01});
                break;
            case 5 :               // [-1.5000e-01,  1.5000e+00, -7.0000e+00,  2.1000e+01, -5.2500e+01,  1.8900e+01,  2.1000e+01, -3.0000e+00,  2.5000e-01]
                diff_array.insert(i0, {-1.5000e-01,  1.5000e+00, -7.0000e+00,  2.1000e+01, -5.2500e+01,  1.8900e+01,  2.1000e+01, -3.0000e+00,  2.5000e-01});
                break;
            case 6 :               // [ 2.5000e-01, -2.4000e+00,  1.0500e+01, -2.8000e+01,  5.2500e+01, -8.4000e+01,  3.9900e+01,  1.2000e+01, -7.5000e-01]
                diff_array.insert(i0, { 2.5000e-01, -2.4000e+00,  1.0500e+01, -2.8000e+01,  5.2500e+01, -8.4000e+01,  3.9900e+01,  1.2000e+01, -7.5000e-01});
                break;
            case 7 :               // [-7.5000e-01,  7.0000e+00, -2.9400e+01,  7.3500e+01, -1.2250e+02,  1.4700e+02, -1.4700e+02,  6.6900e+01,  5.2500e+00]
                diff_array.insert(i0, {-7.5000e-01,  7.0000e+00, -2.9400e+01,  7.3500e+01, -1.2250e+02,  1.4700e+02, -1.4700e+02,  6.6900e+01,  5.2500e+00});
                break;
            case 8 :               // [ 5.2500e+00, -4.8000e+01,  1.9600e+02, -4.7040e+02,  7.3500e+02, -7.8400e+02,  5.8800e+02, -3.3600e+02,  1.1415e+02]
                diff_array.insert(i0, { 5.2500e+00, -4.8000e+01,  1.9600e+02, -4.7040e+02,  7.3500e+02, -7.8400e+02,  5.8800e+02, -3.3600e+02,  1.1415e+02});
                break;
        }
    }

    for (size_t II = 0; II < diff_array.size(); II++) {
        diff_array.at(II) = 
            diff_array.at(II) / ( pow(delta, order_of_deriv) * scale_factor);
    }

}
