#include <vector>
#include "../../differentiation_tools.hpp"
#include "../../constants.hpp"

// differentaition_vector sets the vector needed to compute 
//   a fourth-order derivative on a five-point stencil.
// index specifies the index at which we want the derivative to be computed

void differentiation_vector(
        std::vector<double> & diff_array,   /**< [in] the array into which the differentiation coefficients should be placed */
        const double delta,                 /**< [in] the uniform grid spacing */
        const int index                     /**< [in] at which position in the five-point stencil do we want the derivative */
        ) {

    std::vector<double>::iterator start;
    start = diff_array.begin();

    double scale_factor = 1.;

    if (constants::DiffOrd == 2) {
        scale_factor = 1.;
        // Use second order
        switch (index) {
            case 0 :                   // [-1.5,  2., -0.5]
                diff_array.insert(start, { -1.5,  2., -0.5});
                break;
            case 1 :                   // [-0.5,  0.,  0.5]
                diff_array.insert(start, { -0.5,  0.,  0.5});
                break;
            case 2 :                   // [ 0.5, -2.,  1.5]
                diff_array.insert(start, {  0.5, -2.,  1.5});
                break;
        }
    } else if (constants::DiffOrd == 4) {
        scale_factor = 3.;
        // Use fourth order
        switch (index) {
            case 0 :                 // [ -6.25,  12.,  -9.,    4., -0.75]  /  3
                diff_array.insert(start, {-6.25,  12.,  -9,     4., -0.75 });
                break;
            case 1 :                 // [ -0.75, - 2.5,  4.5, - 1.5, 0.25]  /  3
                diff_array.insert(start, {-0.75, - 2.5,  4.5, - 1.5, 0.25 });
                break;
            case 2 :                 // [  0.25, - 2.,   0.,    2., -0.25]  /  3
                diff_array.insert(start, { 0.25, - 2.,   0.,    2., -0.25 });
                break;
            case 3 :                 // [ -0.25,   1.5, -4.5,   2.5, 0.75]  /  3
                diff_array.insert(start, {-0.25,   1.5, -4.5,   2.5, 0.75 });
                break;
            case 4 :                 // [  0.75, - 4.,   9.,  -12.,  6.25]  /  3
                diff_array.insert(start, { 0.75, - 4.,   9.,  -12.,  6.25 });
                break;
        }
    } else if (constants::DiffOrd == 6) {
        scale_factor = 6.;
        // Use sixth order
        switch (index) {
            case 0 :                  // [-14.7,  36. , -45. ,  40. , -22.5,   7.2,  -1. ]  /  6
                diff_array.insert(start, {-14.7,  36. , -45. ,  40. , -22.5,   7.2,  -1. });
                break;
            case 1 :                  // [ -1. ,  -7.7,  15. , -10. ,   5. ,  -1.5,   0.2]  /  6
                diff_array.insert(start, { -1. ,  -7.7,  15. , -10. ,   5. ,  -1.5,   0.2});
                break;
            case 2 :                  // [  0.2,  -2.4,  -3.5,   8. ,  -3. ,   0.8,  -0.1]  /  6
                diff_array.insert(start, {  0.2,  -2.4,  -3.5,   8. ,  -3. ,   0.8,  -0.1});
                break;
            case 3 :                  // [ -0.1,   0.9,  -4.5,   0. ,   4.5,  -0.9,   0.1]  /  6
                diff_array.insert(start, { -0.1,   0.9,  -4.5,   0. ,   4.5,  -0.9,   0.1});
                break;
            case 4 :                  // [  0.1,  -0.8,   3. ,  -8. ,   3.5,   2.4,  -0.2]  /  6
                diff_array.insert(start, {  0.1,  -0.8,   3. ,  -8. ,   3.5,   2.4,  -0.2});
                break;
            case 5 :                  // [ -0.2,   1.5,  -5. ,  10. , -15. ,   7.7,   1. ]  /  6
                diff_array.insert(start, { -0.2,   1.5,  -5. ,  10. , -15. ,   7.7,   1. });
                break;
            case 6 :                  // [  1. ,  -7.2,  22.5, -40. ,  45. , -36. ,  14.7]  /  6
                diff_array.insert(start, {  1. ,  -7.2,  22.5, -40. ,  45. , -36. ,  14.7});
                break;
        }
    } else if (constants::DiffOrd == 8) {
        scale_factor = 42.;
        // Use eighth order
        //   Not well tested. Do not use.
        switch (index) {
            case 0 :                  // [-1.1415e+02,  3.3600e+02, -5.8800e+02,  7.8400e+02, -7.3500e+02,  4.7040e+02, -1.9600e+02,  4.8000e+01, -5.2500e+00]
                diff_array.insert(start, {-1.1415e+02,  3.3600e+02, -5.8800e+02,  7.8400e+02, -7.3500e+02,  4.7040e+02, -1.9600e+02,  4.8000e+01, -5.2500e+00});
                break;
            case 1 :                  // [-5.2500e+00, -6.6900e+01,  1.4700e+02, -1.4700e+02,  1.2250e+02, -7.3500e+01,  2.9400e+01, -7.0000e+00,  7.5000e-01]
                diff_array.insert(start, {-5.2500e+00, -6.6900e+01,  1.4700e+02, -1.4700e+02,  1.2250e+02, -7.3500e+01,  2.9400e+01, -7.0000e+00,  7.5000e-01});
                break;
            case 2 :                  // [ 7.5000e-01, -1.2000e+01, -3.9900e+01,  8.4000e+01, -5.2500e+01,  2.8000e+01, -1.0500e+01,  2.4000e+00, -2.5000e-01]
                diff_array.insert(start, { 7.5000e-01, -1.2000e+01, -3.9900e+01,  8.4000e+01, -5.2500e+01,  2.8000e+01, -1.0500e+01,  2.4000e+00, -2.5000e-01});
                break;
            case 3 :                  // [-2.5000e-01,  3.0000e+00, -2.1000e+01, -1.8900e+01,  5.2500e+01, -2.1000e+01,  7.0000e+00, -1.5000e+00,  1.5000e-01]
                diff_array.insert(start, {-2.5000e-01,  3.0000e+00, -2.1000e+01, -1.8900e+01,  5.2500e+01, -2.1000e+01,  7.0000e+00, -1.5000e+00,  1.5000e-01});
                break;
            case 4 :                  // [ 1.5000e-01, -1.6000e+00,  8.4000e+00, -3.3600e+01, -2.1000e-13,  3.3600e+01, -8.4000e+00,  1.6000e+00, -1.5000e-01]
                diff_array.insert(start, { 1.5000e-01, -1.6000e+00,  8.4000e+00, -3.3600e+01, -2.1000e-13,  3.3600e+01, -8.4000e+00,  1.6000e+00, -1.5000e-01});
                break;
            case 5 :                  // [-1.5000e-01,  1.5000e+00, -7.0000e+00,  2.1000e+01, -5.2500e+01,  1.8900e+01,  2.1000e+01, -3.0000e+00,  2.5000e-01]
                diff_array.insert(start, {-1.5000e-01,  1.5000e+00, -7.0000e+00,  2.1000e+01, -5.2500e+01,  1.8900e+01,  2.1000e+01, -3.0000e+00,  2.5000e-01});
                break;
            case 6 :                  // [ 2.5000e-01, -2.4000e+00,  1.0500e+01, -2.8000e+01,  5.2500e+01, -8.4000e+01,  3.9900e+01,  1.2000e+01, -7.5000e-01]
                diff_array.insert(start, { 2.5000e-01, -2.4000e+00,  1.0500e+01, -2.8000e+01,  5.2500e+01, -8.4000e+01,  3.9900e+01,  1.2000e+01, -7.5000e-01});
                break;
            case 7 :                  // [-7.5000e-01,  7.0000e+00, -2.9400e+01,  7.3500e+01, -1.2250e+02,  1.4700e+02, -1.4700e+02,  6.6900e+01,  5.2500e+00]
                diff_array.insert(start, {-7.5000e-01,  7.0000e+00, -2.9400e+01,  7.3500e+01, -1.2250e+02,  1.4700e+02, -1.4700e+02,  6.6900e+01,  5.2500e+00});
                break;
            case 8 :                  // [ 5.2500e+00, -4.8000e+01,  1.9600e+02, -4.7040e+02,  7.3500e+02, -7.8400e+02,  5.8800e+02, -3.3600e+02,  1.1415e+02]
                diff_array.insert(start, { 5.2500e+00, -4.8000e+01,  1.9600e+02, -4.7040e+02,  7.3500e+02, -7.8400e+02,  5.8800e+02, -3.3600e+02,  1.1415e+02});
                break;
        }
    }

    for (int II = 0; II < constants::DiffOrd+1; II++) {
        diff_array.at(II) = diff_array.at(II) / ( delta * scale_factor);
    }

}
