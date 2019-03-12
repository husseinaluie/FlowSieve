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

    if (constants::DiffOrd == 2) {
        // Use second order
        if (index == 0) {
            // [-1.5,  2. , -0.5]
            diff_array.at(0) = -1.5;
            diff_array.at(1) =  2. ;
            diff_array.at(2) = -0.5;
        } else if (index == 1) {
            // [-0.5,  0. ,  0.5]
            diff_array.at(0) = -0.5;
            diff_array.at(1) =  0. ;
            diff_array.at(2) =  0.5;
        } else if (index == 2) {
            // [0.5, -2. ,  1.5]
            diff_array.at(0) =  0.5;
            diff_array.at(1) = -2. ;
            diff_array.at(2) =  1.5;
        }
    } else if (constants::DiffOrd == 4) {
        // Use fourth order
        if (index == 0) {
            // [ -6.25,  12.  ,  -9.  ,   4.  ,  -0.75]  /  3
            diff_array.at(0) = -6.25 / 3;
            diff_array.at(1) =  4.      ;
            diff_array.at(2) = -3.      ;
            diff_array.at(3) =  4.   / 3;
            diff_array.at(4) = -0.75 / 3;
        } else if (index == 1) {
            // [ -0.75,  -2.5 ,   4.5 ,  -1.5 ,   0.25]  /  3
            diff_array.at(0) = -0.75 / 3;
            diff_array.at(1) = -2.5  / 3;
            diff_array.at(2) =  4.5  / 3;
            diff_array.at(3) = -1.5  / 3;
            diff_array.at(4) =  0.25 / 3;
        } else if (index == 2) {
            // [  0.25,  -2.  ,   0.  ,   2.  ,  -0.25]  /  3
            diff_array.at(0) =  0.25 / 3;
            diff_array.at(1) = -2.   / 3;
            diff_array.at(2) =  0.      ;
            diff_array.at(3) =  2.   / 3;
            diff_array.at(4) = -0.25 / 3;
        } else if (index == 3) {
            // [ -0.25,   1.5 ,  -4.5 ,   2.5 ,   0.75]  /  3
            diff_array.at(0) = -0.25 / 3;
            diff_array.at(1) =  1.5  / 3;
            diff_array.at(2) = -4.5  / 3;
            diff_array.at(3) =  2.5  / 3;
            diff_array.at(4) =  0.75 / 3;
        } else if (index == 4) {
            // [  0.75,  -4.  ,   9.  , -12.  ,   6.25]  /  3
            diff_array.at(0) =  0.75 / 3;
            diff_array.at(1) = -4.   / 3;
            diff_array.at(2) =  3.      ;
            diff_array.at(3) = -4.      ;
            diff_array.at(4) =  6.25 / 3;
        }
    } else if (constants::DiffOrd == 6) {
        // Use sixth order
        if (index == 0) {
            // [-14.7,  36. , -45. ,  40. , -22.5,   7.2,  -1. ]  /  6
            diff_array.at(0) = -14.7 / 6;
            diff_array.at(1) =  36.  / 6;
            diff_array.at(2) = -45.  / 6;
            diff_array.at(3) =  40.  / 6;
            diff_array.at(4) = -22.5 / 6;
            diff_array.at(5) =   7.2 / 6;
            diff_array.at(6) = - 1.  / 6;
        } else if (index == 1) {
            // [ -1. ,  -7.7,  15. , -10. ,   5. ,  -1.5,   0.2]  /  6
            diff_array.at(0) = - 1.  / 6;
            diff_array.at(1) = - 7.7 / 6;
            diff_array.at(2) =  15.  / 6;
            diff_array.at(3) = -10.  / 6;
            diff_array.at(4) =   5.  / 6;
            diff_array.at(5) = - 1.5 / 6;
            diff_array.at(6) =   0.2 / 6;
        } else if (index == 2) {
            // [  0.2,  -2.4,  -3.5,   8. ,  -3. ,   0.8,  -0.1]  /  6
            diff_array.at(0) =   0.2 / 6;
            diff_array.at(1) = - 2.4 / 6;
            diff_array.at(2) = - 3.5 / 6;
            diff_array.at(3) =   8.  / 6;
            diff_array.at(4) = - 3.  / 6;
            diff_array.at(5) =   0.8 / 6;
            diff_array.at(6) = - 0.1 / 6;
        } else if (index == 3) {
            // [ -0.1,   0.9,  -4.5,   0. ,   4.5,  -0.9,   0.1]  /  6
            diff_array.at(0) = - 0.1 / 6;
            diff_array.at(1) =   0.9 / 6;
            diff_array.at(2) = - 4.5 / 6;
            diff_array.at(3) =   0.  / 6;
            diff_array.at(4) =   4.5 / 6;
            diff_array.at(5) = - 0.9 / 6;
            diff_array.at(6) =   0.1 / 6;
        } else if (index == 4) {
            // [  0.1,  -0.8,   3. ,  -8. ,   3.5,   2.4,  -0.2]  /  6
            diff_array.at(0) =   0.1 / 6;
            diff_array.at(1) = - 0.8 / 6;
            diff_array.at(2) =   3.  / 6;
            diff_array.at(3) = - 8.  / 6;
            diff_array.at(4) =   3.5 / 6;
            diff_array.at(5) =   2.4 / 6;
            diff_array.at(6) = - 0.2 / 6;
        } else if (index == 5) {
            // [ -0.2,   1.5,  -5. ,  10. , -15. ,   7.7,   1. ]  /  6
            diff_array.at(0) = - 0.2 / 6;
            diff_array.at(1) =   1.5 / 6;
            diff_array.at(2) = - 5.  / 6;
            diff_array.at(3) =  10.  / 6;
            diff_array.at(4) = -15.  / 6;
            diff_array.at(5) =   7.7 / 6;
            diff_array.at(6) =   1.  / 6;
        } else if (index == 6) {
            // [  1. ,  -7.2,  22.5, -40. ,  45. , -36. ,  14.7]  /  6
            diff_array.at(0) =   1.0 / 6;
            diff_array.at(1) = - 7.2 / 6;
            diff_array.at(2) =  22.5 / 6;
            diff_array.at(3) = -40.  / 6;
            diff_array.at(4) =  45.  / 6;
            diff_array.at(5) = -36.  / 6;
            diff_array.at(6) =  14.7 / 6;
        }
    }

    for (int II = 0; II < constants::DiffOrd+1; II++) {
        diff_array.at(II) = diff_array.at(II) / delta;
    }

}
