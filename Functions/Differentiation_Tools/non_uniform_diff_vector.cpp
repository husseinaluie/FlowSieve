#include <vector>
#include <algorithm>
#include <math.h>
#include <cassert>
#include "../../differentiation_tools.hpp"
#include "../../constants.hpp"


void non_uniform_diff_vector(
        std::vector<double> & diff_array,
        const std::vector<double> & grid,
        const int Iref,
        const int LB,
        const int UB,
        const int diff_ord) {

    // For now, only works for diff_ord == 2
    assert(diff_ord == 2);

    std::vector<double>::iterator start;
    start = diff_array.begin();

    double scale_factor = 1.;
    double c1, c2, c3;
    double xn1, x0, xp1;

    switch (diff_ord) {
        case 2:
            xn1 = grid.at(LB);
            x0  = grid.at(LB+1);
            xp1 = grid.at(UB);
            if (Iref == LB) {
                c2 = -1. / ( 0.5 * pow(x0  - xn1, 2.) );
                c3 =  1. / ( 0.5 * pow(xp1 - xn1, 2.) );
                c1 = -(c3 + c2);
                scale_factor = (2/(xp1-xn1)) - (2/(x0-xn1));
            } else if (Iref == UB) {
                c1 = -1. / ( 0.5 * pow(xn1 - xp1, 2.) );
                c2 =  1. / ( 0.5 * pow(x0  - xp1, 2.) );
                c3 = -(c2 + c1);
                scale_factor = (2/(x0-xp1)) - (2/(xn1-xp1));
            } else {
                c1 = -1. / ( 0.5 * pow(xn1 - x0, 2.) );
                c3 =  1. / ( 0.5 * pow(xp1 - x0, 2.) );
                c2 = -(c3 + c1);
                scale_factor = (2/(xp1-x0)) - (2/(xn1-x0));
            }
            diff_array.insert(start, {c1, c2, c3});
            break;
    }

    for (size_t II = 0; II < diff_array.size(); II++) {
        diff_array.at(II) = diff_array.at(II) / scale_factor;
    }

}
