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

    int order_of_deriv = 1;

    std::vector<double>::iterator start;
    start = diff_array.begin();

    double scale_factor = 1.;
    double c1, c2, c3, c4;
    double xn1, x0, xp1;

    if (diff_ord == 2) {  // second order accurate
        switch (order_of_deriv) {
            case 1: // first derivative
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
            case 2: // second derivative
                if (Iref == LB) {

                    double dp1 = grid.at(Iref+1) - grid.at(Iref);
                    double dp2 = grid.at(Iref+2) - grid.at(Iref);
                    double dp3 = grid.at(Iref+3) - grid.at(Iref);

                    double Ap1 = ( pow(1./dp1,1) - pow(1./dp3,1) );
                    double Ap2 = ( pow(1./dp2,1) - pow(1./dp3,1) );

                    double Bp1 = ( pow(1./dp1,2) - pow(1./dp3,2) );
                    double Bp2 = ( pow(1./dp2,2) - pow(1./dp3,2) );

                    double Cp1 = ( pow(1./dp1,3) - pow(1./dp3,3) );
                    double Cp2 = ( pow(1./dp2,3) - pow(1./dp3,3) );

                    c1 = - ( (Cp2 / Bp2) - (Cp1 / Bp1) );
                    c2 = - 1. / ( Bp1 * pow(dp1, 3) );
                    c3 =   1. / ( Bp2 * pow(dp2, 3) );
                    c4 = - ( (1./Bp2) - (1./Bp1) ) * pow(1./dp3,3);

                    scale_factor = 0.5 * ( (Ap2/Bp2) - (Ap1/Bp1) );

                } else if (Iref == LB + 1) {

                    double dn1 = grid.at(Iref-1) - grid.at(Iref);
                    double dp1 = grid.at(Iref+1) - grid.at(Iref);
                    double dp2 = grid.at(Iref+2) - grid.at(Iref);

                    double An1 = ( pow(1./dn1,1) - pow(1./dp2,1) );
                    double Ap1 = ( pow(1./dp1,1) - pow(1./dp2,1) );

                    double Bn1 = ( pow(1./dn1,2) - pow(1./dp2,2) );
                    double Bp1 = ( pow(1./dp1,2) - pow(1./dp2,2) );

                    double Cn1 = ( pow(1./dn1,3) - pow(1./dp2,3) );
                    double Cp1 = ( pow(1./dp1,3) - pow(1./dp2,3) );

                    c1 =   1. / ( Bn1 * pow(dn1, 3) );
                    c2 = - ( (Cn1 / Bn1) - (Cp1 / Bp1) );
                    c3 = - 1. / ( Bp1 * pow(dp1, 3) );
                    c4 = - ( (1./Bn1) - (1./Bp1) ) * pow(1./dp2,3);

                    scale_factor = 0.5 * ( (An1/Bn1) - (Ap1/Bp1) );

                } else if (Iref == LB + 2) {

                } else if (Iref == UB) {

                }
                diff_array.insert(start, {c1, c2, c3, c4});
                break;
        }
    }

    for (size_t II = 0; II < diff_array.size(); II++) {
        diff_array.at(II) = diff_array.at(II) / scale_factor;
    }

}
