#include <vector>
#include <string>
#include "../../differentiation_tools.hpp"
#include "../../constants.hpp"
#include "../../functions.hpp"

void longitude_derivative_at_point(
        const std::vector<double*> & deriv_vals,
        const std::vector<const std::vector<double>*> & fields,
        const dataset & source_data,
        const int Itime,
        const int Idepth,
        const int Ilat,
        const int Ilon,
        const std::vector<short int> & mask,
        const int order_of_deriv,
        const int diff_ord
        ){

    // Essentially just a wrapper to pass the longitude grid without having to manually specify it

    spher_derivative_at_point( deriv_vals, fields, source_data.longitude, "lon",
            source_data, Itime, Idepth, Ilat, Ilon, mask, order_of_deriv, diff_ord );

}
