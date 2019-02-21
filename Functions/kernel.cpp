#include "../functions.hpp"

double kernel(const double dist, const double scale) {
    
    double kern;

    if (dist < scale / 2) {
        kern = 1.;
    } else {
        kern = 0.;
    }

    return kern;
}
