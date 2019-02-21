#include "../functions.hpp"

double kernel(const double dist, const double scale) {
    
    if ( dist < (scale / 2) ) {
        return 1.;
    } else {
        return 0.;
    }
}
