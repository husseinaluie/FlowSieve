#include <math.h>
#include "../functions.hpp"

/*
 * \brief Applies periodic adjustment on [-N,2N] to [0,N]
 *
 * The older formulation (( LON % Nlon + Nlon ) % Nlon;) 
 * is actually pretty expensive consider how frequenly it
 * is called [% is not cheap].
 *
 * This form is optimized to be cheaper, but only takes
 * inputs on [-N,2N]
 *
 */
int apply_periodicity( const int I, const int N ){

    if (I < 0)       { return I + N; }
    else if (I >= N) { return I - N; }
    else             { return I;     }

}
