#include <fenv.h>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <math.h>
#include <vector>

#include "../functions.hpp"
#include "../constants.hpp"


// Snyder, J. P. (1987). Map Projections: A Working Manual. https://doi.org/10.3133/pp1395
// page 173
void near_sided_project( 
       double & proj_x,
       double & proj_y,
       const double & pt_lon,
       const double & pt_lat,
       const double & centre_lon,
       const double & centre_lat,
       const double & viewing_elev
        ) {

    const double    cos_phi = cos( pt_lat ),
                    sin_phi = sin( pt_lat ),
                    cos_phi0 = cos( centre_lat ),
                    sin_phi0 = sin( centre_lat ),
                    sin_lam = sin( pt_lon - centre_lon ),
                    cos_lam = cos( pt_lon - centre_lon );

    const double    P = 1 + viewing_elev / constants::R_earth;

    const double    cos_c = sin_phi0 * sin_phi + cos_phi0 * cos_phi * cos_lam;

    if (cos_c < 1/P ) {
        proj_x = 1e100;
        proj_y = 1e100;
        return;
    }
    
    const double k_prime = ( P - 1 ) / ( P - cos_c );

    proj_x = constants::R_earth * k_prime * cos_phi * sin_lam;
    proj_y = constants::R_earth * k_prime * ( cos_phi0 * sin_phi - sin_phi0 * cos_phi * cos_lam );

};

