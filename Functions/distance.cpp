#include <stdio.h>
#include <math.h>    
#include "../constants.hpp"
#include "../functions.hpp"


double arccos_approx( const double arg ) {

    const double x = std::fabs(arg);
    const double y = sqrt( 1 - x );

    // computing arccos is pretty expensive, especially
    // with how often we do it.
    // Instead, you optimized Horner polynomials
    // for three intervals of argument
    // [0,0.1], [0.1,0.8], [0.8,0.95], [0.95,1]
    //
    // Use symmetry arguments for x < 0

    double P = 0;
    if ( x < 0.1 ) {
        const double P0 =  1.0000000000000000,  
                     P1 =  0.2432988906302404, 
                     P2 = -1.0107836664831029, 
                     P3 =  1.8693217962750983, 
                     P4 = -1.4293766834995951,  
                     P5 =  0.4382603965522277;
        P = ((((P5 * y + P4) * y + P3) * y + P2) * y + P1) * y + P0;
    } else if ( x < 0.8 ) {
        const double P0 =  1.00000000000000000, 
                     P1 = -0.00350543916218819, 
                     P2 =  0.12356793324479017,
                     P3 = -0.19727421166830700,  
                     P4 =  0.55520703661267723, 
                     P5 = -0.87579987951110150,
                     P6 =  0.86700479385040619, 
                     P7 = -0.47558043880601630,  
                     P8 =  0.11709890844963859;
        /*
        const double P0 =  1.00854633e+00,
                     P1 = -1.07918165e-01,
                     P2 =  6.76310149e-01,
                     P3 = -1.85341341e+00, 
                     P4 =  3.62732306e+00,
                     P5 = -4.48910688e+00,
                     P6 =  3.49892720e+00,
                     P7 = -1.56128343e+00,
                     P8 =  3.11335834e-01;
        */
        P = (((((((P8 * y + P7) * y + P6) * y + P5) * y + P4) * y + P3) * y + P2) * y + P1) * y + P0;
    } else if ( x < 0.95 ) {
        const double P0 = 1.0000000000000000e+00, 
                     P1 = 3.7291411582717621e-05, 
                     P2 = 8.2768922348078589e-02, 
                     P3 = 3.3258850309590676e-03, 
                     P4 = 9.3148185535342851e-03, 
                     P5 = 1.2524738477274521e-02;
        P = ((((P5 * y + P4) * y + P3) * y + P2) * y + P1) * y + P0;
    } else {
        const double P0 =  1.0000000000000000e+00, 
                     P1 = -3.0731285022289612e-07, 
                     P2 =  8.3348776002847819e-02, 
                     P3 = -2.1426174881856075e-04, 
                     P4 =  1.9703699567281505e-02;
        P = (((P4 * y + P3) * y + P2) * y + P1) * y + P0;
    }
    const double output = sqrt(2.) * y * P;


    if ( arg < 0 ) { return M_PI - output; }
    else           { return output;        }
}

/*!
 * \brief Compute the distance (in metres) between two points in the domain.
 *
 * In spherical coordinates,computes the distance between two points on
 *   a spherical shell along a great circle.
 *   (see https://en.wikipedia.org/wiki/Great-circle_distance)
 * It should avoid floating point issues for the grid scales that
 *   we're considering. 
 *
 * The last two arguments (Llon and Llat) give the physical 
 *   length of the two dimensions. This is used in the case
 *   of periodic Cartesian grids. They are otherwise unused.
 *
 * @param[in]   lon1,lat1   coordinates for the first position
 * @param[in]   lon2,lat2   coordinates for the second position
 * @param[in]   Llon,Llat   physical length of the dimensions
 *
 * @returns returns the distance (in metres) between two points.
 *
 */
double distance(
        const double lon1,
        const double lat1,
        const double lon2,
        const double lat2,
        const double Llon,
        const double Llat,
        const long double cos_lat_1, 
        const long double cos_lat_2, 
        const long double sin_lat_1, 
        const long double sin_lat_2,
        const long double cos_lon_1, 
        const long double cos_lon_2, 
        const long double sin_lon_1, 
        const long double sin_lon_2,
        const bool pre_comp
        ) {

    // If they're the same point, just return zero. Otherwise, might encounter
    // some floating point issues in the inverse trigs (specifically the arccos call,
    // since rounding might push the argument to be slightly larger than 1.)
    if ( ( lon1 == lon2 ) and ( lat1 == lat2 ) ) { return 0.; }

    double distance;

    if (constants::CARTESIAN) {
        // If we're on a Cartesian grid, then just compute the straight-forward
        //   Cartesian distance. Account for periodicity if necessary
        double del_x, del_y;
        if (constants::PERIODIC_X) { del_x = fmin( fabs(lon1 - lon2), Llon - fabs(lon1 - lon2)); }
        else { del_x = lon1 - lon2; }

        if (constants::PERIODIC_Y) { del_y = fmin( fabs(lat1 - lat2), Llat - fabs(lat1 - lat2)); }
        else { del_y = lat1 - lat2; }

        distance = sqrt( pow(del_x, 2) + pow(del_y, 2) );

    } else {
        // If not on a Cartesian grid, then we're on a spherical grid.
        //   Compute distances along great circles.


        // Because of how frequently distances are computed,
        // these sin/cos calculations are cumulatively quite costly.
        // If the appropriate flags have been set, then we will pre-compute these
        // and pass them along
        const long double
            cos_lat1 = pre_comp ? cos_lat_1 : cos( lat1 ),
            cos_lat2 = pre_comp ? cos_lat_2 : cos( lat2 ),
            sin_lat1 = pre_comp ? sin_lat_1 : sin( lat1 ),
            sin_lat2 = pre_comp ? sin_lat_2 : sin( lat2 );

        // Since cos is even and cos(2pi - x) = cos(x), we don't need to worry 
        // about periodicity in computing Delta_lon for cos(Delta_lon)
        // for sin(Delta_lon), sin is odd and sin(2pi - x) = -sin(x), 
        // but since the result is being squared, it's not a concern either
        const long double Delta_lon = lon2 - lon1;
        const long double cos_Delta_lon = 
            pre_comp ? ( cos_lon_1 * cos_lon_2 + sin_lon_1 * sin_lon_2 ) : cos( Delta_lon );

        long double Delta_sigma;

        if (not(constants::USE_HIGH_PRECISION_DISTANCE)) {
            // This is cheaper, and so long as our distances are at least a couple metres, the floating-point issues shouldn't arise.
            const double acos_argument = sin_lat1 * sin_lat2 + cos_lat1 * cos_lat2 * cos_Delta_lon;
            
            // Handle some rounding cases when distance is nearly maximal
            if ( ( acos_argument < -1 ) and ( fabs(acos_argument + 1) < 1e-10 ) ) { 
                // If less than -1 but just a tiny bit
                Delta_sigma = M_PI;
            } else if ( ( acos_argument > 1 ) and ( acos_argument - 1 < 1e-10 ) ) {
                // If greater than +1 but just a tiny bit
                Delta_sigma = 0.;
            } else {
                // Otherwise, compute the arccos using an approximant (computationally much cheaper)
                // Error < 140 mm (versus true arccos calculation)
                Delta_sigma = arccos_approx( acos_argument );
                // Otherwise, use the full arccos calculator
                //Delta_sigma = acos( acos_argument );
            }
        } else {
            // If desired, use the more expensive distance calculator
            long double numer, denom;
            numer =   my_pow<2>(                                    cos_lat2 * sin(Delta_lon)) 
                    + my_pow<2>( cos_lat1 * sin_lat2  -  sin_lat1 * cos_lat2 * cos_Delta_lon );
            numer =  sqrt(numer);

            denom =        sin_lat1 * sin_lat2  +  cos_lat1 * cos_lat2 * cos_Delta_lon ;

            Delta_sigma = atan2(numer, denom);
        }
        distance = constants::R_earth * Delta_sigma;

    }

    return distance;
}
