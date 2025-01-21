#include "../constants.hpp"
#include "../functions.hpp"
#include "../preprocess.hpp"
#include <algorithm>
#include <vector>
#include <omp.h>
#include <math.h>
#include <unordered_map>

#include <Eigen/QR>


std::tuple< double, double, double > midpoint_on_sphere(
        const std::tuple< double, double, double > & P,
        const std::tuple< double, double, double > & Q
        ){

    const double x = ( std::get<0>(P) + std::get<0>(Q) ) / 2.,
                 y = ( std::get<1>(P) + std::get<1>(Q) ) / 2.,
                 z = ( std::get<2>(P) + std::get<2>(Q) ) / 2.;

    const double r = sqrt( x*x + y*y + z*z );

    #if DEBUG > 0
    assert( r > 0 );
    #endif

    std::tuple< double, double, double > midpoint = { x/r, y/r, z/r };
    return midpoint;

}

std::tuple< double, double, double > midpoint_on_sphere(
        const std::tuple< double, double, double > & P,
        const std::tuple< double, double, double > & Q,
        const std::tuple< double, double, double > & R
        ){

    const double x = ( std::get<0>(P) + std::get<0>(Q) + std::get<0>(R) ) / 3.,
                 y = ( std::get<1>(P) + std::get<1>(Q) + std::get<1>(R) ) / 3.,
                 z = ( std::get<2>(P) + std::get<2>(Q) + std::get<2>(R) ) / 3.;

    const double r = sqrt( x*x + y*y + z*z );

    #if DEBUG > 0
    assert( r > 0 );
    #endif

    std::tuple< double, double, double > midpoint = { x/r, y/r, z/r };
    return midpoint;

}

double SphericalTriangleArea( 
            std::tuple< double, double, double > lons,
            std::tuple< double, double, double > lats
        ) {

    // Length of each side of the triangle (in rad)
    const double a = distance( std::get<0>(lons), std::get<0>(lats), 
                               std::get<1>(lons), std::get<1>(lats) 
                               ) / constants::R_earth,
                 b = distance( std::get<1>(lons), std::get<1>(lats), 
                               std::get<2>(lons), std::get<2>(lats) 
                               ) / constants::R_earth,
                 c = distance( std::get<2>(lons), std::get<2>(lats), 
                               std::get<0>(lons), std::get<0>(lats) 
                               ) / constants::R_earth;

    assert( a > 0 );
    assert( b > 0 );
    assert( c > 0 );

    const double cos_a = cos(a), cos_b = cos(b), cos_c = cos(c),
                 sin_a = sin(a), sin_b = sin(b), sin_c = sin(c);

    // Interiour angles of the triangle (in rad)
    //  using arg_? to handle floating point errors
    const double arg_a = ( cos_a - cos_b * cos_c ) / ( sin_b * sin_c ),
                 arg_b = ( cos_b - cos_c * cos_a ) / ( sin_c * sin_a ),
                 arg_c = ( cos_c - cos_a * cos_b ) / ( sin_a * sin_b );
    const double A = acos( (arg_a < -1) ? -1 : (arg_a > 1) ? 1 : arg_a ),
                 B = acos( (arg_b < -1) ? -1 : (arg_b > 1) ? 1 : arg_b ),
                 C = acos( (arg_c < -1) ? -1 : (arg_c > 1) ? 1 : arg_c );

    // Spherical excess of the triangle
    const double E = A + B + C - M_PI;

    if ( (E > M_PI) or std::isnan(E) ) {
        throw std::runtime_error("Invalid cell area");
    }

    return pow( constants::R_earth, 2. ) * E;
}

struct pair_hash {
    std::size_t operator() ( const std::pair<size_t, size_t> &p ) const {
        const size_t A = p.first;
        const size_t B = p.second;
        const size_t hash_value = (A >= B) ? A*A + A*B : A + B * B;
        return hash_value;
    }
};


void BuildPolyhedralGrid(
        dataset & polyhedral_grid,
        const size_t target_num_points
        ) {

    const unsigned int tetra_order = round( log2( target_num_points / 4 ) / 2. );
    const size_t tetrahedral_npts = 4 * pow( 4, tetra_order );
    const size_t tetrahedral_size_error = ( tetrahedral_npts > target_num_points ) ? tetrahedral_npts - target_num_points : target_num_points - tetrahedral_npts;
    //const size_t tetrahedral_size_error = 1;

    const unsigned int octa_order =  round( log2( target_num_points / 8 ) / 2. );
    const size_t octahedral_npts = 8 * pow( 4, octa_order );
    const size_t octahedral_size_error = ( octahedral_npts > target_num_points ) ? octahedral_npts - target_num_points : target_num_points - octahedral_npts;
    //const size_t octahedral_size_error = 1;

    const unsigned int icosa_order = round( log2( target_num_points / 20 ) / 2. );
    const size_t icosahedral_npts = 20 * pow( 4, icosa_order );
    const size_t icosahedral_size_error = ( icosahedral_npts > target_num_points ) ? icosahedral_npts - target_num_points : target_num_points - icosahedral_npts;
    //const size_t icosahedral_size_error = 0;

    unsigned int refine_rounds, base_faces, base_vertices;

    if ( tetrahedral_size_error < octahedral_size_error ) {
        if ( tetrahedral_size_error < icosahedral_size_error ) {
            // Tetrehedral base is best
            refine_rounds = tetra_order;
            base_faces = 4;
            base_vertices = 4;
        } else {
            // Icosahedral base is best
            refine_rounds = icosa_order;
            base_faces = 20;
            base_vertices = 12;
        }
    } else if ( icosahedral_size_error < octahedral_size_error ) {
        // Icosahedral base is best
        refine_rounds = icosa_order;
        base_faces = 20;
        base_vertices = 12;
    } else {
        // Octahedral base is best
        refine_rounds = octa_order;
        base_faces = 8;
        base_vertices = 6;
    }
    const size_t num_faces = base_faces * pow( (size_t)4, refine_rounds ),
                 num_vertices = base_vertices + (base_faces/2) * ( pow( (size_t)4, refine_rounds ) - 1 );

    #if DEBUG >= 1
    fprintf( stdout, "Starting with %d faces and %d vertices, we will refine %d times to get %'zu faces and %'zu vertices.\n", base_faces, base_vertices, refine_rounds, num_faces, num_vertices );
    #endif

    std::vector< std::tuple< double, double, double > > vertices(num_vertices);
    std::vector< std::tuple< size_t, size_t, size_t > > faces(num_faces);


    if (base_vertices == 4) {
        #if DEBUG >= 2
        fprintf( stdout, "BuildPolyhedralGrid: Using a tetrahedral base.\n" );
        #endif
        // Tetrahedra base
        
        // This rotation has a guaranteed point at the pole
        //vertices[0] = std::tuple<double,double,double> {  sqrt(8./9),           0, -1./3 };
        //vertices[1] = std::tuple<double,double,double> { -sqrt(2./9),  sqrt(2./3), -1./3 };
        //vertices[2] = std::tuple<double,double,double> { -sqrt(2./9), -sqrt(2./3), -1./3 };
        //vertices[3] = std::tuple<double,double,double> {           0,           0,  1 };
        
        // This rotation of the base avoids a guaranteed pole (does not guaranteed the absence
        //  of one though)
        vertices[0] = std::tuple<double,double,double> { -1./3,           0,  sqrt(8./9) };
        vertices[1] = std::tuple<double,double,double> { -1./3,  sqrt(2./3), -sqrt(2./9) };
        vertices[2] = std::tuple<double,double,double> { -1./3, -sqrt(2./3), -sqrt(2./9) };
        vertices[3] = std::tuple<double,double,double> {  1,              0,           0 };

        faces[0] = std::tuple< size_t, size_t, size_t > { 0, 1, 2 };
        faces[1] = std::tuple< size_t, size_t, size_t > { 0, 1, 3 };
        faces[2] = std::tuple< size_t, size_t, size_t > { 0, 2, 3 };
        faces[3] = std::tuple< size_t, size_t, size_t > { 1, 2, 3 };
    } else if (base_vertices == 6) {
        #if DEBUG >= 2
        fprintf( stdout, "BuildPolyhedralGrid: Using an octahedral base.\n" );
        #endif
        // Octahedral base
        vertices[0] = std::tuple<double,double,double> {  1,  0,  0 };
        vertices[1] = std::tuple<double,double,double> { -1,  0,  0 };
        vertices[2] = std::tuple<double,double,double> {  0,  1,  0 };
        vertices[3] = std::tuple<double,double,double> {  0, -1,  0 };
        vertices[4] = std::tuple<double,double,double> {  0,  0,  1 };
        vertices[5] = std::tuple<double,double,double> {  0,  0, -1 };

        faces[0] = std::tuple< size_t, size_t, size_t > { 0, 2, 4 };
        faces[1] = std::tuple< size_t, size_t, size_t > { 0, 2, 5 };
        faces[2] = std::tuple< size_t, size_t, size_t > { 0, 3, 4 };
        faces[3] = std::tuple< size_t, size_t, size_t > { 0, 3, 5 };
        faces[4] = std::tuple< size_t, size_t, size_t > { 1, 2, 4 };
        faces[5] = std::tuple< size_t, size_t, size_t > { 1, 2, 5 };
        faces[6] = std::tuple< size_t, size_t, size_t > { 1, 3, 4 };
        faces[7] = std::tuple< size_t, size_t, size_t > { 1, 3, 5 };
    } else if (base_vertices == 12) {
        #if DEBUG >= 2
        fprintf( stdout, "BuildPolyhedralGrid: Using an icosahedral base.\n" );
        #endif
        // Icosahedral base
        const double phi = 0.5 * ( 1. + sqrt(5.) );
        const double p = sqrt( 1. / ( phi*phi + 1. ) );
        const double q = phi * p;

        vertices[0] = std::tuple<double,double,double> { -p,  0,  q };
        vertices[1] = std::tuple<double,double,double> {  p,  0,  q };
        vertices[2] = std::tuple<double,double,double> { -p,  0, -q };
        vertices[3] = std::tuple<double,double,double> {  p,  0, -q };

        vertices[4] = std::tuple<double,double,double> {  0,  q,  p };
        vertices[5] = std::tuple<double,double,double> {  0,  q, -p };
        vertices[6] = std::tuple<double,double,double> {  0, -q,  p };
        vertices[7] = std::tuple<double,double,double> {  0, -q, -p };

        vertices[ 8] = std::tuple<double,double,double> {  q,  p,  0 };
        vertices[ 9] = std::tuple<double,double,double> { -q,  p,  0 };
        vertices[10] = std::tuple<double,double,double> {  q, -p,  0 };
        vertices[11] = std::tuple<double,double,double> { -q, -p,  0 };

        faces[ 0] = std::tuple< size_t, size_t, size_t > {  0,  4,  1 };
        faces[ 1] = std::tuple< size_t, size_t, size_t > {  0,  9,  4 };
        faces[ 2] = std::tuple< size_t, size_t, size_t > {  9,  4,  5 };
        faces[ 3] = std::tuple< size_t, size_t, size_t > {  4,  5,  8 };
        faces[ 4] = std::tuple< size_t, size_t, size_t > {  4,  8,  1 };
        faces[ 5] = std::tuple< size_t, size_t, size_t > {  8, 10,  1 };
        faces[ 6] = std::tuple< size_t, size_t, size_t > {  8,  3, 10 };
        faces[ 7] = std::tuple< size_t, size_t, size_t > {  5,  3,  8 };
        faces[ 8] = std::tuple< size_t, size_t, size_t > {  5,  2,  3 };
        faces[ 9] = std::tuple< size_t, size_t, size_t > {  2,  7,  3 };
        faces[10] = std::tuple< size_t, size_t, size_t > {  7, 10,  3 };
        faces[11] = std::tuple< size_t, size_t, size_t > {  7,  6, 10 };
        faces[12] = std::tuple< size_t, size_t, size_t > {  7, 11,  6 };
        faces[13] = std::tuple< size_t, size_t, size_t > { 11,  0,  6 };
        faces[14] = std::tuple< size_t, size_t, size_t > {  0,  1,  6 };
        faces[15] = std::tuple< size_t, size_t, size_t > {  6,  1, 10 };
        faces[16] = std::tuple< size_t, size_t, size_t > {  9,  0, 11 };
        faces[17] = std::tuple< size_t, size_t, size_t > {  9, 11,  2 };
        faces[18] = std::tuple< size_t, size_t, size_t > {  9,  2,  5 };
        faces[19] = std::tuple< size_t, size_t, size_t > {  7,  2, 11 };
    } else {
        throw std::logic_error("Tried to initialize a polyhedron with an invalid number of vertices.");
    }


    //
    //// Apply the requested rounds of refinement to build the mesh of appropriate size
    //

    for ( unsigned int refine_round = 0; refine_round < refine_rounds; refine_round++ ) {

        std::unordered_map< std::pair<size_t, size_t>, size_t, pair_hash > pair_children;

        // Get the index max for this iteration and previous iteration
        const size_t num_faces_prev_round = base_faces * pow( (size_t)4, refine_round );
        const size_t num_vertices_prev_round = base_vertices + (base_faces/2) * ( pow( (size_t)4, refine_round ) - 1 );

        size_t new_index = num_vertices_prev_round,
               new_face = num_faces_prev_round;
        for ( size_t Iface = 0; Iface < num_faces_prev_round; Iface++ ) {

            // We're going to get the three new vertices made by
            // subdividing each side of the triangular face
            std::vector< size_t > new_vertex_indices(3);
            for ( int II = 0; II < 3; II++ ) {
                const int JJ = (II+1) % 3;
                // The <I> to get must be a compile-time constant... so boo. Hardcode the three cases
                size_t parent_I =     (II == 0) ? std::get<0>(faces[Iface])
                                    : (II == 1) ? std::get<1>(faces[Iface])
                                    :             std::get<2>(faces[Iface]);
                size_t parent_J =     (JJ == 0) ? std::get<0>(faces[Iface])
                                    : (JJ == 1) ? std::get<1>(faces[Iface])
                                    :             std::get<2>(faces[Iface]);

                // Always record parent pairs in increasing order
                std::pair<size_t,size_t> parents = 
                    ( parent_I < parent_J )
                    ? 
                    std::pair<size_t,size_t>(parent_I,parent_J)
                    : 
                    std::pair<size_t,size_t>(parent_J,parent_I)
                    ;

                size_t child_index;
                if ( pair_children.count( parents ) > 0 ) {
                    child_index = pair_children[parents];
                } else {
                    child_index = new_index;
                    new_index++;
                    pair_children[parents] = child_index;

                    std::tuple<double, double, double> new_vertex = midpoint_on_sphere( vertices[parent_I], vertices[parent_J] );

                    vertices[child_index] = new_vertex;
                }
                new_vertex_indices[II] = child_index;
            }

            // And record the four new faces embedded in the previous face
            std::get<0>( faces[new_face] ) = std::get<0>(faces[Iface]);
            std::get<1>( faces[new_face] ) = new_vertex_indices[0];
            std::get<2>( faces[new_face] ) = new_vertex_indices[2];
            new_face++;

            std::get<0>( faces[new_face] ) = std::get<1>(faces[Iface]);
            std::get<1>( faces[new_face] ) = new_vertex_indices[0];
            std::get<2>( faces[new_face] ) = new_vertex_indices[1];
            new_face++;

            std::get<0>( faces[new_face] ) = std::get<2>(faces[Iface]);
            std::get<1>( faces[new_face] ) = new_vertex_indices[1];
            std::get<2>( faces[new_face] ) = new_vertex_indices[2];
            new_face++;

            // and overwrite the original face with this new one
            std::get<0>( faces[Iface] ) = new_vertex_indices[0];
            std::get<1>( faces[Iface] ) = new_vertex_indices[1];
            std::get<2>( faces[Iface] ) = new_vertex_indices[2];
        }
    }



    //
    //// Convert to the dual to get the actual working map
    //
    const size_t num_faces_final = faces.size();
    const size_t size_neighbourhood = constants::ADJACENCY_SIZE+1;
    #if DEBUG >= 2
    fprintf( stdout, "BuildPolyhedralGrid: The final polyhedron has %'zu faces.\n", num_faces_final );
    #endif

    std::vector<double> &longitude = polyhedral_grid.longitude,
                        &latitude  = polyhedral_grid.latitude,
                        &areas     = polyhedral_grid.areas;
    longitude.resize( num_faces_final, 0. );
    latitude.resize( num_faces_final, 0. );
    areas.resize( num_faces_final, 0. );

    std::vector< std::vector<size_t> > &adjacency_indices = polyhedral_grid.adjacency_indices;
    adjacency_indices.resize( num_faces_final );

    std::vector< std::vector<double> > &ddlon_weights = polyhedral_grid.adjacency_ddlon_weights,
                                       &ddlat_weights = polyhedral_grid.adjacency_ddlat_weights,
                                       &d2dlon2_weights = polyhedral_grid.adjacency_d2dlon2_weights,
                                       &d2dlat2_weights = polyhedral_grid.adjacency_d2dlat2_weights;
    ddlon_weights.resize( num_faces_final );
    ddlat_weights.resize( num_faces_final );
    d2dlon2_weights.resize( num_faces_final );
    d2dlat2_weights.resize( num_faces_final );

    for ( size_t Iface = 0; Iface < num_faces_final; Iface++ ) {
        adjacency_indices[Iface].resize( size_neighbourhood );

        ddlon_weights[Iface].resize( size_neighbourhood, 0. );
        ddlat_weights[Iface].resize( size_neighbourhood, 0. );
        d2dlon2_weights[Iface].resize( size_neighbourhood, 0. );
        d2dlat2_weights[Iface].resize( size_neighbourhood, 0. );
    }

    size_t Iface;
    #pragma omp parallel \
    default(none) \
    shared( vertices, faces, longitude, latitude, areas, adjacency_indices ) \
    firstprivate( num_faces_final, size_neighbourhood ) \
    private( Iface )
    {
        #pragma omp for collapse(1) schedule(static)
        for ( Iface = 0; Iface < num_faces_final; Iface++ ) {

            // First, get the coordinates, the convert to spherical
            std::tuple< double, double, double > Cartesian_coordinates = midpoint_on_sphere(
                    vertices[std::get<0>(faces[Iface])], 
                    vertices[std::get<1>(faces[Iface])],
                    vertices[std::get<2>(faces[Iface])]
                    );
            longitude[Iface] = atan2( std::get<1>(Cartesian_coordinates), std::get<0>(Cartesian_coordinates) );
            latitude[Iface]  = (M_PI/2) - acos( std::get<2>(Cartesian_coordinates) );


            // Next, compute the area of the triangular face
            //  the vertices of the face are in Cartesian coordinates, so we need to convert
            //  the first tuple is longitudes, the second is latitudes
            areas[Iface] = SphericalTriangleArea(
                    std::tuple< double, double, double > ( 
                        atan2( std::get<1>( vertices[ std::get<0>(faces[Iface]) ] ), 
                               std::get<0>( vertices[ std::get<0>(faces[Iface]) ] ) ),
                        atan2( std::get<1>( vertices[ std::get<1>(faces[Iface]) ] ), 
                               std::get<0>( vertices[ std::get<1>(faces[Iface]) ] ) ),
                        atan2( std::get<1>( vertices[ std::get<2>(faces[Iface]) ] ), 
                               std::get<0>( vertices[ std::get<2>(faces[Iface]) ] ) )
                        ),
                    std::tuple< double, double, double > (
                        (M_PI/2) - acos( std::get<2>(vertices[ std::get<0>(faces[Iface])] ) ),
                        (M_PI/2) - acos( std::get<2>(vertices[ std::get<1>(faces[Iface])] ) ),
                        (M_PI/2) - acos( std::get<2>(vertices[ std::get<2>(faces[Iface])] ) )
                        )
                    );


            // Next, we build the adjacency matrix

            size_t num_found = 0;
            adjacency_indices[Iface][size_neighbourhood-1] = Iface;
            for ( size_t Jface = 0; Jface < num_faces_final; Jface++ ) {
                // We need to find the neighbours, which unfortunately means
                // checking through all faces to see if they share vertices

                if ( Iface == Jface ) { continue; }

                // Compare the vertices to see if Iface vertices are Jface vertices
                bool has_vertex_zero = (vertices[ std::get<0>(faces[Iface])] == vertices[std::get<0>(faces[Jface])])
                                    or (vertices[ std::get<0>(faces[Iface])] == vertices[std::get<1>(faces[Jface])])
                                    or (vertices[ std::get<0>(faces[Iface])] == vertices[std::get<2>(faces[Jface])]);

                bool has_vertex_one =  (vertices[ std::get<1>(faces[Iface])] == vertices[std::get<0>(faces[Jface])])
                                    or (vertices[ std::get<1>(faces[Iface])] == vertices[std::get<1>(faces[Jface])])
                                    or (vertices[ std::get<1>(faces[Iface])] == vertices[std::get<2>(faces[Jface])]);

                bool has_vertex_two =  (vertices[ std::get<2>(faces[Iface])] == vertices[std::get<0>(faces[Jface])])
                                    or (vertices[ std::get<2>(faces[Iface])] == vertices[std::get<1>(faces[Jface])])
                                    or (vertices[ std::get<2>(faces[Iface])] == vertices[std::get<2>(faces[Jface])]);

                // Check if two of the vertices are equal
                if ( has_vertex_zero + has_vertex_one + has_vertex_two == 2 ) {
                    adjacency_indices[Iface][num_found] = Jface;
                    num_found++;
                }

                // If we found the three neighbours, stop
                if ( num_found == 3 ) { break; }
            }
        }
    }

    // For various reason, we're going to re-arrange to have the south-most point
    // at the beginning
    double min_lat = M_PI/2;
    size_t Iface_min_lat = 0;
    for ( Iface = 0; Iface < num_faces_final; Iface++ ) {

        double local_lat = latitude[Iface];
        if ( local_lat < min_lat ) {
            min_lat = local_lat;
            Iface_min_lat = Iface;
        }

        // If we're close enough to the pole, stop
        if (min_lat < (-89. * M_PI / 180.)  ) {
            break;
        }
    }

    #if DEBUG >= 2
    fprintf(stdout, "Swapping indices 0 and %'zu to have souther-most at the beginning\n", Iface_min_lat);
    #endif

    // Now do the swapping
    std::iter_swap( areas.begin(), areas.begin() + Iface_min_lat );
    std::iter_swap( longitude.begin(), longitude.begin() + Iface_min_lat );
    std::iter_swap( latitude.begin(), latitude.begin() + Iface_min_lat );
    std::iter_swap( adjacency_indices.begin(), adjacency_indices.begin() + Iface_min_lat );

    // Also update their neighbours to have the correct 'address'
    adjacency_indices[0][size_neighbourhood-1] = 0;
    adjacency_indices[Iface_min_lat][size_neighbourhood-1] = Iface_min_lat;
    for ( size_t In = 0; In < 3; In++ ) {
        size_t Iface = adjacency_indices[0][In];
        for ( size_t Jn = 0; Jn < 3; Jn++ ) {
            if ( adjacency_indices[Iface][Jn] == Iface_min_lat ) {
                adjacency_indices[Iface][Jn] = 0;
            }
        }

        Iface = adjacency_indices[Iface_min_lat][In];
        for ( size_t Jn = 0; Jn < 3; Jn++ ) {
            if ( adjacency_indices[Iface][Jn] == 0 ) {
                adjacency_indices[Iface][Jn] = Iface_min_lat;
            }
        }
            
    }

    #if DEBUG >= 2
    fprintf( stdout, "BuildPolyhedralGrid: lat[0] = %f\n", latitude[0] );
    #endif



    // Now that we've set our nearest neighbours, add our neighbours neighbours 
    // until we've got a full adjacency set
    //  This needs to happen after all of the initial neighbours have been set.
    //  [hence a separate loop block]
    #pragma omp parallel \
    default(none) \
    shared( longitude, latitude, adjacency_indices, ddlon_weights, ddlat_weights, \
            d2dlon2_weights, d2dlat2_weights ) \
    firstprivate( num_faces_final, size_neighbourhood ) \
    private( Iface )
    {
        #pragma omp for collapse(1) schedule(static)
        for ( Iface = 0; Iface < num_faces_final; Iface++ ) {
            size_t num_found = 3;
            while ( num_found < size_neighbourhood ) {
                for ( size_t In = 0; In < num_found; In++ ) {
                    size_t Jface = adjacency_indices[Iface][In];

                    // For each neighbour, check if its neighbours
                    // are already one of our neighbours. If not, add it
                    for ( size_t Jn = 0; Jn < 3; Jn++ ) {
                        bool is_new = true;

                        // If our neighbour's neighbour is us, skip
                        if ( adjacency_indices[Jface][Jn] == Iface ) { continue; }

                        // Otherwise, check if it's already one of our neighbours
                        for ( size_t Kn = 0; Kn < num_found; Kn++ ) {
                            // Left: new face index.
                            // Right: index of existing neighbour
                            if ( adjacency_indices[Jface][Jn] == adjacency_indices[Iface][Kn] ) {
                                is_new = false;
                                break;
                            }
                        }
                        if ( is_new ) {
                            adjacency_indices.at(Iface).at(num_found) = adjacency_indices.at(Jface).at(Jn);
                            //adjacency_indices[Iface][num_found] = adjacency_indices[Jface][Jn];
                            num_found++;
                        }
                        if ( num_found == size_neighbourhood ) { break; }
                    }

                    if ( num_found == size_neighbourhood ) { break; }
                }
                if ( num_found == size_neighbourhood ) { break; }
            }


            //
            //// And, at last, compute the differentiation weights
            //

            Eigen::MatrixXd A(size_neighbourhood,6);
            double centre_lat = latitude[Iface];
            if ( std::fabs(centre_lat) < 0.*M_PI/180. ) {
                // straight lat/lon derivatives
                for ( size_t In = 0; In < size_neighbourhood; In++ ) {

                    double local_lon = longitude[ adjacency_indices[Iface][In] ] - longitude[Iface];
                    if (local_lon > M_PI) { local_lon = local_lon - 2*M_PI; }
                    else if (local_lon < -M_PI) { local_lon = local_lon + 2*M_PI; }

                    double local_lat = latitude[ adjacency_indices[Iface][In] ] - latitude[Iface];

                    // We construct the local quadratic approximation
                    A(In,0) = 1.;
                    A(In,1) = local_lon;
                    A(In,2) = local_lat;
                    A(In,3) = local_lon * local_lat;
                    A(In,4) = 0.5 * pow( local_lon, 2. );
                    A(In,5) = 0.5 * pow( local_lat, 2. );
                }
                Eigen::MatrixXd pinv = A.completeOrthogonalDecomposition().pseudoInverse();

                for ( size_t In = 0; In < size_neighbourhood; In++ ) {
                    ddlon_weights[  Iface][In] = pinv(1,In);
                    ddlat_weights[  Iface][In] = pinv(2,In);
                    d2dlon2_weights[Iface][In] = pinv(4,In);
                    d2dlat2_weights[Iface][In] = pinv(5,In);
                }
            } else {
                // use a local Cartesian projection to get derivative weights near the poles
                for ( size_t In = 0; In < size_neighbourhood; In++ ) {

                    double local_lon = longitude[ adjacency_indices[Iface][In] ] - longitude[Iface];
                    if (local_lon > M_PI) { local_lon = local_lon - 2*M_PI; }
                    else if (local_lon < -M_PI) { local_lon = local_lon + 2*M_PI; }

                    double local_lat = latitude[ adjacency_indices[Iface][In] ];

                    double x = cos( local_lat ) * sin( local_lon );
                    double y =   cos( centre_lat ) * sin( local_lat ) 
                               - sin( centre_lat ) * cos( local_lat ) * cos( local_lon );

                    // We construct the local quadratic approximation
                    A(In,0) = 1.;
                    A(In,1) = x;
                    A(In,2) = y;
                    A(In,3) = x * y;
                    A(In,4) = 0.5 * pow( x, 2. );
                    A(In,5) = 0.5 * pow( y, 2. );
                }
                Eigen::MatrixXd pinv = A.completeOrthogonalDecomposition().pseudoInverse();

                for ( size_t In = 0; In < size_neighbourhood; In++ ) {
                    ddlon_weights[  Iface][In] = cos(centre_lat) * pinv(1,In);
                    ddlat_weights[  Iface][In] = pinv(2,In);
                    d2dlon2_weights[Iface][In] = pow( cos(centre_lat), 2. ) * pinv(4,In)
                        + cos(centre_lat)*sin(centre_lat)*pinv(2,In);
                    d2dlat2_weights[Iface][In] = pinv(5,In);
                }
            }
        }

    }

}
