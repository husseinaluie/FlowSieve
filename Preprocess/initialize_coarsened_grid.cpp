#include "../constants.hpp"
#include "../functions.hpp"
#include "../preprocess.hpp"
#include <algorithm>
#include <vector>
#include <omp.h>
#include <math.h>

void initialize_coarsened_grid(
        dataset & coarsened_grid,
        const dataset & reference_grid,
        const int Nlat_coarse,
        const int Nlon_coarse
        ) {

    // Reset the coarsened grid
    coarsened_grid.clear();

    // The coarsened grid is assumed to be a MeshGrid type
    coarsened_grid.full_Ntime  = reference_grid.full_Ntime;
    coarsened_grid.full_Ndepth = reference_grid.full_Ndepth;

    coarsened_grid.Ntime  = reference_grid.Ntime;
    coarsened_grid.Ndepth = reference_grid.Ndepth;

    coarsened_grid.time = reference_grid.time;
    coarsened_grid.depth = reference_grid.depth;

    if ( constants::GRID_TYPE != constants::GridType::MeshGrid ) {
        // If we're not on a mesh grid, then build a polyhedral grid
        // and stop. That grid contains everything we need
        size_t Npts = Nlat_coarse * (size_t) Nlon_coarse;
        BuildPolyhedralGrid( coarsened_grid, Npts );
        Npts = coarsened_grid.longitude.size();
        coarsened_grid.mask.resize( Npts, true );

        // Setting myCounts and myStarts
        coarsened_grid.myCounts.resize(3);
        coarsened_grid.myCounts = { reference_grid.myCounts[0], reference_grid.myCounts[1], Npts };

        coarsened_grid.myStarts.resize(3);
        coarsened_grid.myStarts = { reference_grid.myStarts[0], reference_grid.myStarts[1], 0 };

        return;
    }
    coarsened_grid.Nlat = Nlat_coarse;
    coarsened_grid.Nlon = Nlon_coarse;


    if ( constants::GRID_TYPE == constants::GridType::MeshGrid ) {
        coarsened_grid.myCounts.resize(4);
        coarsened_grid.myCounts = { reference_grid.myCounts[0], reference_grid.myCounts[1], Nlat_coarse, Nlon_coarse };

        coarsened_grid.myStarts.resize(4);
        coarsened_grid.myStarts = { reference_grid.myStarts[0], reference_grid.myStarts[1], 0, 0 };
    } else {
        coarsened_grid.myCounts.resize(3);
        coarsened_grid.myCounts = { reference_grid.myCounts[0], reference_grid.myCounts[1], Nlat_coarse * Nlon_coarse };

        coarsened_grid.myStarts.resize(3);
        coarsened_grid.myStarts = { reference_grid.myStarts[0], reference_grid.myStarts[1], 0 };
    }


    coarsened_grid.mask.resize(Nlon_coarse * (size_t)Nlat_coarse, 1);

    const double dlat = M_PI / Nlat_coarse,
                 dlon = 2 * M_PI / Nlon_coarse;
    const size_t Npts = Nlat_coarse * (size_t) Nlon_coarse;
    size_t index;
    if ( constants::GRID_TYPE == constants::GridType::MeshGrid ) {
        coarsened_grid.latitude.resize(Nlat_coarse, 0.);
        coarsened_grid.longitude.resize(Nlon_coarse, 0.);
        for ( int Ilat = 0; Ilat < Nlat_coarse; Ilat++ ) {
            coarsened_grid.latitude[Ilat] = (Ilat+0.5)*dlat - (M_PI/2.);
        }
        for ( int Ilon = 0; Ilon < Nlon_coarse; Ilon++ ) {
            coarsened_grid.longitude[Ilon] = Ilon*dlon - M_PI;
        }
        coarsened_grid.compute_cell_areas();
    } else {
        coarsened_grid.latitude.resize(Npts, 0.);
        coarsened_grid.longitude.resize(Npts, 0.);
        coarsened_grid.areas.resize(Npts, 0.);
        for ( int Ilat = 0; Ilat < Nlat_coarse; Ilat++ ) {
            for ( int Ilon = 0; Ilon < Nlon_coarse; Ilon++ ) {
                index = Ilat * Nlon_coarse + Ilon;
                coarsened_grid.longitude[index] = Ilon*dlon - M_PI;
                coarsened_grid.latitude[index] = (Ilat+0.5)*dlat - (M_PI/2.);

                coarsened_grid.areas[index] = pow(6371e3, 2.) * dlat * dlon * cos(coarsened_grid.latitude[index]);
            }
        }
    }


    if ( constants::GRID_TYPE == constants::GridType::MeshGrid ) { return; }
    // For mesh grids, we're done.
    //
    // For unstructured grids, build the adjacency / derivative stuff
    //  'unstructured' for the original grid, but the coarsened are all
    //  uniform lat/lon because they may as well be
    coarsened_grid.adjacency_indices.resize( Npts );

    coarsened_grid.adjacency_ddlon_weights.resize( Npts );
    coarsened_grid.adjacency_ddlat_weights.resize( Npts );

    coarsened_grid.adjacency_d2dlon2_weights.resize( Npts );
    coarsened_grid.adjacency_d2dlat2_weights.resize( Npts );

    for ( int Ilat = 0; Ilat < Nlat_coarse; Ilat++ ) {
        for ( int Ilon = 0; Ilon < Nlon_coarse; Ilon++ ) {
            index = Ilat * Nlon_coarse + Ilon;

            coarsened_grid.adjacency_indices[index].resize( coarsened_grid.num_neighbours+1, 0 );

            coarsened_grid.adjacency_ddlon_weights[index].resize( coarsened_grid.num_neighbours+1, 0 );
            coarsened_grid.adjacency_ddlat_weights[index].resize( coarsened_grid.num_neighbours+1, 0 );

            coarsened_grid.adjacency_d2dlon2_weights[index].resize( coarsened_grid.num_neighbours+1, 0 );
            coarsened_grid.adjacency_d2dlat2_weights[index].resize( coarsened_grid.num_neighbours+1, 0 );

            coarsened_grid.adjacency_indices[index][0] = 
                (Ilat + 0) * Nlon_coarse + ( (Ilon - 1 + Nlon_coarse) % Nlon_coarse );
            coarsened_grid.adjacency_indices[index][8] = 
                (Ilat + 0) * Nlon_coarse + ( (Ilon - 0) % Nlon_coarse );
            coarsened_grid.adjacency_indices[index][1] = 
                (Ilat + 0) * Nlon_coarse + ( (Ilon + 1) % Nlon_coarse );

            // First lon derivative
            coarsened_grid.adjacency_ddlon_weights[index][0] = -1. / (2.*dlon);
            coarsened_grid.adjacency_ddlon_weights[index][1] =  1. / (2.*dlon);
            /*
            // If we put this back, need to correct weighting at poles!
            coarsened_grid.adjacency_ddlon_weights[index][0] = -1. / (2.*dlon) * (1./2);
            coarsened_grid.adjacency_ddlon_weights[index][1] =  1. / (2.*dlon) * (1./2);

            coarsened_grid.adjacency_ddlon_weights[index][2] = -1. / (2.*dlon) * (1./4);
            coarsened_grid.adjacency_ddlon_weights[index][4] =  1. / (2.*dlon) * (1./4);

            coarsened_grid.adjacency_ddlon_weights[index][5] = -1. / (2.*dlon) * (1./4);
            coarsened_grid.adjacency_ddlon_weights[index][7] =  1. / (2.*dlon) * (1./4);
            */

            // Second lon derivative
            coarsened_grid.adjacency_d2dlon2_weights[index][0] =  1. / pow(dlon,2.);
            coarsened_grid.adjacency_d2dlon2_weights[index][8] = -2. / pow(dlon,2.);
            coarsened_grid.adjacency_d2dlon2_weights[index][1] =  1. / pow(dlon,2.);
            /*
            // If we put this back, need to correct weighting at poles!
            coarsened_grid.adjacency_d2dlon2_weights[index][0] =  1. / pow(dlon,2.) * (1./2);
            coarsened_grid.adjacency_d2dlon2_weights[index][8] = -2. / pow(dlon,2.) * (1./2);
            coarsened_grid.adjacency_d2dlon2_weights[index][1] =  1. / pow(dlon,2.) * (1./2);

            coarsened_grid.adjacency_d2dlon2_weights[index][2] =  1. / pow(dlon,2.) * (1./4);
            coarsened_grid.adjacency_d2dlon2_weights[index][3] = -2. / pow(dlon,2.) * (1./4);
            coarsened_grid.adjacency_d2dlon2_weights[index][4] =  1. / pow(dlon,2.) * (1./4);

            coarsened_grid.adjacency_d2dlon2_weights[index][5] =  1. / pow(dlon,2.) * (1./4);
            coarsened_grid.adjacency_d2dlon2_weights[index][6] = -2. / pow(dlon,2.) * (1./4);
            coarsened_grid.adjacency_d2dlon2_weights[index][7] =  1. / pow(dlon,2.) * (1./4);
            */

            if ( Ilat == 0 ) {
                coarsened_grid.adjacency_ddlat_weights[index][8] = -1.5 / dlat;// * (1./2);
                coarsened_grid.adjacency_ddlat_weights[index][6] =  2.  / dlat;// * (1./2);
                coarsened_grid.adjacency_ddlat_weights[index][3] = -0.5 / dlat;// * (1./2);

                coarsened_grid.adjacency_d2dlat2_weights[index][8] =  1. / pow(dlat,2.);// * (1./2);
                coarsened_grid.adjacency_d2dlat2_weights[index][6] = -2. / pow(dlat,2.);// * (1./2);
                coarsened_grid.adjacency_d2dlat2_weights[index][3] =  1. / pow(dlat,2.);// * (1./2);
            } else if ( Ilat == Nlat_coarse - 1 ) {
                coarsened_grid.adjacency_ddlat_weights[index][8] =  1.5 / dlat;// * (1./2);
                coarsened_grid.adjacency_ddlat_weights[index][3] = -2.  / dlat;// * (1./2);
                coarsened_grid.adjacency_ddlat_weights[index][6] =  0.5 / dlat;// * (1./2);

                coarsened_grid.adjacency_d2dlat2_weights[index][8] =  1. / pow(dlat,2.);// * (1./2);
                coarsened_grid.adjacency_d2dlat2_weights[index][3] = -2. / pow(dlat,2.);// * (1./2);
                coarsened_grid.adjacency_d2dlat2_weights[index][6] =  1. / pow(dlat,2.);// * (1./2);
            } else {
                coarsened_grid.adjacency_ddlat_weights[index][3] = -1. / (2.*dlat);// * (1./2);
                coarsened_grid.adjacency_ddlat_weights[index][6] =  1. / (2.*dlat);// * (1./2);

                coarsened_grid.adjacency_d2dlat2_weights[index][3] =  1. / pow(dlat,2.);// * (1./2);
                coarsened_grid.adjacency_d2dlat2_weights[index][8] = -2. / pow(dlat,2.);// * (1./2);
                coarsened_grid.adjacency_d2dlat2_weights[index][6] =  1. / pow(dlat,2.);// * (1./2);
            }
            /*
            coarsened_grid.adjacency_ddlat_weights[index][2] = (1./2) * coarsened_grid.adjacency_ddlat_weights[index][3];
            coarsened_grid.adjacency_ddlat_weights[index][5] = (1./2) * coarsened_grid.adjacency_ddlat_weights[index][6];

            coarsened_grid.adjacency_ddlat_weights[index][4] = (1./2) * coarsened_grid.adjacency_ddlat_weights[index][3];
            coarsened_grid.adjacency_ddlat_weights[index][7] = (1./2) * coarsened_grid.adjacency_ddlat_weights[index][6];


            coarsened_grid.adjacency_d2dlat2_weights[index][2] = (1./2) * coarsened_grid.adjacency_d2dlat2_weights[index][3];
            coarsened_grid.adjacency_d2dlat2_weights[index][0] = (1./2) * coarsened_grid.adjacency_d2dlat2_weights[index][8];
            coarsened_grid.adjacency_d2dlat2_weights[index][5] = (1./2) * coarsened_grid.adjacency_d2dlat2_weights[index][6];

            coarsened_grid.adjacency_d2dlat2_weights[index][4] = (1./2) * coarsened_grid.adjacency_d2dlat2_weights[index][3];
            coarsened_grid.adjacency_d2dlat2_weights[index][1] = (1./2) * coarsened_grid.adjacency_d2dlat2_weights[index][8];
            coarsened_grid.adjacency_d2dlat2_weights[index][7] = (1./2) * coarsened_grid.adjacency_d2dlat2_weights[index][6];
            */


            // And finally set the indices
            if ( Ilat > 0 ) {
                coarsened_grid.adjacency_indices[index][2] = 
                    (Ilat - 1) * Nlon_coarse + ( (Ilon - 1 + Nlon_coarse) % Nlon_coarse );
                coarsened_grid.adjacency_indices[index][3] = 
                    (Ilat - 1) * Nlon_coarse + ( (Ilon + 0) % Nlon_coarse );
                coarsened_grid.adjacency_indices[index][4] = 
                    (Ilat - 1) * Nlon_coarse + ( (Ilon + 1) % Nlon_coarse );
            } else {
                coarsened_grid.adjacency_indices[index][2] = 
                    (Ilat + 2) * Nlon_coarse + ( (Ilon - 1 + Nlon_coarse) % Nlon_coarse );
                coarsened_grid.adjacency_indices[index][3] = 
                    (Ilat + 2) * Nlon_coarse + ( (Ilon + 0) % Nlon_coarse );
                coarsened_grid.adjacency_indices[index][4] = 
                    (Ilat + 2) * Nlon_coarse + ( (Ilon + 1) % Nlon_coarse );
            }

            if ( Ilat < Nlat_coarse - 1 ) {
                coarsened_grid.adjacency_indices[index][5] = 
                    (Ilat + 1) * Nlon_coarse + ( (Ilon - 1 + Nlon_coarse) % Nlon_coarse );
                coarsened_grid.adjacency_indices[index][6] = 
                    (Ilat + 1) * Nlon_coarse + ( (Ilon + 0) % Nlon_coarse );
                coarsened_grid.adjacency_indices[index][7] = 
                    (Ilat + 1) * Nlon_coarse + ( (Ilon + 1) % Nlon_coarse );
            } else {
                coarsened_grid.adjacency_indices[index][5] = 
                    (Ilat - 2) * Nlon_coarse + ( (Ilon - 1 + Nlon_coarse) % Nlon_coarse );
                coarsened_grid.adjacency_indices[index][6] = 
                    (Ilat - 2) * Nlon_coarse + ( (Ilon + 0) % Nlon_coarse );
                coarsened_grid.adjacency_indices[index][7] = 
                    (Ilat - 2) * Nlon_coarse + ( (Ilon + 1) % Nlon_coarse );
            }

        }
    }

}
