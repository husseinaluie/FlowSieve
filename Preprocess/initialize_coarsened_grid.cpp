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

    // Set time and depth variables
    coarsened_grid.full_Ntime  = reference_grid.full_Ntime;
    coarsened_grid.full_Ndepth = reference_grid.full_Ndepth;

    coarsened_grid.Ntime  = reference_grid.Ntime;
    coarsened_grid.Ndepth = reference_grid.Ndepth;

    coarsened_grid.time = reference_grid.time;
    coarsened_grid.depth = reference_grid.depth;

    if ( not(constants::GRID_TYPE == constants::GridType::MeshGrid) ) {
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

    // Set the lat/lon dimension sizes
    coarsened_grid.Nlat = Nlat_coarse;
    coarsened_grid.Nlon = Nlon_coarse;

    coarsened_grid.myCounts.resize(4);
    coarsened_grid.myCounts = { reference_grid.myCounts[0], reference_grid.myCounts[1], Nlat_coarse, Nlon_coarse };

    coarsened_grid.myStarts.resize(4);
    coarsened_grid.myStarts = { reference_grid.myStarts[0], reference_grid.myStarts[1], 0, 0 };

    // Make a new mask variable of appropriate size.
    const size_t Npts = Nlat_coarse * (size_t) Nlon_coarse;
    coarsened_grid.mask.resize( Npts, 1 );

    const double dlat = M_PI / (Nlat_coarse + 1),
                 dlon = 2 * M_PI / Nlon_coarse;

    // Compute the new (uniform) grid
    coarsened_grid.latitude.resize(Nlat_coarse, 0.);
    for ( int Ilat = 0; Ilat < Nlat_coarse; Ilat++ ) {
        coarsened_grid.latitude[Ilat] = (Ilat+1)*dlat - (M_PI/2.);
    }

    coarsened_grid.longitude.resize(Nlon_coarse, 0.);
    for ( int Ilon = 0; Ilon < Nlon_coarse; Ilon++ ) {
        coarsened_grid.longitude[Ilon] = Ilon*dlon - M_PI;
    }

    // Compute cell areas
    coarsened_grid.compute_cell_areas();

}
