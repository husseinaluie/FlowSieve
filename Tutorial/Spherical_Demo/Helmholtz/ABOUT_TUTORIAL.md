[TOC]
\pagetutorials3

This tutorial shows how to use the Helmholtz projection code.

Additionally, it uses a coarsen-refine approach that is very helpful for higher resolution settings.
Multiple coarsenings can also be used, which can substantially speed up convergence on a high resolution grid.

The general idea
1. coarsen the velocity field to a lower resolution
2. compute the Helmholtz decomposition on the coarse data
3. interpolate the coarse result onto the fine grid
4. compute the Helmholtz decomposition on the fine data, using the results from the coarse data as a seed or 'guess' value

This tutorial includes:
- a python script to generate a sample dataset: `generate_data.py`
- sample submit script for a SLURM-scheduled computing cluster: `submit_all_Helmholtz_steps.sh`
- this tutorial outline

## About the sample data
 
 - `generate_data.py` script creates a file (`velocity_sample.nc`) that contains zonal- and meridional- velocity components
 - full spherical domain, 0.5 degree resolution, points at poles removed
 - velocity field is a collection of eddies of varying sizes, randomly placed throughout the globe

## What to do

Note: Steps 3-6 are all included in the sample submit script

1. Compile `Case_Files/toroidal_projection.x`, `Case_Files/potential_projection.x`, `Case_Files/coarsen_grid.x`, `Case_Files/refine_Helmholtz_seed.x` (see notes below) and copy into this directory.
2. Create the sample dataset (`python generate_data.py`)
3. Create a coarsened grid version
4. Run the Helmholtz projection on the coarsen data
5. Refine the coarse projection results to make a seed for the full resolution
6. Run the Helmholtz projection on the full resolution data


### Parallelization

For this particular example, it is not really worth parallelizing (it is a small enough problem), but notes on how you could parallelize are included below.

#### OpenMPI

OpenMPI parallelization is only used to divide the time and depth dimensions (since communication overhead in lat/lon would be costly).
As this example only has one point in time and depth, we cannot use more than one MPI process.

### Openmp

The Helmholtz decomposition routine has minimal OpenMP optimization (the free version of ALGLIB does not support it).
In general, extra OpenMP threads has pretty low return.

### Notes when compiling

Make sure that the variables in `constants.hpp` are set appropraitely. These include:
- `CARTESIAN = false`
- `PERIODIC_X = true`
- `PERIODIC_Y = false`
- The grid is uniform, and so `UNIFORM_LON_GRID` and `UNIFORM_LAT_GRID` should be set accordingly. However, these are strictly optimization flags, and d not impact the output.
- `FULL_LON_SPAN = true`
- `CAST_TO_SINGLE` and `CAST_TO_INT` simply modify the data type (i.e. precision) used to store the outputs, and can be set however you wish
