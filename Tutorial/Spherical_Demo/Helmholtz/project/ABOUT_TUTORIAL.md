[TOC]
\pagetutorials3

This tutorial shows how to use the Helmholtz projection code.

Additionally, it uses a coarsen-refine approach that is very helpful for higher resolution settings.
Multiple coarsenings can also be used, which can substantially speed up convergence on a high resolution grid.

For instance, a 1/12 degree grid can be well projected in approximately 24 hours by using coarsened grids of 1 degree and 1/3 degree resolution.
In constrast, several days / a few weeks of directly projecting the 1/12 degree grid is still likely to not have sufficienly captured the largest scales.

The general idea
1. coarsen the velocity field to a lower resolution
2. compute the Helmholtz decomposition on the coarse data
3. interpolate the coarse result onto the fine grid
4. compute the Helmholtz decomposition on the fine data, using the results from the coarse data as a seed or 'guess' value

This tutorial includes:
- sample submit script for a SLURM-scheduled computing cluster: `submit_all_Helmholtz_steps.sh`
- a python script to produce some plots of the results
- this tutorial outline

## What to do

Note: After step 1, all of the subsequent steps are included in `submit_all_Helmholtz_steps.sh`

1. Compile `Case_Files/Helmholtz_projection.x`, `Case_Files/coarsen_grid_linear.x`, `Case_Files/refine_Helmholtz_seed.x` (see notes below) and copy into this directory.
2. Create the sample dataset
3. Create a coarsened grid version of the sample data
4. Run the Helmholtz projection on the coarsen data
5. Refine the coarse projection results to make a seed for the full resolution
6. Run the Helmholtz projection on the full resolution data
7. Run the python analysis script to make some plots of the projection


### Parallelization

For this particular example, it is not really worth parallelizing (it is a small enough problem), but notes on how you could parallelize are included below.

#### OpenMPI

OpenMPI parallelization is only used to divide the time and depth dimensions (since communication overhead in lat/lon would be costly).
As this example only has one point in time and depth, we cannot use more than one MPI process.

### Openmp

The Helmholtz decomposition routine has minimal OpenMP optimization (the free version of ALGLIB does not support it).
The other `for` loops (such as computing velocities, setting up the least squares problem, etc) are parellelized with OpenMP, but in general, extra OpenMP threads has pretty low return for the Helmholtz projection scripts.

### Notes when compiling

When compiling, use the `constants.hpp` provided in this directory.
