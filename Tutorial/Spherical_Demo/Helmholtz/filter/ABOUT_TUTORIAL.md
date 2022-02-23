[TOC]
\pagetutorials3

This part of the tutorial shows how to filter the Helmholtz-decomposed flow.

This tutorial includes:
- sample submit script for a SLURM-scheduled computing cluster: `submit_filtering.sh`
- a python script to produce some plots of the results
- this tutorial outline

## What to do

1. Compile `Case_Files/coarse_grain_helmholtz.x` (see notes below) and copy into this directory.
2. Run the coarse-graining routine, following the command structure outlined in the submit script.
3. Run the python analysis script to make some plots of the projection.


### Parallelization

For this particular example, it is not really worth parallelizing (it is a small enough problem), but notes on how you could parallelize are included below.

#### OpenMPI

OpenMPI parallelization is only used to divide the time and depth dimensions (since communication overhead in lat/lon would be costly).
As this example only has one point in time and depth, we cannot use more than one MPI process.

### Openmp

Unlike the Helmholtz decomposition code, the coarse-graining code takes full advantage of OpenMP parallelization in the lat-lon dimensions.
While this is a coarse grid and does not require high computing power, we can take advantage of a full compute node-worth of processors, if we desired.

### Notes when compiling

Make sure that the variables in `constants.hpp` are set appropraitely. These include:
- `CARTESIAN = false`
- `PERIODIC_X = true`
- `PERIODIC_Y = false`
- The grid is uniform, and so `UNIFORM_LON_GRID` and `UNIFORM_LAT_GRID` should be set accordingly. However, these are strictly optimization flags, and d not impact the output.
- `FILTER_OVER_LAND = true`
- `EXTEND_DOMAIN_TO_POLES = false`
- `FULL_LON_SPAN = true`
- `CAST_TO_SINGLE` and `CAST_TO_INT` simply modify the data type (i.e. precision) used to store the outputs, and can be set however you wish
- `MINIMAL_OUTPUT = false` to get extra output variables
- `APPLY_POSTPROCESS = false` (this will be explored in the postprocessing part of the tutorial)
- `DO_OKUBOWEISS_ANALYSIS = false`, this will be handled in another tutorial
