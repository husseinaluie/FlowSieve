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

When compiling, use the `constants.hpp` provided in this directory.
