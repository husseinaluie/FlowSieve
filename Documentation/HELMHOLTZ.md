# Helmholtz Decomposition {#helmholtz1}
[TOC]

# Helmholtz Routines {#helmholtz1-0}

---

The main executables for apply a Helmholtz decomposition are `toroidal_projection` and `potential_projection`, which respectively compute the streamfunction (giving the divergence-free velocity) and potential (giving the irrotational velocity) fields.

This is done using a sparse least-squares solver using the ALGLIB library.

The main variables controlling the convergence are the maximum number of iterations and the target relative tolerance.
If the maximum number of iterations are reach, or if the error goes below the tolerance, the solver halts and outputs the computed terms.

## Decomposing High Resolution Velocities {#helmholtz1-1}

In the case of high-resolution velocities, two logistical concerns arise.
1. Memory: the least-squares systems scales like N^2, so if you double the number of points in each horizontal dimension, the problem becomes 16 times larger.
2. Small-scale bias. Since the projection solves a Laplace problem, it is implicitely weighted to favour converging small scales before large scales. This means that on a high-resolution grid, it can take a very long time for the largest scales to converge.

A very practical work-around for this problem is to do coarsen/refine approach.
For example, suppose you have global data at 1/20 degree resolution.
Helmholtz-decomposition that grid directly would be very inefficient and computationally expensive.
Instead, coarsen the data to multiple lower-resolution grids: 1 degree, 1/4 degree, 1/12 degree, for instance.
Begin by first decomposing the 1 degree data, the result of which can be used as a seed to speed up the 1/4 degree decomposition.
The 1/4 can then be fed into the 1/12 degree, and so on.

These steps can be summarized as:
1. Coarsen the to a grid that can be solved efficienty (e.g. 1 degree)
2. Apply Helmholtz decomposition.
3. Using the coarse result, create a 'seed' file for a high-resolution Helmholtz decomposition
4. Repeat steps 2 and 3 until the final resolution is achieved.

Two auxiliary executables are provided for this purpose.
* `coarsen_grid` takes in velocity data and produces another data file on a coarse lat/lon grid (user specifies the coarsening factor as a command-line input)
* `refine_Helmholtz_seed` takes in the Helmholtz outputs from one grid and interpolates (linear interpolation) onto a finer grid. The result is then output to a file that can be read in by the main Helmholtz decomposition routines.
