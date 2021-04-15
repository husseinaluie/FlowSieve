[TOC]

This tutorial shows the basic functionality of the coarse-graining codebase.

It includes:
- a python script to generate a sample dataset: `generate_data.py`
- a python script to generate a region definitions file: `define_geographic_regions.py`
- a python script to interpret the coarse-graining results: `process_results.py`
- a sample submit script for a SLURM-scheduled computing cluster: `sample_submit.sh`
- this tutorial outline

## About the sample data
 
 - `generate_data.py` script creates a file (`velocity_sample.nc`) that contains zonal- and meridional- velocity components
 - full spherical domain, 0.5 degree resolution, points at poles removed
 - velocity field is a collection of eddies of varying sizes, randomly placed throughout the globe

## What to do

1. Compile `Case_Files/coarse_grain.x` (see notes below) and copy into this directory.
2. Create the sample dataset (`python generate_data.py`)
3. Create the regions definitions file (`python define_geographic_regions.py`)
4. Run coarse_grain.x ( see the sample submit script for details )
5. Run the python plotting routine (`python process_results.py`)

When specifying filtering scales, consider a wide sweep. It can also be beneficial to use logarithmically-spaced scales, for plotting purposes.
Python can be helpful for this. For example, `numpy.logspace( np.log(50e3), np.log(2000e3), 10 )` would produce 10 logarithmically-spaced
filter scales between 50km and 2000km.

Hint: to print filter scales to only three significant digits, the `numpy.format_float_scientific` function can help.
> scales = numpy.logspace( np.log(50e3), np.log(2000e3), 10 )
> [print( numpy.format_float_scientific( scale, precision = 2 ) ) for scale in scales]
The result of which is: 5.e+04 7.53e+04 1.13e+05 1.71e+05 2.58e+05 3.88e+05 5.85e+05 8.81e+05 1.33e+06 2.00e+06

### Parallelization

While our dataset is not particularly high resolution (only 358 x 720 grid points), we are considering a large sweep of filter scales,
so it might be nice to speed things up with a bit of parallelization.

#### OpenMPI

OpenMPI parallelization is only used to divide the time and depth dimensions (since communication overhead in lat/lon would be costly).
As this example only has one point in time and depth, we cannot use more than one MPI process.

### Openmp

Openmp (threading) is used to parallelize everything within an MPI process (i.e. lat/lon, as well as the allotment of time/depth points to the MPI process).
As a result, Openmp parallelizations are _always_ available for coarse-graining.

To specify the number of processors, simply set the OMP_NUM_THREADS environment variable. (see the sample submit script for details)

### Notes when compiling

Make sure that the variables in `constants.hpp` are set appropraitely. These include:
- `CARTESIAN = false`
- `PERIODIC_X = true`
- `PERIODIC_Y = false`
- `COMP_BC_TRANSFERS = false`, this requires pressure and density, which are not included in this tutorial
- Since there are no 'land' areas in this tutorial, the choice of land treatment (set by `DEFORM_AROUND_LAND` and `FILTER_OVER_LAND`) has no effect.
- The grid is uniform, and so `UNIFORM_LON_GRID` and `UNIFORM_LAT_GRID` should be set accordingly. However, these are strictly optimization flags, and d not impact the output.
- `FULL_LON_SPAN = true`
- `MINIMAL_OUTPUT` and `NO_FULL_OUTPUTS` both to `true`. We focus on the post-process outputs in this tutorial.
- `CAST_TO_SINGLE` and `CAST_TO_INT` simply modify the data type (i.e. precision) used to store the outputs, and can be set however you wish
- `DO_TIMING = true` prints some summary information about how long different parts of the code took. Turn on or off as you wish.
- `APPLY_POSTPROCESS = true`, this will generate the outputs used in this tutorial
- `DO_OKUBOWEISS_ANALYSIS = false`, this will be handled in another tutorial
