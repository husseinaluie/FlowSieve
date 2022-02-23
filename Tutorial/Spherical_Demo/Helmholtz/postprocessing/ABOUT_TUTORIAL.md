[TOC]
\pagetutorials3

This tutorial shows the postprocessing functionality of the coarse-graining routine.

It includes:
- a python script to generate a region definitions file: `define_geographic_regions.py`
- a python script to interpret the coarse-graining results: `plot_results.py`
- a sample submit script for a SLURM-scheduled computing cluster: `submit_filter_for_postprocessing.sh`
- a bash and python script to merge the output files: merge_postprocess_files.sh and  merge_postprocess_results.py
- this tutorial outline

## What to do

1. Compile `Case_Files/coarse_grain_helmholtz.x` (see notes below) and copy into this directory.
3. Create the regions definitions file (`python define_geographic_regions.py`)
4. Run coarse_grain_helmholtz.x ( see the sample submit script for details )
5. Merge the individual postprocess files into one file that contains all scales ( `./merge_postprocess_files.sh` )
6. Run the python plotting routine (`python plot_results.py`)

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
- Since there are no 'land' areas in this tutorial, the choice of land treatment (set by `DEFORM_AROUND_LAND` and `FILTER_OVER_LAND`) has no effect.
- The grid is uniform, and so `UNIFORM_LON_GRID` and `UNIFORM_LAT_GRID` should be set accordingly. However, these are strictly optimization flags, and d not impact the output.
- `FILTER_OVER_LAND = true`
- `EXTEND_DOMAIN_TO_POLES = false`
- `FULL_LON_SPAN = true`
- `MINIMAL_OUTPUT` and `NO_FULL_OUTPUTS` both to `true`. We focus on the post-process outputs in this tutorial.
- `CAST_TO_SINGLE` and `CAST_TO_INT` simply modify the data type (i.e. precision) used to store the outputs, and can be set however you wish
- `APPLY_POSTPROCESS = true`, this will generate the outputs used in this tutorial
- `POSTPROCESS_DO_ZONAL_MEANS = true`, also produce zonal means in addition to area-means
- `POSTPROCESS_DO_TIME_MEANS = false`, do not produce time-means (only one time-point)
- `DO_OKUBOWEISS_ANALYSIS = false`, this will be handled in another tutorial
