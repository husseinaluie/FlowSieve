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
5. Merge the individual postprocess files into one file that contains all scales ( see bash scripts for usage details )
6. Run the python plotting routine (`python plot_results.py` or juse the Jupyter notebook)

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

When compiling, use the `constants.hpp` provided in this directory.
