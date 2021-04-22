[TOC]
\page tutorials2

This tutorial shows the basic functionality of the coarse-graining codebase.

It includes:
- a python script to generate a sample dataset: `generate_data.py`
- a python script to interpret the coarse-graining results: `process_results.py`
- this tutorial outline

## About the sample data
 
 - `generate_data.py` script creates a file (`velocity_sample.nc`) that contains x- and y- velocity components
 - doubly-periodic Cartesian grid
 - domain is 400km in x by 150km in y
 - velocity field is a 'wavy jet' with noise

## What to do

1. Compile `Case_Files/coarse_grain.x` (see notes below) and copy into this directory.
2. Create the sample dataset (`python generate_data.py`)
3. Run coarse_grain.x (`./coarse_grain.x --input_file velocity_sample.nc --filter_scales "1e3 15e3 50e3 100e3"`). 
4. Run the python plotting routine (`python process_results.py`)

When specifying filtering scales, only pick a few. Other tutorials will cover how to handle outputs with many filter scales. A suggestion is: `1e3, 15e3, 50e3, 100e3`.
While the codebase is heavily parallelized, it unnecessary for this example and will be covered in later tutorials, and so we simply run in serial here.

### Notes when compiling

Make sure that the variables in `constants.hpp` are set appropraitely. These include:
- `CARTESIAN = true`
- `PERIODIC_X = true`
- `PERIODIC_Y = true`
- `COMP_BC_TRANSFERS = false`, this requires pressure and density, which are not included in this tutorial
- Since there are no 'land' areas in this tutorial, the choice of land treatment (set by `DEFORM_AROUND_LAND` and `FILTER_OVER_LAND`) has no effect.
- The grid is uniform, and so `UNIFORM_LON_GRID` and `UNIFORM_LAT_GRID` should be set accordingly. However, these are strictly optimization flags, and d not impact the output.
- `FULL_LON_SPAN = true`
- `MINIMAL_OUTPUT` and `NO_FULL_OUTPUTS` both to `false`. We want the full output files (file size is not a concern, since the sample is small)
- `CAST_TO_SINGLE` and `CAST_TO_INT` simply modify the data type (i.e. precision) used to store the outputs, and can be set however you wish
- `DO_TIMING = true` prints some summary information about how long different parts of the code took. Turn on or off as you wish.
- `APPLY_POSTPROCESS = false`, this will be handled in another tutorial
- `DO_OKUBOWEISS_ANALYSIS = false`, this will be handled in another tutorial
