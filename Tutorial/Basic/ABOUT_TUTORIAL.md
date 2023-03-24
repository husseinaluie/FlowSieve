# Basic Tutorial {#tutorials2}

\brief This tutorial shows the basic functionality of the coarse-graining codebase.

# Tutorial
[TOC]

---

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

When compiling, use the `constants.hpp` file provided in this directory.
