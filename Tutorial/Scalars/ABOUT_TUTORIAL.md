# Small Tutorial for Working on a Sphere {#tutorialsScalar}
[TOC]

This tutorial walks through analysing a scalar field with coarse_grain_scalars.x

---

The generated data here is low resolution (1-degree grid spacing) and fairly boring, but it runs in ~1 minute on a single processor, and illustrates the general idea.

1. First, compile `coarse_grain_scalars.x` and copy it into this directory
2. Run the `run_all_steps.sh` script (`./run_all_steps.sh`)
3. Afterwards, you can run the jupyter notebook to analyse the results.


### Notes when compiling

When compiling the executables, set the `constants.hpp` variables as follows:
- `CARTESIAN = false`
- `PERIODIC_X = true`
- `PERIODIC_Y = false`
- `DEFORM_AROUND_LAND = false`
- `FILTER_OVER_LAND = true`
- `UNIFORM_LON_GRID = true`
- `UNIFORM_LAT_GRID = true`
- `EXTEND_DOMAIN_TO_POLES = false`
- `FULL_LON_SPAN = true`
- `MINIMAL_OUTPUT = false`
- `NO_FULL_OUTPUTS = false`
- `CAST_TO_SINGLE = true`
- `CAST_TO_INT = false`
- `APPLY_POSTPROCESS = false`
