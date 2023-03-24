# Small Tutorial for Working on a Sphere {#tutorials4}
[TOC]

This tutorial walks through analysing a velocity field using Helmholtz decomposition with coarse-graining.

---

The generated data here is very low resolution (2-degree grid spacing), but runs in ~5 minutes on a single processor.

1. First, compile `Helmholtz_projection.x` and `coarse_grain_helmholtz.x` and copy them into this directory
2. Run the `run_all_steps.sh` script (`./run_all_steps.sh`)
3. Afterwards, you can run the jupyter notebook to analyse the results.

---

There is a known issue where the python scripts that run after calling `coarse_grain_helmholtz.x` in the `run_all_steps.sh` script fail to run withour printing an error. 
If the coarse-graining successfully runs (i.e. directory `outputs` is created and populated with `postprocessing_*.nc` files) but the `RESULTS_*.nc` files are not created, then you will need to manually call the three python scripts listed at the end of `run_all_steps.sh` script.


### Notes when compiling

When compiling the executables, set the `constants.hpp` with the provided file.
