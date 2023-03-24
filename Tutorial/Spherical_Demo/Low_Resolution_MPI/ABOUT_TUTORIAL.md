# Small Tutorial for Working on a Sphere {#tutorials5}
[TOC]

This tutorial highlights the MPI-functionality of FlowSieve by analysing a depth-varying velocity field.
For a non-MPI tutorial, please [see this tutorial](\ref tutorials4).

---

The generated data here is very low resolution (2-degree grid spacing), and is intended to be run on multiple processors.
Running in 40 processors, the runtime is ~10 minutes.

1. First, compile `Helmholtz_projection.x` and `coarse_grain_helmholtz.x` and copy them into this directory
2. If on a SLURM scheduler, use the provided submit script. Alternately, modify `run_all_steps.sh` to specify the number of MPI processors to use, and run the `run_all_steps.sh` script (`./run_all_steps.sh`)
3. Afterwards, you can run the jupyter notebook to analyse the results.

---

There is a known issue where the python scripts that run after calling `coarse_grain_helmholtz.x` in the `run_all_steps.sh` script fail to run withour printing an error. 
If the coarse-graining successfully runs (i.e. directory `outputs` is created and populated with `postprocessing_*.nc` files) but the `RESULTS_*.nc` files are not created, then you will need to manually call the three python scripts listed at the end of `run_all_steps.sh` script.


### Notes when compiling

When compiling the executables, set the `constants.hpp` with the provided file.
