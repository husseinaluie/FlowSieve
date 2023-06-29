# Small Tutorial for Working on a Sphere {#tutorials5}
[TOC]

This tutorial highlights the MPI-functionality of FlowSieve by analysing a depth-varying velocity field.
For a non-MPI tutorial, please [see this tutorial](\ref tutorials4).

---

The generated data here is very low resolution (2-degree grid spacing), and is intended to be run on multiple processors.
Running on 24 processors, the runtime is ~10 minutes.

1. First, compile `Helmholtz_projection.x` and `coarse_grain_helmholtz.x` and copy them into this directory
2. If on a SLURM scheduler, use the provided submit script. Alternately, modify `run_all_steps.sh` to specify the number of MPI processors to use, and run the `run_all_steps.sh` script (`./run_all_steps.sh`)
3. Afterwards, you can run the jupyter notebook to analyse the results.

---

## Reducing the MPI-requirement

24 processors is a fairly heavy requirement if you are not running on a computing cluster.
You can simply run on fewer processors (highest efficiency if the number of processors divides evenly into 48 - the number of vertical levels), but at the cost of increasing the runtime.

To reduce the processor cost without increasing runtime, you can decrease the number of vertical levels proportionately. E.g. you can reduce the vertical levels to 12 in order to run on 6 processors in a similar amount of time.

To adjust the number of vertical levels, you can adjust line 13 of `generate_data_sphere.py`, which reads `Nlon, Nlat, Ndepth = int(360//2), int(180//2), 48`. The last number, `48`, specifies the number of vertical levels.

When running the code, you can use any number of MPI ranks up to the number of verticals levels, but the most efficient use of processors occurs when the number of MPI ranks divides evenly into the number of vertical levels.


---

There is a known issue where the python scripts that run after calling `coarse_grain_helmholtz.x` in the `run_all_steps.sh` script fail to run withour printing an error. 
If the coarse-graining successfully runs (i.e. directory `outputs` is created and populated with `postprocessing_*.nc` files) but the `RESULTS_*.nc` files are not created, then you will need to manually call the three python scripts listed at the end of `run_all_steps.sh` script.


### Notes when compiling

When compiling the executables, set the `constants.hpp` with the provided file.
