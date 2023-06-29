# Tutorial for Filtering Scalars {#tutorialsScalar}
[TOC]

This tutorial walks through analysing a scalar field with coarse_grain_scalars.x

---

The generated data here is low resolution (1-degree grid spacing) and fairly boring, but it runs in ~1 minute on a single processor, and illustrates the general idea.

1. First, compile `coarse_grain_scalars.x` and copy it into this directory
2. Run the `run_all_steps.sh` script (`./run_all_steps.sh`)
3. Afterwards, you can run the jupyter notebook to analyse the results.


### Notes when compiling

When compiling the executables, set the `constants.hpp` with the provided file.
