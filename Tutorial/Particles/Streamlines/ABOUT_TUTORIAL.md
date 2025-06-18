# Tutorial for Working Generating Streamlines {#tutorialsParticles1}
[TOC]

This tutorial walks through generating streamlines for a global velocity field

---

In this directory there is a python script `generate_data.py` that creates a sample velocity field.

The steps for this tutorial are
1. `python generate_data.py`
2. Using the `constants.hpp` file provide, compile `particles.x` (move to the main FlowSieve directory, copy this `constants.hpp` file there, and then call `make Case_File/particles.x`. Copy that `particles.x` here.
3. Run the script. You can use `./run_particles.sh` to run directly, or use the supplied SLURM submit script. [Note that you may need to adjust the SLURM flags according to your system.]
