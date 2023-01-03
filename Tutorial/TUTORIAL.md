# Tutorials {#tutorials1}

# FlowSieve Tutorials
[TOC]

---

## Basic

Introduces the basic workings on a sample dataset.
Nothing fancy, but introduces the basics of running the code and plotting some of the results.

Tutorial [details are here](\ref tutorials2).

## Scalars

Introduces the basic functionality of `coarse_grain_scalars` by filtering a sample density field.

Tutorial [details are here](\ref tutorialsScalar)

## Spherical Demo

These tutorial works through the various steps involved in working with spherical data.

* The tutorial [details are here](\ref tutorials4).
 * This is a very small case that can run in ~5 minutes on one processor.
* The tutorial [details are here](\ref tutorials3).
 * There are three components, that go through Helmholtz decomposition, filtering, and using on-line post-processing tools.


---

## Conda Environment

For convenience and reproducibility, an Anaconda environment file is included in the main Tutorial directory.
To create the environment, run `conda env create -f environment.yml` (this may take several minutes, but should only need to be done once).
Barring errors, this will create an environment named `FlowSieve-tutorial-env`.
To activate the environment, simply call `conda activate FlowSieve-tutorial-env`.

Alternatively, if you don't use Anaconda to manage python packages, the environment file lists the package dependencies needed to run the tutorial.
The specific package version are sufficient, but may not be strictly necessary (i.e. other version combinations may work, but have not been tested).
