[TOC]
# Coarse Graining

This repository stores source code for running coarse graining procedures on netcdf files.

---

## Tutorial

A series of [basic tutorials](./Tutorial/TUTORIAL.md) (\ref tutorials1) are provided to outline both various usage cases as well as how to use / process the outputs.


---

## Methods

Some details regarding underlying methods are discussed [on this page](./Documentation/METHODS.md) (\ref methods1) (warning, math content).

### Helmholtz Decomposition

For notes about the Helmholtz decomposition, [go to this page](./Documentation/HELMHOLTZ.md) (\ref helmholtz1).

---

## Compilation / Installation

For notes on installation, please see [this page](./Documentation/INSTALL.md) (\ref install1).

---

## Postprocessing

Post-processing (such as region-averaging, Okubo-Weiss histogram binning, time-averaging, etc) can be enabled and run on-line
by setting the **APPLY_POSTPROCESS** flag in constants.hpp to **true**.

This will produce an additional output file for each filtering scale.

Various geographic regions of interest can be provided in a netcdf file.

---

## Command-line Arguments

1. `--version`
 * Calling `./coarse_grain.x --version` prints a summary of the constants / variables used when compiling


### Specifying Filtering Scales

When specifying filtering scales, consider a wide sweep. It can also be beneficial to use logarithmically-spaced scales, for plotting purposes.
Python can be helpful for this. For example, `numpy.logspace( np.log(50e3), np.log(2000e3), 10 )` would produce 10 logarithmically-spaced
filter scales between 50km and 2000km.

Hint: to print filter scales to only three significant digits, the `numpy.format_float_scientific` function can help.
> import numpy
> number_of_scales = 10
> smallest_scale = 50e3
> largest_scale  = 2000e3
> scales = numpy.logspace( numpy.log(smallest_scale), numpy.log(largest_scale), number_of_scales )
> [print( numpy.format_float_scientific( scale, precision = 2 ) ) for scale in scales]

---

## Known Issues

Some known issues (with solutions where available) are [given on this page](./Documentation/ISSUES.md) (\ref issues1).

---

## Technical Matters

### DEBUG flag

Setting the debug flag in the Makefile specifies how much information is printed
during runtime. 

This list may not be quite up-to-date. Rule of thumb:
 * Use _DEBUG <= -2_ silences netcdf errors
 * Use _DEBUG = 0_ for normal production runs
 * Use _DEBUG = 1_ if you want to keep track on the progress of a longer production run
 * Use _DEBUG = 2_ if you're running into some issues and want to narrow it down a bit
 * Going beyond this is really only necessary / useful if you're running into some fatal errors that you can't pinpoint
 * Setting _DEBUG_ to be negative is generally not advised. Setting to 0 shouldn't produce overly much output, and certainly not enough to hamper performance. If you're trying to silence errors, make sure you understand _why_ the errors are happening, and that you're really okay with ignoring them.


### Function Map

See the function map for [filtering] to get an overview of the function dependencies.
[filtering]: @ref filtering() "the main filtering function"
