[TOC]
# Coarse Graining

This repository stores source code for running coarse graining procedures on netcdf files.

---

## Compilation / Installation

1. Copy the appropriate file from `./Systems/` to `./system.mk`
 * The system directory contains files that declare system compiler information.
 * You may need to create a new file to correspond to your system
2. Call `make all` in the root directory to build the main executable.
3. Call `make tests` to compile the unit test routines.
 * The resulting executables are then stored in `Tests/`

### Interpolator

1. The interpolator requires the ALGLIB package. 
 * Compile via `make ALGLIB`
 * only needs to be done once (unless you call `make hardclean`)
2. Call `make Case_Files/interpolator.x`

### Helmholtz Decomposition

[Go to this page](./HELMHOLTZ.md)

### Postprocessing

Post-processing (such as region-averaging, Okubo-Weiss histogram binning, time-averaging, etc) can be enabled and run on-line
by setting the **APPLY_POSTPROCESS** flag in constants.hpp to **true**.

This will produce an additional output file for each filtering scale.

### System File

The `Systems` directory contains a few sample system files. To prepare a system file for your machine, there are a few steps.
1. Copy the most appropraite system file from `Systems` to `system.mk` into the main directory.
2. Compilers: Indicate your chosen compilers (CXX for C++ and MPICXX for C++ with openmpi)
3. Links: these may need to change depending on your compiler (e.g. `-fopenmp` to `-qopenmp` for intel compilers)
4. EXTRA_OPT_FLAGS: this is for other optimizations that you may not want on by default. For example, `-ip -ipo` for intel compilers.
5. Depending on how your libraries are set up, you may need to add various library and include directories to `LIB_DIRS` and `INC_DIRS`

---

## Command-line Arguments

1. `--version`
 * Calling `./coarse_grain.x --version` prints a summary of the constants / variables used when compiling

---

## Function Map

See the function map for [filtering] to get an overview of the function dependencies.
[filtering]: @ref filtering() "the main filtering function"

---

## Makefile

### make arguments
* `make clean`
 * This removes all object files (\*.o) in the source tree
* `make hardclean`
 * In addition to the `clean` removals, also removes the executables, symbol table (dSYM)
* `make Case_Files/<filename>.x`
 * Makes the specified executable. `Case_Files` contains the main files that can be compiled into executables.
* `make docs`
 * Makes the doxygen-produced documentation. **Note: `doxygen` and `dot` must be installed and on the path**
* `make cleandocs`
 * Removes the previous documentation build
* `make tests`
 * Compiles the test routines in `Tests/`

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
