# Installing {#install1}

\brief A guide to installing and compiling FlowSieve.

# Installation
[TOC]

---

## Requirements

* netcdf4 (tested on versions 4.7 and 4.8)
  * requires parallel installation
  * hdf5-mpi required
* mpi (openmpi or intel mpi, not mpich)
* compilers: either gcc or icpc (intel compiler)
  * For gcc, tested on gcc8 and gcc9. When tested, FlowSieve compiled _without_ warnings or errors. 
    * Note that FlowSieve currently **does not support for gcc10+**.
  * If Intel compilers are available, they are recommended, since the additional compiler optimizations provide a minor performance improvement. 
  * Intended to be compiled with c++14 standard (this is specified in the compiler flags in the `system.mk` file)
* curl

The following, while not required to run FlowSieve itself, are required to run the tutorial. These requirements are outline in `Tutorials/environment.yml` to facility create an anaconda environment. Note that these _exact_ package versions may not be required, but they are sufficient. That is, other package versions may work, but have not been tested.
* python=3.7
* numpy=1.17
* matplotlib=3.1
* netCDF4=1.4
* scipy=1.3
* pip

---

## Initial (one-time) set-up

First, you will need to copy the FlowSieve source files onto your HPC system.
The source files are available through git (https://github.com/husseinaluie/FlowSieve), allowing you to directly copy/fork the repository, or download a zipped version.

### System File

Next, you will need to create a `system.mk` file in the main code directory that specifies the compilers and provides the compiler flags needed to link to the appropriate libraries.

Some system files are already provided in the `Systems` directory for a few computing environments, and may be a helpful reference.
Where appropraite, they specify which modules are needed to compile and run the code (note: these modules need to be loaded both when you compile the code and also when you run the code).

If you create a system file for a computing environment that was not previously included, please consider submitting it to the repository so that others can benefit from it.

The `Systems` directory contains a few sample system files. To prepare a system file for your machine, there are a few steps.
1. Copy the most appropraite system file from `Systems` to `system.mk` into the main directory.
2. Compilers: Indicate your chosen compilers (CXX for C++ and MPICXX for C++ with mpi)
3. Links: these may need to change depending on your compiler (e.g. `-fopenmp` to `-qopenmp` for intel compilers)
4. EXTRA_OPT_FLAGS: this is for other optimizations that you may not want on by default. For example, `-ip -ipo` for intel compilers.
5. Depending on how your libraries are set up, you may need to add various library and include directories to `LIB_DIRS` and `INC_DIRS`

If you are unsure of how to make an appropriate `system.mk` file, please feel free to contact the FlowSieve developers for support.

# Compiling

A makefile is included to simplify compiling. 
As a first test of your install, try `make clean; make Case_Files/coarse_grain.x`.
If the executable compiles successfully, then a good next step is to check out the [tutorials that are provided](\ref tutorials1) to familiarize yourself with the code usage.

## make arguments
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

## ALGLIB

Some of the executables (such as interpolator.x and helmholtz_projection.x) require the use of a third-party library, [ALGLIB](https://www.alglib.net).
To compile this library, simply call `make ALGLIB` (should only need to be done once, unless you call `make hardclean`, change system libraries, etc.).
