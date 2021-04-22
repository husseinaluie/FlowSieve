[TOC]
\page install1

# Installing

## Initial (one-time) set-up

You will need to create a `system.mk` file in the main code directory that specifies the compilers and provides the compiler flags needed to link to the appropriate libraries.

Some system files are already provided in the `Systems` directory for a few computing environments.
Where appropraite, they specify which modules are needed to compile and run the code (note: these modules need to be loaded both when you compile the code and also when you run the code).

If you create a system file for a computing environment that was not previously included, please consider submitting it to the repository so that others can benefit from it.

### System File

The `Systems` directory contains a few sample system files. To prepare a system file for your machine, there are a few steps.
1. Copy the most appropraite system file from `Systems` to `system.mk` into the main directory.
2. Compilers: Indicate your chosen compilers (CXX for C++ and MPICXX for C++ with mpi)
3. Links: these may need to change depending on your compiler (e.g. `-fopenmp` to `-qopenmp` for intel compilers)
4. EXTRA_OPT_FLAGS: this is for other optimizations that you may not want on by default. For example, `-ip -ipo` for intel compilers.
5. Depending on how your libraries are set up, you may need to add various library and include directories to `LIB_DIRS` and `INC_DIRS`

# Compiling

A makefile is included to simplify compiling. 
As a first test of your install, try `make clean; make Case_Files/coarse_grain.x`.
If the executable compiles successfully, then a good next step is to check out the [tutorials that are provided](./TUTORIALS) to familiarize yourself with the code usage.

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
* `make tests`
 * Compiles the test routines in `Tests/`

## ALGLIB

Some of the executables (such as interpolator.x, toroidal_projection.x, and potential_projection.x) require the use of a third-party library, ALGLIB.
To compile this library, simply call `make ALGLIB` (should only need to be done once, unless you cal `make hardclean`).
