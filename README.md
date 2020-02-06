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

### Postprocessing

While the core executable gives you the option to run post-processing routines online, you do also have the option to run them after the fact.
These scripts can be found in the _Postprocess_ directory. 

* Postprocesss/integrator.cpp
 * compile via _make Postprocess/integrator.x_
 * run via _integrator.x <file-to-process.nc> <file-to-create.nc>_ 

**!!THIS IS VERY OUTDATED!! Use on-line routines if possible until this is updated.**

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
 * In addition to the `clean` removals, also removes the executables, symbol table (dSYM), and documentation
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
 * Use _DEBUG = 0_ for normal production runs
 * Use _DEBUG = 1_ if you want to keep track on the progress of a longer production run
 * Use _DEBUG = 2_ if you're running into some issues and want to narrow it down a bit
 * Going beyond this is really only necessary / useful if you're running into some fatal errors that you can't pinpoint
 * Setting _DEBUG_ to be negative is generally not advised. Setting to 0 shouldn't produce overly much output, and certainly not enough to hamper performance. If you're trying to silence errors, make sure you understand _why_ the errors are happening, and that you're really okay with ignoring them.

In particular:

* DEBUG <= -2
  * This setting **silences all netcdf errors**
* DEBUG <= -1
  * This setting gives **no** standard output (printing to screen).
* DEBUG >= 0
  * This is the default setting
  * Prints the date/time of compilation
  * Prints a  note when the main script finishes.
  * Prints the progress through Time and Depth during filtering (filtering.cpp)
* DEBUG >= 1
  * Prints when starting to read inputs.
  * Prints when starting to compute cell areas
  * Prints when starting to convert spherical velocities to Cartesian (filtering.cpp)
  * Prints when creating the output file (filtering.cpp)
  * Prints a 10-dot sequence showing progress through Lat/Lon points
* DEBUG >= 2
  * Prints when finished computing cell areas (compute_areas.cpp)
  * Prints when entering the main filtering loop sequence (filtering.cpp)
  * Prints when the output file is initialized (initialize_output_file.cpp)
  * Prints when the output is written (write_to_output.cpp)
  * Prints when the output vorticity is written (write_vorticity.cpp)
* DEBUG >= 3
  * **[THIS WILL PRINT A LOT]**
  * Prints the progress through Latitude and Longitude during filtering (filtering.cpp) 
* DEBUG >= 4
  * **[THIS WILL PRINT A LOT]**
  * Prints linenumber checks in the main filtering loop to show progress.
* DEBUG >= 5
  * **[THIS WILL PRINT A LOT]**
  * Prints *every* index conversion (Index.cpp)
* DEBUG >= 6
  * **[THIS WILL PRINT A LOT]**
  * Prints every kernel calculation (kernel.cpp)
  * Prints every distance calculation (distance.cpp)
