# Coarse Graining

This repository stores source code for running coarse graining procedures on netcdf files.

## Makefile

### make arguments
* `make clean`
  * This removes all object files (\*.o) in the source tree
* `make hardclean`
  * In addition to the `clean` removals, also removes the executables and symbol table (dSYM)
* `make all`
  * Makes coarse_grain.x
* `make coarse_grain.x`
  * Makes coarse_grain.x

### DEBUG flag

Setting the debug flag in the Makefile specifies how much information is printed
during runtime. In particular:

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
* DEBUG >= 2
  * Prints when finished computing cell areas (compute_areas.cpp)
  * Prints when entering the main filtering loop sequence (filtering.cpp)
  * Prints when the output file is initialized (initialize_output_file.cpp)
  * Prints when the output is written (write_to_output.cpp)
  * Prints when the output vorticity is written (write_vorticity.cpp)
* DEBUG >= 3
  * **[THIS WILL PRINT A LOT]**
  * Prints *every* index conversion (Index.cpp)
* DEBUG >= 4
  * **[THIS WILL PRINT A LOT]**
  * Prints the progress through Latitude and Longitude during filtering (filtering.cpp) 
  * Prints whenever a velocity conversion happens (Spher->Cart and Cart->Spher) 
* DEBUG >= 5
  * **[THIS WILL PRINT A LOT]**
  * Prints every kernel calculation (kernel.cpp)
  * Prints every distance calculation (distance.cpp)
