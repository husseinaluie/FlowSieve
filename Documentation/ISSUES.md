# Solutions to Common Problems {#issues1}

\brief A description of some know issues / bugs, with solutions where available.

# Known Issues
[TOC]

---

## I/O

### Fails to read in array / zero-dimensional variables (when should have more dimensions)

This appears to be a result of the netcdf format.
The FlowSieve codebase requires netcdf-4 format.

You can check the format of your input file by calling `ncdump -k <filename>`.
If the returned value is not "netCDF-4", then you can convert it by calling `nccopy -k netCDF-4 original.nc nc4_version.nc`



## ALGLIB-related

### alglib::ap_error
When running some of the alglib routines (interpolator, Helmholtz projections, etc), you may encounter the error message 
"terminate called after throwing an instance of 'alglib::ap_error'".
So far, this seems to come up for two reasons
1. not enough memory (a somewhat confusing form of memory error)
2. as the result of NaN values. 

#### Memory

Running the Helmholtz interpolator on large grids can be expensize, and there appears to be more
MPI memory overhead / inefficient splitting than expected [ hunting this down is a to-do ].
Currents options are: use more memory ( may require more OpenMP threads ), or manually split
the input files across time / depth and run the Helmholtz solver on the problems seperately.

Aplogies for the inconvenience.

#### NaNs

There are two main reasons that have caused this so far.

For the interpolator, this often means:
1. The input fields (potential temparture, salinity, sea level anomaly) have different mask,
which is not accounted for in the input file. (The code has been modified so as to mostly avoid this).
2. The input fields are using unexpected units or have erroneous values. Please check variable values.
(e.g. temperature should be in Celsius, salinity values should not be negative)

The output file pre_interp.nc (requires DEBUG >= 1) shows the values computed before interpolation.
Checking those values may help to clarify the issue.


For the Helmholtz projection, this typically means something funky happened in matrix creation, often involving the poles


### Bus error / bad file descriptor

The main time that I have encountered these errors has been a result of insufficient memory, particularly when creating the sparse matrices.
Running out of available memory during the sparse matrix creation seems to cause funky pointer issues that lead to bus errors.

On large grids, the Helmholtz projectors can be memory hogs (even when using sparse systems).

On HPC systems, you can often get more memory by requesting more threads ( sometimes called tasks per cpu ).
While this is inefficient since the solver is not threaded (at least, not the free version that we use), is can solve the memory problem.

If more memory is not a realistic option and you are close to having enough memory, reducing the order of the differentiation scheme will reduce the size of the differentiation stencil,
thereby reducing the number of non-zero points in the sparse matrices.

And, of course, using only one time/depth point per MPI rank will reduce memory used in storage.
