# The following modules were last used
#
#  Modules: 
#  module load intel/2020.1
#  module load intel-mpi/2019.8.254
#  module load intel-mkl/2022.0
#  module load netcdf-c/4.7.4
#
#  Note: load python/anaconda modules *before* 
#  		 loading those listed above, otherwise 
#  		 netcdf libraries may have conflicts

# Specify compilers
CXX    ?= icpc
MPICXX  ?= mpicxx

# Linking flags for netcdf
LINKS:=-lnetcdf -lhdf5_hl -lhdf5 -lz -fopenmp -lm

# Default compiler flags
CFLAGS:=-Wall -std=c++14

# Debug flags
DEBUG_FLAGS:=-g3 -ggdb
DEBUG_LDFLAGS:=-g3 -ggdb

# Basic optimization flags
OPT_FLAGS:=-O3

# Extra optimization flags (intel inter-process optimizations)
EXTRA_OPT_FLAGS:=

# Specify optimization flags for ALGLIB
ALGLIB_OPT_FLAGS:=-O3 -DAE_CPU=AE_INTEL

# Modules are automatically on lib dir
NETCDF_LIBS="-L/nasa/netcdf/4.4.1.1_mpt/lib"
NETCDF_INCS=`nc-config --cflags`

LIB_DIRS:=${NETCDF_LIBS}
INC_DIRS:=${NETCDF_INCS}
