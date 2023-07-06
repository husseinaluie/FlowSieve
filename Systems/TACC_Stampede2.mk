# The following modules were last used
#
#  Modules: 
#  module load intel/19.1.1 
#  module load impi/19.0.9 
#  module load phdf5/1.10.4
#  module load parallel-netcdf/4.6.2
#
#  Note: load python/anaconda modules *before* 
#  		 loading those listed above, otherwise 
#  		 netcdf libraries may have conflicts

# Specify compilers
CXX     = icpc
MPICXX  = mpiicpc

# Linking flags for netcdf
LINKS:=-lnetcdf_c++4 -lnetcdf -lhdf5_hl -lhdf5 -lz -qopenmp -lm

# Default compiler flags
CFLAGS:=-Wall -std=c++14

# Debug flags
DEBUG_FLAGS:=-g3 -ggdb
DEBUG_LDFLAGS:=-g3 -ggdb

# Basic optimization flags
OPT_FLAGS:=-O3 -fp-model fast=2

# Extra optimization flags (intel inter-process optimizations)
EXTRA_OPT_FLAGS:=-ip -ipo

# Specify optimization flags for ALGLIB
ALGLIB_OPT_FLAGS:=-O3 -DAE_CPU=AE_INTEL

# Modules are automatically on lib dir
NETCDF_LIBS=-L${TACC_NETCDF_LIB}
NETCDF_INCS=-I${TACC_NETCDF_INC}

HDF_LIBS=-L${TACC_HDF5_LIB}
HDF_INCS=-I${TACC_HDF5_INC}

LIB_DIRS:=${NETCDF_LIBS} ${HDF_LIBS}
INC_DIRS:=${NETCDF_INCS} ${HDF_INCS}
