# The following modules were last used
#
#  Modules: 
#  module load comp/intel/2021.2.0   
#  module load mpi/impi/2021.2.0   
#  module load netcdf4/4.8.0
#
#  Note: load python modules *before* loading
#        those listed above, otherwise netcdf
#        libraries may have conflicts


# Specify compilers
CXX     ?= icc
MPICXX  ?= mpicxx

# Linking flags for netcdf
LINKS:=-lnetcdf -lhdf5_hl -lhdf5 -lz -lcurl -qopenmp -lm

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
NETCDF_LIBS="-L/usr/local/other/netcdf4/4.8.0/impi/2021.2.0/lib -L/usr/local/other/hdf5/1.12.0/intel-2021.2.0_impi-2021.2.0//lib -L/usr/local/intel/oneapi/2021/mpi/2021.2.0/lib"
NETCDF_INCS=`nc-config --cflags`

LIB_DIRS:=${NETCDF_LIBS}
INC_DIRS:=${NETCDF_INCS}
