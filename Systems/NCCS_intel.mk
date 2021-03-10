# The following modules were last used
#
#  Modules: 
#  module load comp/intel/19.1.3.304
#  module load mpi/impi/19.1.3.304
#  module load netcdf4/4.7.4

# Specify compilers
CXX     ?= icc
MPICXX  ?= mpicxx

# Linking flags for netcdf
LINKS:=-lnetcdf -lhdf5_hl -lhdf5 -lz -lcurl -qopenmp -lm

# Default compiler flags
CFLAGS:=-Wall -std=c++14

# Debug flags
DEBUG_FLAGS:=-g
DEBUG_LDFLAGS:=-g

# Basic optimization flags
OPT_FLAGS:=-O3 -fp-model fast=2

# Extra optimization flags (intel inter-process optimizations)
EXTRA_OPT_FLAGS:=-ip -ipo

# Modules are automatically on lib dir
NETCDF_LIBS=`nc-config --cxx4libs`
NETCDF_INCS=`nc-config --cxx4flags`

LIB_DIRS:=${NETCDF_LIBS}
INC_DIRS:=${NETCDF_INCS}
