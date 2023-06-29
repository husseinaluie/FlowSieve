# Barkla is the hpc cluster at the University of Liverpool, UK 
# The following modules need to be loaded to compile FlowSieve on Barkla
#
# libs/netcdf_mpi/4.7.4/gcc-9.3.0+openmpi-3.1.6+hdf5_mpi-1.10.5
# apps/hdf5_mpi/1.10.5/gcc-9.3.0+openmpi-3.1.6
# mpi/openmpi/3.1.6/gcc-9.3.0
# libs/gcc/9.3.0

# Specify compilers
CXX    ?= g++
MPICXX ?= mpicxx

# Linking flags for netcdf
LINKS:=-lnetcdf -lhdf5_hl -lhdf5 -lm -ldl -lz -fopenmp

# Default compiler flags
CFLAGS:=-Wall -std=c++14

# Debug flags
DEBUG_FLAGS:=-g
DEBUG_LDFLAGS:=-g

# Basic optimization flags
OPT_FLAGS:=-O3

# Extra optimization flags (intel inter-process optimizations)
EXTRA_OPT_FLAGS:=

# Specify optimization flags for ALGLIB
ALGLIB_OPT_FLAGS:=-O3

# Modules are automatically on lib dir
NETCDF_LIBS="-L/opt/gridware/depots/e2b91392/el7/pkg/libs/netcdf_mpi/4.7.4/gcc-9.3.0+openmpi-3.1.6+hdf5_mpi-1.10.5/lib"
NETCDF_INCS=`nc-config --cflags`

HDF5_LIBS="-L/opt/gridware/depots/e2b91392/el7/pkg/apps/hdf5_mpi/1.10.5/gcc-9.3.0+openmpi-3.1.6/lib"
HDF5_INCS="-I/opt/gridware/depots/e2b91392/el7/pkg/apps/hdf5_mpi/1.10.5/gcc-9.3.0+openmpi-3.1.6/include"

LIB_DIRS:=${NETCDF_LIBS} ${HDF5_LIBS}
INC_DIRS:=${NETCDF_INCS} ${HDF5_INCS}
