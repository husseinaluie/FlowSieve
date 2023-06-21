# The following modules were last used
#
#  Modules: 
#
# module load cpu
# module load PrgEnv-gnu
# module load cray-mpich/8.1.21
# module load cray-hdf5-parallel/1.12.2.3
# module load cray-netcdf-hdf5parallel/4.9.0.3
#

# Specify compilers
CXX     ?= g++
MPICXX  ?= mpic++

# Linking flags for netcdf
LINKS:=-lnetcdf -lhdf5_hl -lhdf5 -lm -ldl -lz -fopenmp -lm

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
ALGLIB_OPT_FLAGS:=-O3

# Modules are automatically on lib dir
NETCDF_INCS=-I/opt/cray/pe/netcdf-hdf5parallel/4.9.0.3/gnu/9.1/include -I/opt/cray/pe/hdf5-parallel/1.12.2.3/gnu/9.1/include
NETCDF_LIBS=-L/opt/cray/pe/netcdf-hdf5parallel/4.9.0.3/gnu/9.1/lib -L/opt/cray/pe/hdf5-parallel/1.12.2.3/gnu/9.1/lib

LIB_DIRS:=${NETCDF_LIBS}
INC_DIRS:=${NETCDF_INCS}
