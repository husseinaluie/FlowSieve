# The following modules were last used
#  netcdf/4.7.0/b1 (loading this loads the rest)
#  hdf5/1.10.5/b2
#  curl/7.56.1
#  openmpi/2.1.6/b2
#  gcc/8.2.0/b1

# Specify compilers
CXX     = g++
MPICXX  = mpicxx

# Linking flags for netcdf
LINKS:=-lnetcdf -lhdf5_hl -lhdf5 -lm -ldl -lz -lcurl -fopenmp

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
LIB_DIRS:=-L/software/netcdf/4.7.0/b1/lib
INC_DIRS:=-I/software/netcdf/4.7.0/b1/include
