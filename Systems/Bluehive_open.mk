# The following modules were last used
#  openmpi/2.0.1/b1  
#  hdf5/1.8.19/b1     
#  netcdf/4.3.3.1
#  gcc/8.2.0/b1
#  fftw3/3.3.6/b1

# Specify compilers
CXX     ?= g++
MPICXX  ?= mpicxx

# Linking flags for netcdf
LINKS:=-lnetcdf -lhdf5_hl -lhdf5 -lz -lcurl -fopenmp

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
LIB_DIRS:=
INC_DIRS:=
