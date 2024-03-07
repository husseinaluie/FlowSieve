##
## This system file is for TACC Stampede3
##	using GCC compilers
##
# The modules used for this system file were:
# module load gcc/13.2.0 impi/21.11 pnetcdf/4.9.2 phdf5/1.14.3

# Specify compilers
CXX    ?= g++
MPICXX ?= mpicxx

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

# Link in netcdf
GCC_LIBS=-L${TACC_GCC_LIB}
GCC_INCS=-I${TACC_GCC_INC}
NETCDF_LIBS=-L${TACC_PNETCDF_LIB}
NETCDF_INCS=-I${TACC_PNETCDF_INC}
HDF_LIBS=-L${TACC_PHDF5_LIB}
HDF_INCS=-I${TACC_PHDF5_INC}

LIB_DIRS:=${NETCDF_LIBS} ${HDF_LIBS} ${GCC_LIBS}
INC_DIRS:=${NETCDF_INCS} ${HDF_INCS} ${GCC_INCS}
