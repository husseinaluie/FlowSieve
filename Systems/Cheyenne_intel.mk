# the following works on NCAR's cheyenne (tested on Oct 25, 2022):
# $ module reset
# $ module swap netcdf netcdf-mpi
# then
# $ make clean
# $ make Case_Files/coarse_grain.x

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

# Specify optimization flags for ALGLIB
ALGLIB_OPT_FLAGS:=-O3 -DAE_CPU=AE_INTEL

# Modules are automatically on lib dir
NETCDF_LIBS:=${NCAR_LIBS_NETCDF}
NETCDF_INCS:=${NCAR_INCS_NETCDF}

LIB_DIRS:=${NETCDF_LIBS}
INC_DIRS:=${NETCDF_INCS}

