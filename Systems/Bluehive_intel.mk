# The following modules were last used
#  impi/2018
#  intel/2018 (automatically loaded with impi/2018)
#  netcdf/4.4.1.1/b1
#  hdf5/1.8.17/b4  (automatically loaded with netcdf/4.4.1.1/b1)
#  module unload impi/2017  (this is loaded with netcdf too, but we don't want it)

# Specify compilers
CXX     ?= icc
MPICXX  ?= mpiicc

# Linking flags for netcdf
LINKS:=-lnetcdf -lhdf5_hl -lhdf5 -lz -lcurl -qopenmp

# Default compiler flags
CFLAGS:=-Wall 

# Debug flags
DEBUG_FLAGS:=-g
DEBUG_LDFLAGS:=-g

# Basic optimization flags
OPT_FLAGS:=-O3 -fp-model fast=2

# Extra optimization flags (intel inter-process optimizations)
EXTRA_OPT_FLAGS:=-ip -ipo

# Modules are automatically on lib dir
LIB_DIRS:=
INC_DIRS:=
